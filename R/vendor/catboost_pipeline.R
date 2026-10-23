library(pROC)
library(PRROC)
library(catboost)

# MAIN TRAINING FUNCTION #####

train_model <- function(data, param_config) {
  # Validate model type and target compatibility
  validate_model_config(data, param_config)
  
  # Setup and validation
  file_paths <- setup_file_paths(param_config$results_dir)
  
  if (results_exist(file_paths$results_file)) {
    return()
    # return(load_existing_results(file_paths))
  }
  
  # Data preparation
  data_splits <- prepare_data_splits(
    data, 
    param_config$target, 
    param_config$fold_col,
    test_fold = param_config$test_fold,
    valid_ratio = param_config$valid_ratio
  )
  if (is.null(data_splits)) return(NULL)

  # Model training
  model_result <- train_catboost_model(data_splits, param_config$setting, param_config$model_type)

  # Evaluation and feature analysis
  evaluation <- evaluate_model(model_result$model, data_splits, param_config$target, param_config$model_type)
  feature_importance <- extract_feature_importance(model_result$model)
  shap_analysis <- calculate_shap_values(model_result$model, data_splits, param_config$model_type)
  predictions <- create_predictions_dataframe(data_splits, evaluation, param_config$model_type)
  
  # Save results
  saveRDS(model_result, file_paths$results_file)
  saveRDS(feature_importance, file_paths$feature_imp_file)
  saveRDS(shap_analysis, file_paths$shap_file)
  saveRDS(predictions, file_paths$predictions_file)
  saveRDS(evaluation, file_paths$evaluation_file)
  catboost.save_model(model_result$model, file_paths$model_file)
}

# VALIDATION FUNCTIONS #####

validate_model_config <- function(data, param_config) {
  # Check required parameters
  required_params <- c("model_type", "target", "fold_col", "test_fold", "valid_ratio", "setting", "results_dir")
  missing_params <- setdiff(required_params, names(param_config))
  if (length(missing_params) > 0) {
    stop("Missing required parameters: ", paste(missing_params, collapse = ", "))
  }
  
  # Validate model_type
  valid_types <- c("binary", "multiclass", "regression")
  if (!param_config$model_type %in% valid_types) {
    stop("model_type must be one of: ", paste(valid_types, collapse = ", "))
  }
  
  # Check target column exists
  if (!param_config$target %in% names(data)) {
    stop("Target column '", param_config$target, "' not found in data")
  }
  
  # Validate target compatibility with model type
  target_values <- data[[param_config$target]]
  unique_values <- unique(target_values[!is.na(target_values)])
  n_unique <- length(unique_values)
  
  if (param_config$model_type == "binary") {
    if (n_unique != 2) {
      stop("Binary classification requires exactly 2 unique target values. Found ", n_unique, 
           " unique values: ", paste(sort(unique_values), collapse = ", "))
    }
    if (!all(unique_values %in% c(0, 1))) {
      cat("Warning: Binary target values are not 0/1. Found: ", paste(sort(unique_values), collapse = ", "), "\n")
    }
  } else if (param_config$model_type == "multiclass") {
    if (n_unique < 3) {
      stop("Multiclass classification requires at least 3 unique target values. Found ", n_unique, 
           " unique values: ", paste(sort(unique_values), collapse = ", "))
    }
    if (n_unique > 50) {
      cat("Warning: Large number of classes (", n_unique, "). Consider if this should be regression instead.\n")
    }
    if (!is.numeric(target_values) && !is.factor(target_values) && !is.character(target_values)) {
      stop("Multiclass target must be numeric, factor, or character")
    }
  } else if (param_config$model_type == "regression") {
    if (!is.numeric(target_values)) {
      stop("Regression target must be numeric. Found type: ", class(target_values)[1])
    }
    if (n_unique < 10) {
      cat("Warning: Few unique values (", n_unique, ") for regression. Consider if this should be classification instead.\n")
    }
  }
}

# FILE MANAGEMENT FUNCTIONS #####

setup_file_paths <- function(results_dir) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  list(
    results_file = file.path(results_dir, "results.rds"),
    feature_imp_file = file.path(results_dir, "feature_importance.rds"),
    shap_file = file.path(results_dir,"shap_values.rds"),
    predictions_file = file.path(results_dir, "predictions.rds"),
    model_file = file.path(results_dir, "model.cbm"),
    evaluation_file = file.path(results_dir, "evaluation.rds")
  )
}

results_exist <- function(results_file) {
  file.exists(results_file)
}

load_existing_results <- function(file_paths) {
  cat("Results already exist, loading from cache...\n")
  list(
    results = readRDS(file_paths$results_file),
    feature_importance = readRDS(file_paths$feature_imp_file),
    shap_values = readRDS(file_paths$shap_file),
    predictions = readRDS(file_paths$predictions_file)
  )
}

# DATA PREPARATION FUNCTIONS #####

create_catboost_pools <- function(training_data, validation_data, test_data, target_col, fold_col = NULL) {
  # Columns to exclude from features
  exclude_cols <- target_col
  if (!is.null(fold_col)) {
    exclude_cols <- c(exclude_cols, fold_col)
  }
  
  list(
    train_pool = catboost.load_pool(
      data = training_data[, !names(training_data) %in% exclude_cols, drop = FALSE],
      label = training_data[[target_col]]
    ),
    validation_pool = catboost.load_pool(
      data = validation_data[, !names(validation_data) %in% exclude_cols, drop = FALSE],
      label = validation_data[[target_col]]
    ),
    test_pool = catboost.load_pool(
      data = test_data[, !names(test_data) %in% exclude_cols, drop = FALSE],
      label = test_data[[target_col]]
    )
  )
}

prepare_data_splits <- function(data, target_col, fold_col, test_fold = 1, valid_ratio = 0.2) {
  # Input validation
  if (!target_col %in% names(data)) {
    warning(paste("Target column", target_col, "not found in data"))
    return(NULL)
  }
  
  if (!fold_col %in% names(data)) {
    warning(paste("Fold column", fold_col, "not found in data"))
    return(NULL)
  }
  
  if (nrow(data) == 0) {
    warning("Data is empty")
    return(NULL)
  }
  
  if (valid_ratio <= 0 || valid_ratio >= 1) {
    warning("valid_ratio must be between 0 and 1")
    return(NULL)
  }
  
  # Get unique folds
  unique_folds <- unique(data[[fold_col]])
  
  # Check if test_fold exists
  if (!test_fold %in% unique_folds) {
    warning(paste("Test fold", test_fold, "not found in data"))
    return(NULL)
  }
  
  # Check if we have enough folds
  if (length(unique_folds) < 2) {
    warning("Need at least 2 folds for train/validation/test split")
    return(NULL)
  }
  
  # Create test data
  test_data <- data[data[[fold_col]] == test_fold, ]
  
  # Get remaining data (non-test folds)
  remaining_data <- data[data[[fold_col]] != test_fold, ]
  remaining_folds <- unique_folds[unique_folds != test_fold]
  
  if (nrow(remaining_data) == 0) {
    warning("No data remaining after removing test fold")
    return(NULL)
  }
  
  # Split remaining data into training and validation based on valid_ratio
  set.seed(123)  # For reproducible splits
  remaining_data$temp_id <- 1:nrow(remaining_data)
  
  # Sample validation indices based on valid_ratio
  n_valid <- floor(nrow(remaining_data) * valid_ratio)
  valid_indices <- sample(remaining_data$temp_id, n_valid)
  
  # Create training and validation data
  validation_data <- remaining_data[remaining_data$temp_id %in% valid_indices, ]
  training_data <- remaining_data[!remaining_data$temp_id %in% valid_indices, ]
  
  # Remove temporary ID column
  training_data$temp_id <- NULL
  validation_data$temp_id <- NULL
  remaining_data$temp_id <- NULL
  
  # Check if all splits have data
  if (nrow(training_data) == 0 || nrow(validation_data) == 0 || nrow(test_data) == 0) {
    warning("One or more data splits are empty")
    return(NULL)
  }
  
  # Return data splits instead of pools
  list(
    training_data = training_data,
    validation_data = validation_data,
    test_data = test_data,
    target_col = target_col,
    fold_col = fold_col,
    fold_info = list(
      test_fold = test_fold,
      remaining_folds = remaining_folds,
      valid_ratio = valid_ratio,
      train_ratio = 1 - valid_ratio,
      total_folds = length(unique_folds),
      n_train = nrow(training_data),
      n_validation = nrow(validation_data),
      n_test = nrow(test_data)
    )
  )
}

# MODEL TRAINING FUNCTIONS ######

get_default_params <- function(model_type) {
  base_params <- list(
    iterations = 1000,
    task_type = "GPU",
    devices = "0",
    boosting_type = "Plain",
    gpu_ram_part = 0.95,
    od_type = "Iter",
    od_wait = 50,
    learning_rate = 0.03,
    depth = 6
  )
  
  if (model_type == "binary") {
    base_params$loss_function <- "Logloss"
    base_params$custom_metric <- c("Recall", "Precision", "F1", "Accuracy")
  } else if (model_type == "multiclass") {
    base_params$loss_function <- "MultiClass"
    base_params$custom_metric <- c("MultiClass", "TotalF1", "Accuracy")
  } else if (model_type == "regression") {
    base_params$loss_function <- "RMSE"
    base_params$custom_metric <- c("RMSE", "MAE", "R2")
  }
  
  return(base_params)
}

apply_class_weights <- function(params, training_data, target, setting, model_type) {
  if (model_type == "regression") {
    # No class weights for regression
    params$class_weights <- NULL
    cat("Setting:", setting, "- No class weights (regression)\n")
    return(params)
  }
  
  if (setting %in% c("tuning_with_weights", "default_with_weights")) {
    class_counts <- table(training_data[[target]])
    if (length(class_counts) > 1) {
      class_weights <- as.numeric(nrow(training_data) / (length(class_counts) * class_counts))
      params$class_weights <- class_weights
      cat("Setting:", setting, "- Class weights:", class_weights, "\n")
    }
  } else {
    params$class_weights <- NULL
    cat("Setting:", setting, "- No class weights\n")
  }
  params
}

perform_hyperparameter_tuning <- function(params, train_pool, validation_pool, validation_data, target, model_type) {
  cat("Performing hyperparameter tuning...\n")
  
  learning_rates <- c(0.01, 0.03, 0.05, 0.1)
  depths <- c(3, 6, 8)
  
  best_metric <- ifelse(model_type == "regression", Inf, 0)
  best_params <- params
  
  for (lr in learning_rates) {
    for (depth in depths) {
      tuning_params <- params
      tuning_params$learning_rate <- lr
      tuning_params$depth <- depth
      tuning_params$iterations <- 300  # Reduced iterations for tuning
      
      cv_model <- catboost.train(train_pool, params = tuning_params, test_pool = validation_pool)
      
      if (model_type == "regression") {
        cv_preds <- catboost.predict(cv_model, validation_pool)
        cv_metric <- sqrt(mean((cv_preds - validation_data[[target]])^2))  # RMSE
        cat("  LR:", lr, "Depth:", depth, "RMSE:", cv_metric, "\n")
        
        if (cv_metric < best_metric) {
          best_metric <- cv_metric
          best_params$learning_rate <- lr
          best_params$depth <- depth
        }
      } else if (model_type == "binary") {
        cv_preds <- catboost.predict(cv_model, validation_pool, prediction_type = "Probability")
        cv_roc_obj <- roc(validation_data[[target]], cv_preds)
        cv_metric <- auc(cv_roc_obj)
        cat("  LR:", lr, "Depth:", depth, "AUC:", cv_metric, "\n")
        
        if (cv_metric > best_metric) {
          best_metric <- cv_metric
          best_params$learning_rate <- lr
          best_params$depth <- depth
        }
      } else if (model_type == "multiclass") {
        cv_preds <- catboost.predict(cv_model, validation_pool, prediction_type = "Class")
        cv_metric <- mean(cv_preds == validation_data[[target]])  # Accuracy
        cat("  LR:", lr, "Depth:", depth, "Accuracy:", cv_metric, "\n")
        
        if (cv_metric > best_metric) {
          best_metric <- cv_metric
          best_params$learning_rate <- lr
          best_params$depth <- depth
        }
      }
    }
  }
  
  cat("Best hyperparameters - LR:", best_params$learning_rate, "Depth:", best_params$depth, "\n")
  best_params
}

train_catboost_model <- function(data_split_pool, setting, model_type) {
  # Create pools from data
  pools <- create_catboost_pools(
    data_split_pool$training_data, 
    data_split_pool$validation_data, 
    data_split_pool$test_data, 
    data_split_pool$target_col, 
    data_split_pool$fold_col
  )
  
  params <- get_default_params(model_type)
  params <- apply_class_weights(params, data_split_pool$training_data, data_split_pool$target_col, setting, model_type)
  
  # Hyperparameter tuning for appropriate settings
  if (setting %in% c("tuning_with_weights", "tuning_no_weights")) {
    best_params <- perform_hyperparameter_tuning(params, pools$train_pool, 
                                                 pools$validation_pool, 
                                                 data_split_pool$validation_data,
                                                 data_split_pool$target_col,
                                                 model_type)
  } else {
    best_params <- params
    cat("Using default parameters - LR:", best_params$learning_rate, "Depth:", best_params$depth, "\n")
  }
  
  # Train final model oh
  set.seed(123)
  best_params$iterations <- 1000  # Reset to full iterations
  model <- catboost.train(pools$train_pool, params = best_params, test_pool = pools$validation_pool)
  
  list(model = model, best_params = best_params)
}

# MODEL EVALUATION FUNCTIONS　####

evaluate_model <- function(model, data_splits, target, model_type) {
  pools <- create_catboost_pools(
    data_splits$training_data, 
    data_splits$validation_data, 
    data_splits$test_data, 
    data_splits$target_col, 
    data_splits$fold_col
  )
  
  if (model_type == "regression") {
    # Regression predictions
    predicted_values <- catboost.predict(model, pools$test_pool)
    predicted_values_val <- catboost.predict(model, pools$validation_pool)
    
    # Calculate regression metrics
    metrics <- calculate_regression_metrics(predicted_values, data_splits$test_data[[target]])
    
    return(list(
      predicted_values = predicted_values,
      predicted_values_val = predicted_values_val,
      metrics = metrics
    ))
    
  } else if (model_type == "binary") {
    # Binary classification predictions
    predicted_probabilities <- catboost.predict(model, pools$test_pool, prediction_type = "Probability")
    predicted_probabilities_val <- catboost.predict(model, pools$validation_pool, prediction_type = "Probability")
    predicted_classes <- ifelse(predicted_probabilities > 0.5, 1, 0)
    
    # Calculate binary classification metrics
    metrics <- calculate_binary_classification_metrics(predicted_classes, data_splits$test_data[[target]], 
                                                       predicted_probabilities)
    
    return(list(
      predicted_probabilities = predicted_probabilities,
      predicted_probabilities_val = predicted_probabilities_val,
      predicted_classes = predicted_classes,
      metrics = metrics
    ))
    
  } else if (model_type == "multiclass") {
    # Multiclass predictions
    predicted_probabilities <- catboost.predict(model, pools$test_pool, prediction_type = "Probability")
    predicted_probabilities_val <- catboost.predict(model, pools$validation_pool, prediction_type = "Probability")
    predicted_classes <- catboost.predict(model, pools$test_pool, prediction_type = "Class")
    
    # Calculate multiclass metrics
    metrics <- calculate_multiclass_metrics(predicted_classes, data_splits$test_data[[target]], 
                                            predicted_probabilities)
    
    return(list(
      predicted_probabilities = predicted_probabilities,
      predicted_probabilities_val = predicted_probabilities_val,
      predicted_classes = predicted_classes,
      metrics = metrics
    ))
  }
}

# REGRESSION METRICS #####

calculate_regression_metrics <- function(predicted_values, actual_values) {
  n <- length(actual_values)
  residuals <- actual_values - predicted_values
  
  # Basic metrics
  mae <- mean(abs(residuals))
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  
  # R-squared and adjusted R-squared
  ss_res <- sum(residuals^2)
  ss_tot <- sum((actual_values - mean(actual_values))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  
  # Mean absolute percentage error (handling division by zero)
  mape <- ifelse(any(actual_values == 0), NA, 
                 mean(abs(residuals / actual_values)) * 100)
  
  # Median absolute error
  median_ae <- median(abs(residuals))
  
  # Mean bias (systematic error)
  mean_bias <- mean(residuals)
  
  # Standard deviation of residuals
  residual_sd <- sd(residuals)
  
  list(
    n_samples = n,
    mae = mae,
    mse = mse,
    rmse = rmse,
    r_squared = r_squared,
    mape = mape,
    median_ae = median_ae,
    mean_bias = mean_bias,
    residual_sd = residual_sd,
    min_residual = min(residuals),
    max_residual = max(residuals)
  )
}

# BINARY CLASSIFICATION METRICS #####

calculate_binary_classification_metrics <- function(predicted_classes, actual_classes, predicted_probabilities) {
  # Convert to factors
  predicted_classes <- factor(predicted_classes, levels = c(0, 1))
  y_test <- factor(actual_classes, levels = c(0, 1))
  
  # Confusion matrix
  conf_matrix <- table(Predicted = predicted_classes, Actual = y_test)
  
  # Extract confusion matrix components
  confusion_components <- extract_confusion_matrix_components(conf_matrix)
  
  # Basic metrics
  accuracy <- mean(predicted_classes == y_test)
  
  # ROC metrics
  roc_metrics <- calculate_roc_metrics(actual_classes, predicted_probabilities)
  
  # PR metrics
  pr_metrics <- calculate_pr_metrics(actual_classes, predicted_probabilities)
  
  # Binary classification metrics
  binary_metrics <- calculate_binary_metrics(confusion_components)
  
  # Advanced metrics
  advanced_metrics <- calculate_advanced_metrics(roc_metrics$auc, pr_metrics$pr_auc, 
                                                 binary_metrics, actual_classes)
  
  # Model calibration metrics
  calibration_metrics <- calculate_calibration_metrics(predicted_probabilities, actual_classes)
  
  # Class distribution metrics
  class_metrics <- calculate_class_distribution_metrics(actual_classes, predicted_classes)
  
  # Combine all metrics
  c(
    list(
      accuracy = accuracy,
      confusion_matrix = conf_matrix,
      n_samples = length(actual_classes)
    ),
    roc_metrics,
    pr_metrics,
    binary_metrics,
    advanced_metrics,
    calibration_metrics,
    class_metrics
  )
}

# MULTICLASS METRICS #####

calculate_multiclass_metrics <- function(predicted_classes, actual_classes, predicted_probabilities) {
  # Convert to factors with same levels
  all_classes <- sort(unique(c(actual_classes, predicted_classes)))
  predicted_classes <- factor(predicted_classes, levels = all_classes)
  actual_classes <- factor(actual_classes, levels = all_classes)
  
  # Confusion matrix
  conf_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)
  
  # Overall accuracy
  accuracy <- mean(predicted_classes == actual_classes)
  
  # Per-class metrics
  n_classes <- length(all_classes)
  class_metrics <- list()
  
  # Calculate per-class precision, recall, F1
  for (i in 1:n_classes) {
    class_name <- all_classes[i]
    
    # True positives, false positives, false negatives for this class
    tp <- conf_matrix[i, i]
    fp <- sum(conf_matrix[i, ]) - tp
    fn <- sum(conf_matrix[, i]) - tp
    
    # Per-class metrics
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    f1 <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
    
    class_metrics[[paste0("class_", class_name, "_precision")]] <- precision
    class_metrics[[paste0("class_", class_name, "_recall")]] <- recall
    class_metrics[[paste0("class_", class_name, "_f1")]] <- f1
  }
  
  # Macro averages (unweighted average across classes)
  precisions <- sapply(1:n_classes, function(i) {
    tp <- conf_matrix[i, i]
    fp <- sum(conf_matrix[i, ]) - tp
    ifelse(tp + fp > 0, tp / (tp + fp), 0)
  })
  recalls <- sapply(1:n_classes, function(i) {
    tp <- conf_matrix[i, i]
    fn <- sum(conf_matrix[, i]) - tp
    ifelse(tp + fn > 0, tp / (tp + fn), 0)
  })
  f1s <- sapply(1:n_classes, function(i) {
    tp <- conf_matrix[i, i]
    fp <- sum(conf_matrix[i, ]) - tp
    fn <- sum(conf_matrix[, i]) - tp
    precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
    recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
    ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
  })
  
  macro_precision <- mean(precisions)
  macro_recall <- mean(recalls)
  macro_f1 <- mean(f1s)
  
  # Micro averages (weighted by class frequency)
  total_tp <- sum(diag(conf_matrix))
  total_fp <- sum(conf_matrix) - total_tp
  total_fn <- total_fp  # In multiclass, total FP = total FN
  
  micro_precision <- ifelse(total_tp + total_fp > 0, total_tp / (total_tp + total_fp), 0)
  micro_recall <- ifelse(total_tp + total_fn > 0, total_tp / (total_tp + total_fn), 0)
  micro_f1 <- ifelse(micro_precision + micro_recall > 0, 
                     2 * micro_precision * micro_recall / (micro_precision + micro_recall), 0)
  
  # Weighted averages (weighted by class support)
  class_supports <- as.numeric(table(actual_classes))
  weighted_precision <- sum(precisions * class_supports) / sum(class_supports)
  weighted_recall <- sum(recalls * class_supports) / sum(class_supports)
  weighted_f1 <- sum(f1s * class_supports) / sum(class_supports)
  
  # ROC/PR metrics and Kappa
  roc_pr_kappa_metrics <- calculate_multiclass_roc_pr_kappa(predicted_probabilities, actual_classes, conf_matrix, all_classes)
  
  # Combine all metrics
  c(
    list(
      accuracy = accuracy,
      confusion_matrix = conf_matrix,
      n_samples = length(actual_classes),
      n_classes = n_classes,
      macro_precision = macro_precision,
      macro_recall = macro_recall,
      macro_f1 = macro_f1,
      micro_precision = micro_precision,
      micro_recall = micro_recall,
      micro_f1 = micro_f1,
      weighted_precision = weighted_precision,
      weighted_recall = weighted_recall,
      weighted_f1 = weighted_f1
    ),
    class_metrics,
    roc_pr_kappa_metrics
  )
}

# HELPER FUNCTIONS #####

calculate_multiclass_roc_pr_kappa <- function(predicted_probabilities, actual_classes, conf_matrix, all_classes) {
  n_classes <- length(all_classes)
  
  # Initialize per-class metrics
  class_roc_auc <- numeric(n_classes)
  class_pr_auc <- numeric(n_classes)
  class_normalized_pr_auc <- numeric(n_classes)
  
  # Calculate per-class ROC and PR metrics
  for (i in 1:n_classes) {
    class_name <- all_classes[i]
    
    # Create binary classification: class i vs all others
    binary_actual <- as.numeric(actual_classes == class_name)
    binary_probs <- predicted_probabilities[, i]
    
    # ROC-AUC
    tryCatch({
      roc_obj <- roc(binary_actual, binary_probs, direction = "auto")
      class_roc_auc[i] <- as.numeric(auc(roc_obj))
    }, error = function(e) {
      class_roc_auc[i] <<- NA
    })
    
    # PR-AUC and Normalized PR-AUC
    tryCatch({
      positives <- binary_probs[binary_actual == 1]
      negatives <- binary_probs[binary_actual == 0]
      
      if (length(positives) > 0 && length(negatives) > 0) {
        pr_curve <- pr.curve(scores.class0 = positives, scores.class1 = negatives)
        class_pr_auc[i] <- pr_curve$auc.integral
        
        # Normalized PR-AUC
        positive_rate <- mean(binary_actual)
        if (positive_rate > 0 && positive_rate < 1) {
          class_normalized_pr_auc[i] <- (class_pr_auc[i] - positive_rate) / (1 - positive_rate)
        } else {
          class_normalized_pr_auc[i] <- NA
        }
      } else {
        class_pr_auc[i] <- NA
        class_normalized_pr_auc[i] <- NA
      }
    }, error = function(e) {
      class_pr_auc[i] <<- NA
      class_normalized_pr_auc[i] <<- NA
    })
  }
  
  # Calculate macro averages
  macro_roc_auc <- mean(class_roc_auc, na.rm = TRUE)
  macro_pr_auc <- mean(class_pr_auc, na.rm = TRUE)
  macro_normalized_pr_auc <- mean(class_normalized_pr_auc, na.rm = TRUE)
  
  # Calculate weighted averages
  class_supports <- as.numeric(table(actual_classes))
  weighted_roc_auc <- ifelse(all(!is.na(class_roc_auc)), 
                             sum(class_roc_auc * class_supports, na.rm = TRUE) / sum(class_supports), 
                             NA)
  weighted_pr_auc <- ifelse(all(!is.na(class_pr_auc)), 
                            sum(class_pr_auc * class_supports, na.rm = TRUE) / sum(class_supports), 
                            NA)
  weighted_normalized_pr_auc <- ifelse(all(!is.na(class_normalized_pr_auc)), 
                                       sum(class_normalized_pr_auc * class_supports, na.rm = TRUE) / sum(class_supports), 
                                       NA)
  
  # Cohen's Kappa
  kappa <- calculate_multiclass_kappa(conf_matrix)
  
  # Create per-class metric lists
  class_roc_metrics <- setNames(as.list(class_roc_auc), paste0("class_", all_classes, "_roc_auc"))
  class_pr_metrics <- setNames(as.list(class_pr_auc), paste0("class_", all_classes, "_pr_auc"))
  class_normalized_pr_metrics <- setNames(as.list(class_normalized_pr_auc), paste0("class_", all_classes, "_normalized_pr_auc"))
  
  # Combine all ROC/PR/Kappa metrics
  c(
    list(
      # Aggregate ROC metrics
      macro_roc_auc = macro_roc_auc,
      weighted_roc_auc = weighted_roc_auc,
      
      # Aggregate PR metrics
      macro_pr_auc = macro_pr_auc,
      weighted_pr_auc = weighted_pr_auc,
      
      # Aggregate normalized PR metrics
      macro_normalized_pr_auc = macro_normalized_pr_auc,
      weighted_normalized_pr_auc = weighted_normalized_pr_auc,
      
      # Cohen's Kappa
      kappa = kappa
    ),
    # Per-class metrics
    class_roc_metrics,
    class_pr_metrics,
    class_normalized_pr_metrics
  )
}

calculate_multiclass_kappa <- function(conf_matrix) {
  # Calculate Cohen's Kappa for multiclass classification
  n <- sum(conf_matrix)
  
  if (n == 0) return(NA)
  
  # Observed accuracy
  observed_accuracy <- sum(diag(conf_matrix)) / n
  
  # Expected accuracy by chance
  row_marginals <- rowSums(conf_matrix) / n
  col_marginals <- colSums(conf_matrix) / n
  expected_accuracy <- sum(row_marginals * col_marginals)
  
  # Cohen's Kappa
  if (expected_accuracy >= 1) {
    return(NA)  # Perfect agreement by chance
  }
  
  kappa <- (observed_accuracy - expected_accuracy) / (1 - expected_accuracy)
  return(kappa)
}

calculate_roc_metrics <- function(actual_classes, predicted_probabilities) {
  tryCatch({
    roc_obj <- roc(actual_classes, predicted_probabilities, direction = "auto")
    auc_value <- as.numeric(auc(roc_obj))
    
    # Calculate confidence interval for AUC
    auc_ci <- ci.auc(roc_obj, conf.level = 0.95)
    
    list(
      auc = auc_value,
      auc_lower_ci = as.numeric(auc_ci[1]),
      auc_upper_ci = as.numeric(auc_ci[3])
    )
  }, error = function(e) {
    list(auc = NA, auc_lower_ci = NA, auc_upper_ci = NA)
  })
}

calculate_pr_metrics <- function(actual_classes, predicted_probabilities) {
  tryCatch({
    positives <- predicted_probabilities[actual_classes == 1]
    negatives <- predicted_probabilities[actual_classes == 0]
    
    if (length(positives) > 0 && length(negatives) > 0) {
      pr_curve <- pr.curve(scores.class0 = positives, scores.class1 = negatives, curve = TRUE)
      pr_auc <- pr_curve$auc.integral
      
      # Calculate baseline (random classifier) PR AUC
      positive_rate <- mean(actual_classes)
      baseline_pr_auc <- positive_rate
      
      # Normalized PR AUC
      normalized_pr_auc <- ifelse(positive_rate >= 1, NA, (pr_auc - positive_rate) / (1 - positive_rate))
      
      list(
        pr_auc = pr_auc,
        baseline_pr_auc = baseline_pr_auc,
        normalized_pr_auc = normalized_pr_auc,
        positive_rate = positive_rate
      )
    } else {
      list(pr_auc = NA, baseline_pr_auc = NA, normalized_pr_auc = NA, positive_rate = mean(actual_classes))
    }
  }, error = function(e) {
    list(pr_auc = NA, baseline_pr_auc = NA, normalized_pr_auc = NA, positive_rate = mean(actual_classes))
  })
}

calculate_binary_metrics <- function(components) {
  TP <- components$TP; TN <- components$TN; FP <- components$FP; FN <- components$FN
  
  # Basic rates
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)  # Recall, TPR
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)  # TNR
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)    # PPV
  npv <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)          # Negative Predictive Value
  
  # Derived metrics
  recall <- sensitivity
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), NA)
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  # Additional rates
  fpr <- ifelse((TN + FP) > 0, FP / (TN + FP), NA)  # False Positive Rate
  fnr <- ifelse((TP + FN) > 0, FN / (TP + FN), NA)  # False Negative Rate
  fdr <- ifelse((TP + FP) > 0, FP / (TP + FP), NA)  # False Discovery Rate
  
  # Likelihood ratios
  lr_positive <- ifelse(fpr > 0, sensitivity / fpr, NA)     # Positive Likelihood Ratio
  lr_negative <- ifelse(specificity > 0, fnr / specificity, NA)  # Negative Likelihood Ratio
  
  # Diagnostic odds ratio
  dor <- ifelse(lr_negative > 0, lr_positive / lr_negative, NA)
  
  # Matthews Correlation Coefficient
  mcc_denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- ifelse(mcc_denominator > 0, (TP * TN - FP * FN) / mcc_denominator, NA)
  
  list(
    sensitivity = sensitivity, specificity = specificity, precision = precision, npv = npv,
    recall = recall, f1_score = f1_score, balanced_accuracy = balanced_accuracy,
    fpr = fpr, fnr = fnr, fdr = fdr,
    lr_positive = lr_positive, lr_negative = lr_negative, dor = dor, mcc = mcc
  )
}

calculate_advanced_metrics <- function(auc, pr_auc, binary_metrics, actual_classes) {
  # Somers' D (Informedness): 2*AUC - 1
  somers_d <- ifelse(is.na(auc), NA, 2 * auc - 1)
  
  # Youden's J statistic (Informedness)
  youden_j <- ifelse(is.na(binary_metrics$sensitivity) || is.na(binary_metrics$specificity), 
                     NA, binary_metrics$sensitivity + binary_metrics$specificity - 1)
  
  # Markedness
  markedness <- ifelse(is.na(binary_metrics$precision) || is.na(binary_metrics$npv),
                       NA, binary_metrics$precision + binary_metrics$npv - 1)
  
  # Cohen's Kappa
  kappa <- calculate_cohens_kappa(actual_classes, binary_metrics)
  
  # Geometric mean of sensitivity and specificity
  gmean <- ifelse(is.na(binary_metrics$sensitivity) || is.na(binary_metrics$specificity),
                  NA, sqrt(binary_metrics$sensitivity * binary_metrics$specificity))
  
  # Fowlkes-Mallows index
  fm_index <- ifelse(is.na(binary_metrics$sensitivity) || is.na(binary_metrics$precision),
                     NA, sqrt(binary_metrics$sensitivity * binary_metrics$precision))
  
  list(
    somers_d = somers_d, youden_j = youden_j, markedness = markedness,
    kappa = kappa, gmean = gmean, fm_index = fm_index
  )
}

calculate_cohens_kappa <- function(actual_classes, binary_metrics) {
  # Calculate from sensitivity, specificity, and class distribution
  pos_rate <- mean(actual_classes)
  neg_rate <- 1 - pos_rate
  
  # Observed agreement
  observed_agreement <- binary_metrics$sensitivity * pos_rate + binary_metrics$specificity * neg_rate
  
  # Expected agreement (assuming independence)
  pred_pos_rate <- binary_metrics$sensitivity * pos_rate + (1 - binary_metrics$specificity) * neg_rate
  pred_neg_rate <- 1 - pred_pos_rate
  
  expected_agreement <- pos_rate * pred_pos_rate + neg_rate * pred_neg_rate
  
  # Cohen's Kappa
  ifelse(expected_agreement < 1, 
         (observed_agreement - expected_agreement) / (1 - expected_agreement), 
         NA)
}

calculate_calibration_metrics <- function(predicted_probabilities, actual_classes) {
  # Brier Score
  brier_score <- mean((predicted_probabilities - actual_classes)^2)
  
  # Brier Score components (reliability, resolution, uncertainty)
  n_bins <- 10
  bin_edges <- seq(0, 1, length.out = n_bins + 1)
  
  reliability <- 0
  resolution <- 0
  n_total <- length(actual_classes)
  overall_rate <- mean(actual_classes)
  
  for (i in 1:n_bins) {
    # Find samples in this bin
    if (i == n_bins) {
      in_bin <- predicted_probabilities >= bin_edges[i] & predicted_probabilities <= bin_edges[i + 1]
    } else {
      in_bin <- predicted_probabilities >= bin_edges[i] & predicted_probabilities < bin_edges[i + 1]
    }
    
    if (sum(in_bin) > 0) {
      bin_prob <- mean(predicted_probabilities[in_bin])
      bin_rate <- mean(actual_classes[in_bin])
      bin_count <- sum(in_bin)
      
      # Reliability (calibration error)
      reliability <- reliability + (bin_count / n_total) * (bin_prob - bin_rate)^2
      
      # Resolution
      resolution <- resolution + (bin_count / n_total) * (bin_rate - overall_rate)^2
    }
  }
  
  uncertainty <- overall_rate * (1 - overall_rate)
  
  list(
    brier_score = brier_score,
    reliability = reliability,
    resolution = resolution,
    uncertainty = uncertainty,
    brier_skill_score = 1 - brier_score / uncertainty
  )
}

calculate_class_distribution_metrics <- function(actual_classes, predicted_classes) {
  # Class distributions
  actual_pos <- sum(actual_classes == 1)
  actual_neg <- sum(actual_classes == 0)
  pred_pos <- sum(predicted_classes == 1)
  pred_neg <- sum(predicted_classes == 0)
  
  # Class imbalance ratio
  imbalance_ratio <- ifelse(actual_neg > 0, actual_pos / actual_neg, NA)
  
  # Prediction bias
  prediction_bias <- ifelse(actual_pos > 0, pred_pos / actual_pos, NA)
  
  list(
    actual_positive_count = actual_pos,
    actual_negative_count = actual_neg,
    predicted_positive_count = pred_pos,
    predicted_negative_count = pred_neg,
    imbalance_ratio = imbalance_ratio,
    prediction_bias = prediction_bias
  )
}

extract_confusion_matrix_components <- function(conf_matrix) {
  if (nrow(conf_matrix) == 2 && ncol(conf_matrix) == 2) {
    list(TN = conf_matrix[1, 1], FP = conf_matrix[1, 2], 
         FN = conf_matrix[2, 1], TP = conf_matrix[2, 2])
  } else {
    list(
      TN = ifelse("0" %in% rownames(conf_matrix) && "0" %in% colnames(conf_matrix), conf_matrix["0", "0"], 0),
      TP = ifelse("1" %in% rownames(conf_matrix) && "1" %in% colnames(conf_matrix), conf_matrix["1", "1"], 0),
      FN = ifelse("0" %in% colnames(conf_matrix) && "1" %in% rownames(conf_matrix), conf_matrix["1", "0"], 0),
      FP = ifelse("1" %in% colnames(conf_matrix) && "0" %in% rownames(conf_matrix), conf_matrix["0", "1"], 0)
    )
  }
}

# FEATURE ANALYSIS FUNCTIONS #######

extract_feature_importance <- function(model) {
  # Use model-based importance for global feature utility (complements SHAP's local explanations)
  importance <- catboost.get_feature_importance(model, type = "FeatureImportance")
  
  data.frame(
    feature = rownames(importance),
    importance = importance[,1],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


calculate_shap_values <- function(model, data_splits, model_type = "binary") {
  pools <- create_catboost_pools(
    data_splits$training_data, 
    data_splits$validation_data, 
    data_splits$test_data, 
    data_splits$target_col, 
    data_splits$fold_col
  )
  
  # Calculate SHAP values
  shap_values <- catboost.get_feature_importance(model, pool = pools$test_pool, type = "ShapValues")
  
  # Get feature names from the pool
  feature_names <- names(data_splits$test_data)[!names(data_splits$test_data) %in% c(data_splits$target_col, data_splits$fold_col)]
  
  if (model_type == "multiclass") {
    # For multiclass: shap_values is [n_samples, n_classes, n_features]
    # Remove bias values (last feature for each class)
    n_samples <- dim(shap_values)[1]
    n_features <- dim(shap_values)[3] - 1  # Remove bias
    n_classes <- dim(shap_values)[2]
    
    # Remove bias column (last column in each class)
    shap_matrix <- shap_values[, , 1:n_features, drop = FALSE]
    
    # Get class names/numbers
    unique_classes <- sort(unique(data_splits$test_data[[data_splits$target_col]]))
    
    # Create a list to store SHAP values for each class
    shap_list <- list()
    
    for (class_idx in 1:n_classes) {
      class_name <- unique_classes[class_idx]
      
      # Extract SHAP values for this class
      class_shap <- shap_matrix[, class_idx, , drop = FALSE]
      class_shap_df <- as.data.frame(class_shap[, 1, ])
      
      # Set column names
      colnames(class_shap_df) <- feature_names
      
      # Add metadata
      class_shap_df$sample_id <- rownames(data_splits$test_data)
      class_shap_df$actual <- data_splits$test_data[[data_splits$target_col]]
      class_shap_df$class <- class_name
      
      shap_list[[paste0("class_", class_name)]] <- class_shap_df
    }
    
    # Combine all classes into one dataframe with class identifier
    shap_df <- do.call(rbind, shap_list)
    
    return(list(
      shap_by_class = shap_list,
      shap_combined = shap_df,
      n_classes = n_classes,
      class_names = unique_classes
    ))
    
  } else {
    # For binary classification and regression: shap_values is [n_samples, n_features]
    # Remove bias column (last column)
    shap_matrix <- shap_values[, -ncol(shap_values), drop = FALSE]
    
    # Create SHAP dataframe
    shap_df <- as.data.frame(shap_matrix)
    colnames(shap_df) <- feature_names
    
    # Add metadata
    shap_df$sample_id <- rownames(data_splits$test_data)
    shap_df$actual <- data_splits$test_data[[data_splits$target_col]]
    
    return(shap_df)
  }
}

# PREDICTION DATAFRAME FUNCTIONS #######

create_predictions_dataframe <- function(data_splits, evaluation, model_type) {
  if (model_type == "regression") {
    # Test predictions
    test_predictions <- data.frame(
      sample_id = rownames(data_splits$test_data),
      actual = data_splits$test_data[[data_splits$target_col]],
      predicted_value = evaluation$predicted_values,
      split = "test",
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    # Validation predictions
    validation_predictions <- data.frame(
      sample_id = rownames(data_splits$validation_data),
      actual = data_splits$validation_data[[data_splits$target_col]],
      predicted_value = evaluation$predicted_values_val,
      split = "valid",
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
  } else {
    # Classification predictions
    test_predictions <- data.frame(
      sample_id = rownames(data_splits$test_data),
      actual = data_splits$test_data[[data_splits$target_col]],
      predicted_class = as.character(evaluation$predicted_classes),
      split = "test",
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    validation_predictions <- data.frame(
      sample_id = rownames(data_splits$validation_data),
      actual = data_splits$validation_data[[data_splits$target_col]],
      predicted_class = "NA",  # Would need validation predictions
      split = "valid",
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    # Add probabilities for binary classification
    if (model_type == "binary") {
      test_predictions$predicted_prob <- evaluation$predicted_probabilities
      validation_predictions$predicted_prob <- evaluation$predicted_probabilities_val
    } else if (model_type == "multiclass") {
      # For multiclass, could add probability matrix here if needed
      test_predictions$predicted_prob_max <- apply(evaluation$predicted_probabilities, 1, max)
      validation_predictions$predicted_prob_max <- apply(evaluation$predicted_probabilities_val, 1, max)
    }
  }
  
  dplyr::bind_rows(test_predictions, validation_predictions)
}

# Model Loading and Prediction Function #######
load_trained_model <- function(results_dir) {
  file_paths <- setup_file_paths(results_dir)
  if (!file.exists(file_paths$model_file)) {
    stop("Model file not found")
  }
  catboost.load_model(file_paths$model_file)
}

predict_new_data <- function(model_path_or_object, new_data, target_col = NULL, fold_col = NULL, model_type = "binary") {
  if (is.character(model_path_or_object)) {
    model <- catboost.load_model(model_path_or_object)
  } else {
    model <- model_path_or_object
  }
  
  exclude_cols <- c(target_col, fold_col)
  feature_data <- new_data[, !names(new_data) %in% exclude_cols, drop = FALSE]
  
  new_pool <- catboost.load_pool(data = feature_data)
  
  if (model_type == "regression") {
    catboost.predict(model, new_pool)
  } else if (model_type == "binary") {
    catboost.predict(model, new_pool, prediction_type = "Probability")
  } else if (model_type == "multiclass") {
    catboost.predict(model, new_pool, prediction_type = "Class")
  }
}

library(dplyr)
library(tidyr)

# MAIN COLLECTION AND SUMMARY FUNCTION #####
collect_cv_results <- function(param_config, dat_for_ml) {
  #"""
  #Collect results from all cross-validation folds and create summary statistics
  #
  #Args:
  #  param_config: Configuration list containing OUTLOC, fold_col, model_type
  #  dat_for_ml: Original dataset used for training
  #
  #Returns:
  #  List containing summarized metrics, feature importance, and SHAP values
  #"""
  
  # Get unique folds
  unique_folds <- unique(dat_for_ml[[param_config$fold_col]])
  
  # Initialize collections
  all_evaluations <- list()
  all_feature_importance <- list()
  all_shap_values <- list()
  all_predictions <- list()
  
  # Collect results from each fold
  cat("Collecting results from", length(unique_folds), "folds...\n")
  
  for (fold_indx in unique_folds) {
    fold_dir <- paste0(param_config$OUTLOC, "submodels/model_", fold_indx, "/")
    
    # Check if fold results exist
    if (!dir.exists(fold_dir)) {
      warning(paste("Fold", fold_indx, "directory not found:", fold_dir))
      next
    }
    
    # Load fold results
    fold_results <- load_fold_results(fold_dir, fold_indx)
    
    if (!is.null(fold_results$evaluation)) {
      all_evaluations[[paste0("fold_", fold_indx)]] <- fold_results$evaluation
    }
    
    if (!is.null(fold_results$feature_importance)) {
      fold_results$feature_importance$fold <- fold_indx
      all_feature_importance[[paste0("fold_", fold_indx)]] <- fold_results$feature_importance
    }
    
    if (!is.null(fold_results$shap_values)) {
      all_shap_values[[paste0("fold_", fold_indx)]] <- add_fold_to_shap(fold_results$shap_values, fold_indx, param_config$model_type)
    }
    
    if (!is.null(fold_results$predictions)) {
      fold_results$predictions$fold <- fold_indx
      all_predictions[[paste0("fold_", fold_indx)]] <- fold_results$predictions
    }
  }
  
  # Create summaries
  cat("Creating summary statistics...\n")
  
  evaluation_summary <- summarize_evaluations(all_evaluations, param_config$model_type)
  feature_importance_summary <- summarize_feature_importance(all_feature_importance)
  shap_summary <- summarize_shap_values(all_shap_values, param_config$model_type)
  #predictions_combined <- combine_predictions(all_predictions)
  
  # Save combined results
  save_cv_summaries(
    evaluation_summary, 
    feature_importance_summary, 
    shap_summary, 
    # predictions_combined,
    param_config$OUTLOC
  )
  
  return(list(
    evaluation_summary = evaluation_summary,
    feature_importance_summary = feature_importance_summary,
    shap_summary = shap_summary,
    #predictions_combined = predictions_combined,
    raw_results = list(
      evaluations = all_evaluations,
      feature_importance = all_feature_importance,
      shap_values = all_shap_values,
      predictions = all_predictions
    )
  ))
}

# HELPER FUNCTIONS FOR LOADING FOLD RESULTS #####
load_fold_results <- function(fold_dir, fold_indx) {
  # """Load results from a single fold directory"""
  
  file_paths <- list(
    evaluation = file.path(fold_dir, "evaluation.rds"),
    feature_importance = file.path(fold_dir, "feature_importance.rds"),
    shap_values = file.path(fold_dir, "shap_values.rds"),
    predictions = file.path(fold_dir, "predictions.rds")
  )
  
  results <- list()
  
  # Load evaluation
  if (file.exists(file_paths$evaluation)) {
    results$evaluation <- readRDS(file_paths$evaluation)
  } else {
    warning(paste("Evaluation file not found for fold", fold_indx))
  }
  
  # Load feature importance
  if (file.exists(file_paths$feature_importance)) {
    results$feature_importance <- readRDS(file_paths$feature_importance)
  } else {
    warning(paste("Feature importance file not found for fold", fold_indx))
  }
  
  # Load SHAP values
  if (file.exists(file_paths$shap_values)) {
    results$shap_values <- readRDS(file_paths$shap_values)
  } else {
    warning(paste("SHAP values file not found for fold", fold_indx))
  }
  
  # Load predictions
  if (file.exists(file_paths$predictions)) {
    results$predictions <- readRDS(file_paths$predictions)
  } else {
    warning(paste("Predictions file not found for fold", fold_indx))
  }
  
  return(results)
}

add_fold_to_shap <- function(shap_values, fold_indx, model_type) {
  #"""Add fold identifier to SHAP values"""
  
  if (model_type == "multiclass" && is.list(shap_values) && "shap_combined" %in% names(shap_values)) {
    shap_values$shap_combined$fold <- fold_indx
    # Also add to individual class dataframes
    for (class_name in names(shap_values$shap_by_class)) {
      shap_values$shap_by_class[[class_name]]$fold <- fold_indx
    }
  } else if (is.data.frame(shap_values)) {
    shap_values$fold <- fold_indx
  }
  
  return(shap_values)
}

# EVALUATION SUMMARY FUNCTIONS #####
summarize_evaluations <- function(all_evaluations, model_type) {
  # """Summarize evaluation metrics across all folds"""
  
  if (length(all_evaluations) == 0) {
    warning("No evaluation results found")
    return(NULL)
  }
  
  if (model_type == "regression") {
    return(summarize_regression_metrics(all_evaluations))
  } else if (model_type == "binary") {
    return(summarize_binary_metrics(all_evaluations))
  } else if (model_type == "multiclass") {
    return(summarize_multiclass_metrics(all_evaluations))
  }
}

summarize_regression_metrics <- function(all_evaluations) {
  # """Summarize regression metrics"""
  
  # Extract numeric metrics
  metric_names <- c("mae", "mse", "rmse", "r_squared", "mape", "median_ae", 
                    "mean_bias", "residual_sd", "min_residual", "max_residual")
  
  metrics_df <- data.frame(fold = character(), stringsAsFactors = FALSE)
  
  for (metric in metric_names) {
    metrics_df[[metric]] <- numeric()
  }
  
  for (fold_name in names(all_evaluations)) {
    fold_metrics <- all_evaluations[[fold_name]]$metrics
    
    new_row <- data.frame(fold = fold_name, stringsAsFactors = FALSE)
    for (metric in metric_names) {
      new_row[[metric]] <- ifelse(metric %in% names(fold_metrics), 
                                  fold_metrics[[metric]], NA)
    }
    metrics_df <- rbind(metrics_df, new_row)
  }
  
  # Calculate summary statistics
  summary_stats <- data.frame(
    metric = metric_names,
    mean = sapply(metric_names, function(m) mean(metrics_df[[m]], na.rm = TRUE)),
    sd = sapply(metric_names, function(m) sd(metrics_df[[m]], na.rm = TRUE)),
    min = sapply(metric_names, function(m) min(metrics_df[[m]], na.rm = TRUE)),
    max = sapply(metric_names, function(m) max(metrics_df[[m]], na.rm = TRUE)),
    q25 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.25, na.rm = TRUE)),
    median = sapply(metric_names, function(m) median(metrics_df[[m]], na.rm = TRUE)),
    q75 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  return(list(
    fold_metrics = metrics_df,
    summary_stats = summary_stats
  ))
}

summarize_binary_metrics <- function(all_evaluations) {
  #"""Summarize binary classification metrics"""
  
  # Key metrics to summarize
  metric_names <- c("accuracy", "auc", "pr_auc", "sensitivity", "specificity", 
                    "precision", "recall", "f1_score", "balanced_accuracy",
                    "mcc", "youden_j", "brier_score","kappa","somers_d","normalized_pr_auc")
  
  metrics_df <- data.frame(fold = character(), stringsAsFactors = FALSE)
  
  for (metric in metric_names) {
    metrics_df[[metric]] <- numeric()
  }
  
  for (fold_name in names(all_evaluations)) {
    fold_metrics <- all_evaluations[[fold_name]][["metrics"]]
    
    new_row <- data.frame(fold = fold_name, stringsAsFactors = FALSE)
    for (metric in metric_names) {
      new_row[[metric]] <- ifelse(metric %in% names(fold_metrics), 
                                  fold_metrics[[metric]], NA)
    }
    metrics_df <- rbind(metrics_df, new_row)
  }
  
  # Calculate summary statistics
  summary_stats <- data.frame(
    metric = metric_names,
    mean = sapply(metric_names, function(m) mean(metrics_df[[m]], na.rm = TRUE)),
    sd = sapply(metric_names, function(m) sd(metrics_df[[m]], na.rm = TRUE)),
    min = sapply(metric_names, function(m) min(metrics_df[[m]], na.rm = TRUE)),
    max = sapply(metric_names, function(m) max(metrics_df[[m]], na.rm = TRUE)),
    q25 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.25, na.rm = TRUE)),
    median = sapply(metric_names, function(m) median(metrics_df[[m]], na.rm = TRUE)),
    q75 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  return(list(
    fold_metrics = metrics_df,
    summary_stats = summary_stats
  ))
}

summarize_multiclass_metrics <- function(all_evaluations) {
  #"""Summarize multiclass classification metrics"""
  
  # Key metrics to summarize - now includes ROC/PR/Kappa metrics
  metric_names <- c("accuracy", "macro_precision", "macro_recall", "macro_f1",
                    "micro_precision", "micro_recall", "micro_f1",
                    "weighted_precision", "weighted_recall", "weighted_f1",
                    "macro_roc_auc", "weighted_roc_auc", 
                    "macro_pr_auc", "weighted_pr_auc",
                    "macro_normalized_pr_auc", "weighted_normalized_pr_auc",
                    "kappa")
  
  metrics_df <- data.frame(fold = character(), stringsAsFactors = FALSE)
  
  for (metric in metric_names) {
    metrics_df[[metric]] <- numeric()
  }
  
  # Collect per-class metrics dynamically
  per_class_metrics <- list()
  
  for (fold_name in names(all_evaluations)) {
    fold_metrics <- all_evaluations[[fold_name]][["metrics"]]
    
    new_row <- data.frame(fold = fold_name, stringsAsFactors = FALSE)
    for (metric in metric_names) {
      new_row[[metric]] <- ifelse(metric %in% names(fold_metrics), 
                                  fold_metrics[[metric]], NA)
    }
    metrics_df <- rbind(metrics_df, new_row)
    
    # Collect per-class metrics (class_X_precision, class_X_recall, etc.)
    class_metric_names <- names(fold_metrics)[grepl("^class_.*_(precision|recall|f1|roc_auc|pr_auc|normalized_pr_auc)$", names(fold_metrics))]
    
    for (class_metric in class_metric_names) {
      if (!class_metric %in% names(per_class_metrics)) {
        per_class_metrics[[class_metric]] <- list()
      }
      per_class_metrics[[class_metric]][[fold_name]] <- fold_metrics[[class_metric]]
    }
  }
  
  # Calculate summary statistics for main metrics
  summary_stats <- data.frame(
    metric = metric_names,
    mean = sapply(metric_names, function(m) mean(metrics_df[[m]], na.rm = TRUE)),
    sd = sapply(metric_names, function(m) sd(metrics_df[[m]], na.rm = TRUE)),
    min = sapply(metric_names, function(m) min(metrics_df[[m]], na.rm = TRUE)),
    max = sapply(metric_names, function(m) max(metrics_df[[m]], na.rm = TRUE)),
    q25 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.25, na.rm = TRUE)),
    median = sapply(metric_names, function(m) median(metrics_df[[m]], na.rm = TRUE)),
    q75 = sapply(metric_names, function(m) quantile(metrics_df[[m]], 0.75, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  
  # Calculate summary statistics for per-class metrics
  per_class_summary <- NULL
  if (length(per_class_metrics) > 0) {
    per_class_summary_list <- list()
    
    for (class_metric_name in names(per_class_metrics)) {
      values <- unlist(per_class_metrics[[class_metric_name]])
      
      per_class_summary_list[[class_metric_name]] <- data.frame(
        metric = class_metric_name,
        mean = mean(values, na.rm = TRUE),
        sd = sd(values, na.rm = TRUE),
        min = min(values, na.rm = TRUE),
        max = max(values, na.rm = TRUE),
        q25 = quantile(values, 0.25, na.rm = TRUE),
        median = median(values, na.rm = TRUE),
        q75 = quantile(values, 0.75, na.rm = TRUE),
        n_folds = length(values[!is.na(values)]),
        stringsAsFactors = FALSE
      )
    }
    
    per_class_summary <- do.call(rbind, per_class_summary_list)
    rownames(per_class_summary) <- NULL
    
    # Sort by metric type and class
    per_class_summary <- per_class_summary[order(per_class_summary$metric), ]
  }
  
  # Create extended metrics dataframe with per-class metrics
  extended_metrics_df <- metrics_df
  for (class_metric_name in names(per_class_metrics)) {
    extended_metrics_df[[class_metric_name]] <- NA
    
    for (fold_name in names(per_class_metrics[[class_metric_name]])) {
      fold_idx <- which(extended_metrics_df$fold == fold_name)
      if (length(fold_idx) > 0) {
        extended_metrics_df[[class_metric_name]][fold_idx] <- per_class_metrics[[class_metric_name]][[fold_name]]
      }
    }
  }
  
  return(list(
    fold_metrics = metrics_df,                    # Main aggregate metrics by fold
    extended_fold_metrics = extended_metrics_df,  # All metrics (including per-class) by fold
    summary_stats = summary_stats,                # Summary of main aggregate metrics
    per_class_summary = per_class_summary         # Summary of per-class metrics
  ))
}
# FEATURE IMPORTANCE SUMMARY FUNCTIONS #####
summarize_feature_importance <- function(all_feature_importance) {
  # """Summarize feature importance across all folds"""
  
  if (length(all_feature_importance) == 0) {
    warning("No feature importance results found")
    return(NULL)
  }
  
  # Combine all feature importance dataframes
  combined_importance <- do.call(rbind, all_feature_importance)
  
  # Calculate summary statistics by feature
  importance_summary <- combined_importance %>%
    group_by(feature) %>%
    summarise(
      mean_importance = mean(importance, na.rm = TRUE),
      sd_importance = sd(importance, na.rm = TRUE),
      min_importance = min(importance, na.rm = TRUE),
      max_importance = max(importance, na.rm = TRUE),
      q25_importance = quantile(importance, 0.25, na.rm = TRUE),
      median_importance = median(importance, na.rm = TRUE),
      q75_importance = quantile(importance, 0.75, na.rm = TRUE),
      cv_importance = sd_importance / mean_importance,  # Coefficient of variation
      n_folds = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_importance))
  
  return(list(
    fold_importance = combined_importance,
    summary_stats = importance_summary
  ))
}

# SHAP VALUES SUMMARY FUNCTIONS #####
summarize_shap_values <- function(all_shap_values, model_type) {
  #"""Summarize SHAP values across all folds"""
  
  if (length(all_shap_values) == 0) {
    warning("No SHAP values found")
    return(NULL)
  }
  
  if (model_type == "multiclass") {
    return(summarize_multiclass_shap(all_shap_values))
  } else {
    return(summarize_binary_regression_shap(all_shap_values))
  }
}

summarize_binary_regression_shap <- function(all_shap_values) {
  #"""Summarize SHAP values for binary classification or regression"""
  
  # Combine all SHAP dataframes
  combined_shap <- do.call(rbind, all_shap_values)
  
  # Get feature columns (exclude metadata columns)
  meta_cols <- c("sample_id", "actual", "fold")
  feature_cols <- setdiff(names(combined_shap), meta_cols)
  
  # Calculate summary statistics by feature
  shap_summary <- combined_shap %>%
    select(all_of(feature_cols)) %>%
    pivot_longer(cols = everything(), names_to = "feature", values_to = "shap_value") %>%
    group_by(feature) %>%
    summarise(
      mean_shap = mean(shap_value, na.rm = TRUE),
      sd_shap = sd(shap_value, na.rm = TRUE),
      mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
      median_shap = median(shap_value, na.rm = TRUE),
      q25_shap = quantile(shap_value, 0.25, na.rm = TRUE),
      q75_shap = quantile(shap_value, 0.75, na.rm = TRUE),
      min_shap = min(shap_value, na.rm = TRUE),
      max_shap = max(shap_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_shap))
  
  return(list(
    fold_shap = combined_shap,
    summary_stats = shap_summary
  ))
}

summarize_multiclass_shap <- function(all_shap_values) {
  #"""Summarize SHAP values for multiclass classification"""
  
  # Extract combined SHAP dataframes
  combined_shap_list <- lapply(all_shap_values, function(x) x$shap_combined)
  combined_shap <- do.call(rbind, combined_shap_list)
  
  # Get feature columns
  meta_cols <- c("sample_id", "actual", "class", "fold")
  feature_cols <- setdiff(names(combined_shap), meta_cols)
  
  # Summary by feature and class
  shap_summary <- combined_shap %>%
    select(all_of(c(feature_cols, "class"))) %>%
    pivot_longer(cols = all_of(feature_cols), names_to = "feature", values_to = "shap_value") %>%
    group_by(feature, class) %>%
    summarise(
      mean_shap = mean(shap_value, na.rm = TRUE),
      sd_shap = sd(shap_value, na.rm = TRUE),
      mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
      median_shap = median(shap_value, na.rm = TRUE),
      q25_shap = quantile(shap_value, 0.25, na.rm = TRUE),
      q75_shap = quantile(shap_value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(class, desc(mean_abs_shap))
  
  # Overall summary across all classes
  overall_summary <- combined_shap %>%
    select(all_of(feature_cols)) %>%
    pivot_longer(cols = everything(), names_to = "feature", values_to = "shap_value") %>%
    group_by(feature) %>%
    summarise(
      mean_abs_shap = mean(abs(shap_value), na.rm = TRUE),
      sd_abs_shap = sd(abs(shap_value), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_shap))
  
  return(list(
    fold_shap = combined_shap,
    summary_by_class = shap_summary,
    overall_summary = overall_summary
  ))
}

# PREDICTIONS COMBINATION FUNCTION #####
combine_predictions <- function(all_predictions) {
  #"""Combine predictions from all folds"""
  
  if (length(all_predictions) == 0) {
    warning("No predictions found")
    return(NULL)
  }
  
  combined_predictions <- do.call(rbind, all_predictions)
  
  return(combined_predictions)
}

# SAVE FUNCTIONS #####
save_cv_summaries <- function(evaluation_summary, feature_importance_summary, 
                              shap_summary, output_dir) {
  # """Save all summary results to files"""
  
  cv_summary_dir <- file.path(output_dir, "cv_summary")
  dir.create(cv_summary_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save evaluation summaries
  if (!is.null(evaluation_summary)) {
    saveRDS(evaluation_summary, file.path(cv_summary_dir, "evaluation_summary.rds"))
    write.csv(evaluation_summary$summary_stats, 
              file.path(cv_summary_dir, "evaluation_summary.csv"), 
              row.names = FALSE)
  }
  
  # Save feature importance summaries
  if (!is.null(feature_importance_summary)) {
    saveRDS(feature_importance_summary, file.path(cv_summary_dir, "feature_importance_summary.rds"))
    write.csv(feature_importance_summary$summary_stats, 
              file.path(cv_summary_dir, "feature_importance_summary.csv"), 
              row.names = FALSE)
  }
  
  # Save SHAP summaries
  if (!is.null(shap_summary)) {
    saveRDS(shap_summary, file.path(cv_summary_dir, "shap_summary.rds"))
    if ("summary_stats" %in% names(shap_summary)) {
      write.csv(shap_summary$summary_stats, 
                file.path(cv_summary_dir, "shap_summary.csv"), 
                row.names = FALSE)
    }
    if ("overall_summary" %in% names(shap_summary)) {
      write.csv(shap_summary$overall_summary, 
                file.path(cv_summary_dir, "shap_overall_summary.csv"), 
                row.names = FALSE)
    }
  }
  
  cat("CV summary results saved to:", cv_summary_dir, "\n")
}