#' Train a classifier on perturbed expression data and predict classes for new samples
#'
#' @description
#' `get_classifier()` trains a supervised classification model to discriminate compound
#' groups using perturbed gene expression profiles, and returns the predicted class
#' for one or more new samples provided in `test_matrix`.
#'
#' The function subsets expression data by compound groups (from `...`) and optionally
#' by `probes`, `dose`, and/or `time`, fits the requested model on the subset, then
#' predicts classes for the supplied test expression profile(s).
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`).
#'   These define the class labels for model training.
#' @param ge_matrix Gene expression matrix with probes/genes in rows and samples (barcodes)
#'   in columns.
#' @param metadata A data frame or tibble containing sample metadata corresponding to
#'   `ge_matrix`. Must include the fields required by `get_subset()` (e.g., barcode,
#'   compound name, dose level, time level), and group membership must be derivable
#'   from `...`.
#' @param probes Optional vector of probe/gene identifiers to subset rows of `ge_matrix`.
#'   If `NULL`, all probes are used.
#' @param dose Optional character vector of dose levels to subset the data. If `NULL`,
#'   all dose levels are retained.
#' @param time Optional character vector of time points to subset the data. If `NULL`,
#'   all time points are retained.
#' @param model Character string specifying the classification algorithm to use.
#'   One of `"lr"` (logistic regression), `"svm"` (support vector machine with RBF kernel),
#'   `"rf"` (random forest), or `"knn"` (k-nearest neighbors). Default is `"lr"`.
#' @param test_matrix Expression profile(s) to classify. Either:
#' \itemize{
#'   \item a named numeric vector (one sample), where names are probe/gene identifiers, or
#'   \item a numeric matrix (one or more samples) with samples in rows and probes/genes in columns.
#' }
#' Column names must match the predictors used to train the model (after subsetting by
#' `probes`, `dose`, and `time`). If `test_matrix` is a matrix, row names are used as
#' sample IDs when present.
#' @param error_call The environment used for error reporting. Default is
#'   `rlang::caller_env()`.
#'
#' @details
#' The function validates inputs, subsets data using `get_subset()`, and builds a
#' modeling table by transposing the subset expression matrix (samples as rows,
#' probes/genes as columns) and attaching the group label. A classification model is
#' fit on all available training data and used to predict the class for `test_matrix`.
#' Only class predictions are returned (no resampling or performance metrics).
#'
#' @return A tibble with one row per test sample and two columns:
#' \itemize{
#'   \item `sample_id`: sample identifier (row names of `test_matrix` if provided, otherwise
#'   `"test_sample1"` for vector input or `"test_sample<k>"` for matrix input without row names)
#'   \item `predicted_class`: predicted class label
#' }
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 45, n_ee = 5, n_com = c(10,10))
#' test_sample <- sample(1 : nrow(sim_data$metadata), 10)
#' test_mat <- t(sim_data$expression[,test_sample])
#' gr <- list(A = paste0("Compound", 1:10), B = paste0("Compound", 11:20))
#' get_classifier(gr,
#'   ge_matrix = sim_data$expression[, -test_sample],
#'   metadata = sim_data$metadata[-test_sample, ],
#'   probes = paste0("DE", 1:20),
#'   test_matrix = test_mat)
#'
#' @seealso
#' [get_subset()] for data subsetting and preprocessing; [test_group()] for group parsing;
#' [parsnip::logistic_reg()], [parsnip::svm_rbf()], [parsnip::rand_forest()],
#' [parsnip::nearest_neighbor()] for supported models.
#'
#' @export
get_classifier <- function(...,
                           ge_matrix,
                           test_matrix,
                           metadata,
                           probes = NULL,
                           dose = NULL,
                           time = NULL,
                           model = c("lr", "svm", "rf", "knn"),
                           error_call = rlang::caller_env()) {
  comps_group <- test_group(...)
  test_data(ge_matrix, metadata)

  dt <- get_subset(comps_group,
                   ge_matrix = ge_matrix,
                   metadata = metadata,
                   probes = probes,
                   dose = dose,
                   time = time,
                   error_call = error_call)

  expr_tbl <- tibble::as_tibble(t(dt$expression)) %>%
    dplyr::mutate(group = as.factor(dt$metadata$group))

  model <- test_input(model, auto_input = TRUE)

  if (identical(model, "lr")) {
    mod <- parsnip::logistic_reg() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("glm")
  } else if (identical(model, "svm")) {
    mod <- parsnip::svm_rbf() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("kernlab")
  } else if (identical(model, "rf")) {
    mod <- parsnip::rand_forest() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("ranger", importance = "impurity")
  } else if (identical(model, "knn")) {
    mod <- parsnip::nearest_neighbor(neighbors = 4) %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("kknn")
  }

  mod_wf <- workflows::workflow() %>%
    workflows::add_formula(group ~ .) %>%
    workflows::add_model(mod)

  final_fit <- mod_wf %>%
    parsnip::fit(data = expr_tbl)

  predictors <- setdiff(names(expr_tbl), "group")

  if (is.null(test_matrix)) {
    rlang::abort("`test_matrix` must be provided.", call = error_call)
  }

  if (is.numeric(test_matrix) && is.vector(test_matrix)) {
    if (is.null(names(test_matrix))) {
      rlang::abort(
        "`test_matrix` as a vector must be a *named* numeric vector.",
        call = error_call
      )
    }
    new_df <- as.data.frame(t(test_matrix))
    sample_id <- "test_sample1"

  } else if (is.matrix(test_matrix)) {
    new_df <- as.data.frame(test_matrix)
    sample_id <- rownames(test_matrix)

  } else {
    rlang::abort(
      "`test_matrix` must be a named numeric vector or a numeric matrix.",
      call = error_call
    )
  }

  if (is.null(colnames(new_df))) {
    rlang::abort(
      "`test_matrix` must have column names matching training predictors.",
      call = error_call
    )
  }

  if (is.null(sample_id)) {
    sample_id <- paste0("test_sample", seq_len(nrow(new_df)))
  }

  missing_cols <- setdiff(predictors, colnames(new_df))
  if (length(missing_cols) > 0) {
    rlang::abort(
      paste0(
        "Missing predictors in `test_matrix`: ",
        paste(missing_cols, collapse = ", ")
      ),
      call = error_call
    )
  }

  new_df <- new_df[, predictors, drop = FALSE]

  pred_class <- predict(final_fit, new_data = new_df, type = "class")

  pred_tab <- tibble::tibble(
    sample_id = sample_id,
    predicted_class = pred_class$.pred_class
  )

  return(pred_tab)
}



#' Train and cross-validate a toxicity classifier on perturbed expression profiles
#'
#' @description
#' `tox_eval()` builds a supervised classification model to discriminate compound groups
#' using perturbed gene expression profiles. The function subsets expression data by compound
#' groups and optionally by `probes`, `dose`, and/or `time`, fits a specified classifier, and
#' evaluates predictive performance using repeated v-fold cross-validation. Performance is
#' summarized with an ROC curve and AUC estimated from cross-validation predictions.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`).
#'   These define the class labels for model training.
#' @param ge_matrix Gene expression matrix with probes/genes in rows and samples (barcodes)
#'   in columns.
#' @param metadata A data frame or tibble containing sample metadata corresponding to
#'   `ge_matrix`. Must include the fields required by `get_subset()` (e.g., barcode,
#'   compound name, dose level, time level) and allow group membership to be derived from
#'   `...`.
#' @param probes Optional vector of probe/gene identifiers used to subset rows of `ge_matrix`.
#'   If `NULL`, all probes are used.
#' @param dose Optional character vector of dose levels to subset the data. If `NULL`,
#'   all dose levels are retained.
#' @param time Optional character vector of time points to subset the data. If `NULL`,
#'   all time points are retained.
#' @param nfold Integer. Number of folds used in v-fold cross-validation. Default is `10`.
#' @param nrep Integer. Number of repeats for v-fold cross-validation. Default is `2`.
#' @param model Character string specifying the classification algorithm to use.
#'   One of `"lr"` (logistic regression), `"svm"` (support vector machine with RBF kernel),
#'   `"rf"` (random forest), or `"knn"` (k-nearest neighbors). Default is `"lr"`.
#' @param error_call The environment used for error reporting. Default is
#'   `rlang::caller_env()`.
#'
#' @details
#' The function first subsets data using `get_subset()`. The resulting expression matrix is
#' transposed to a sample-by-feature table and combined with the class label (`group`).
#' A model is specified with **parsnip** and fit/evaluated using **tune** repeated v-fold
#' cross-validation (`rsample::vfold_cv()`).
#'
#' ROC and AUC are computed from the saved resampling predictions using **ROCR**. The current
#' implementation assumes a *binary* classification problem; labels are formed by treating the
#' first factor level of `group` as the negative class and the second level as the positive
#' class.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `perf`: a ROCR performance object for the ROC curve (TPR vs FPR)
#'   \item `auc`: a list of AUC values (`ROCR::performance(..., "auc")@y.values`) from the
#'   resampling predictions
#' }
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 45, n_ee = 5, n_com = c(10,10))
#' gr <- list(A = paste0("Compound", 1:10), B = paste0("Compound", 11:20))
#' res <- tox_eval(gr,
#'   ge_matrix = sim_data$expression,
#'   metadata = sim_data$metadata,
#'   probes = paste0("DE", 1:20),
#'   nfold = 5,
#'   nrep = 10)
#'
#' @seealso
#' [get_subset()] for data subsetting and preprocessing; [test_group()] for group parsing;
#' [rsample::vfold_cv()] for resampling;
#' [tune::fit_resamples()] for cross-validated fitting;
#' [ROCR::performance()] and [ROCR::prediction()] for ROC/AUC calculation;
#' [parsnip::logistic_reg()], [parsnip::svm_rbf()], [parsnip::rand_forest()],
#' [parsnip::nearest_neighbor()] for supported models.
#'
#' @export
tox_eval <- function(...,
                   ge_matrix,
                   metadata,
                   probes = NULL,
                   dose = NULL,
                   time = NULL,
                   nfold = 10,
                   nrep = 2,
                   model = "lr",
                   error_call = caller_env()) {
  comps_gr <- test_group(...)
  test_input(model, c("lr", "svm", "rf", "knn"))
  dt <- get_subset(comps_gr,
                   ge_matrix = ge_matrix,
                   metadata = metadata,
                   probes = probes,
                   dose = dose, time = time,
                   error_call = error_call)
  expr_tbl <- tibble::as_tibble(t(dt$expression)) %>%
    dplyr::mutate(group = as.factor(dt$metadata$group))
  if (identical(model, "lr")) {
    mod <- parsnip::logistic_reg() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("glm")
  } else if (identical(model, "rf")) {
    mod <- parsnip::rand_forest() %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("ranger", importance = "impurity")
  } else if (identical(model, "svm")) {
    mod <- parsnip::svm_rbf() %>%
           parsnip::set_mode("classification") %>%
           parsnip::set_engine("kernlab")
  } else if (identical(model, "knn")) {
    mod <- parsnip::nearest_neighbor(neighbors = 4) %>%
      parsnip::set_mode("classification") %>%
      parsnip::set_engine("kknn")
  }
  mod_wf <- workflows::workflow() %>%
            workflows::add_formula(group ~ .) %>%
            workflows::add_model(mod)
  folds <-  rsample::vfold_cv(expr_tbl, v = nfold, repeats = nrep, strata = group)
  mod_fit <- mod_wf %>%
    tune::fit_resamples(resamples = folds,
                  control = tune::control_resamples(save_pred = TRUE,
                                              verbose = TRUE)) %>%
    suppressMessages()

  labels <- lapply(mod_fit$.predictions, function(x)ifelse(x$group == levels(x$group)[1], 1, 0))
  predictions <- lapply(mod_fit$.predictions, function(x) dplyr::pull(x[,3]))
  pred <- ROCR::prediction(predictions, labels)
  auc_roc <- ROCR::performance(pred, measure = "auc")
  auc <- auc_roc@y.values
  perf <- ROCR::performance(pred, "tpr", "fpr")
  return(list(perf, auc))
}


#' Plot ROC curves for cross-validated toxicity classification
#'
#' @description
#' `tox_roc()` generates a receiver operating characteristic (ROC) plot summarizing
#' cross-validated prediction performance of a toxicity classification model. The
#' function calls [tox_eval()] to train and evaluate the classifier using repeated
#' v-fold cross-validation, then visualizes individual and average ROC curves together
#' with the mean AUC and its standard deviation.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see
#'   [test_group()]). These define the class labels used for model training and
#'   evaluation.
#' @param ge_matrix Gene expression matrix with probes/genes in rows and samples
#'   (barcodes) in columns.
#' @param metadata A data frame or tibble containing sample metadata corresponding
#'   to `ge_matrix`. Must include the fields required by [get_subset()] (e.g., barcode,
#'   compound name, dose level, time level).
#' @param probes Optional vector of probe/gene identifiers to subset rows of
#'   `ge_matrix`. If `NULL`, all probes are used.
#' @param dose Optional character vector of dose levels to subset the data. If `NULL`,
#'   all dose levels are retained.
#' @param time Optional character vector of time points to subset the data. If `NULL`,
#'   all time points are retained.
#' @param nfold Integer. Number of folds used in v-fold cross-validation. Default is `5`.
#' @param nrep Integer. Number of repeats for v-fold cross-validation. Default is `10`.
#' @param model Character string specifying the classification algorithm to use.
#'   One of `"lr"` (logistic regression), `"svm"` (support vector machine with RBF kernel),
#'   `"rf"` (random forest), or `"knn"` (k-nearest neighbors). Default is `"lr"`.
#' @param error_call The environment used for error reporting. Default is
#'   `rlang::caller_env()`.
#'
#' @details
#' The function relies on [tox_eval()] to compute ROC performance objects and AUC values
#' from repeated cross-validation. ROC curves from individual resamples are plotted as
#' dashed lines, while the average ROC curve (threshold-averaged) is overlaid as a solid
#' line. The legend reports the mean AUC and its standard deviation across resamples.
#'
#' The current implementation assumes a binary classification problem.
#'
#' @return
#' A `ggplot` object containing the ROC curve visualization. The plot can be further
#' customized using standard **ggplot2** methods.
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 45, n_ee = 5, n_com = c(10,10))
#' gr <- list(A = paste0("Compound", 1:10), B = paste0("Compound", 11:20))
#' roc_plot <- tox_roc(gr,
#'   ge_matrix = sim_data$expression,
#'   metadata = sim_data$metadata,
#'   probes = paste0("DE", 1:20),
#'   nfold = 5,
#'   nrep = 10)
#'
#' @seealso
#' [tox_eval()] for cross-validated model evaluation;
#' [ROCR::performance()] for ROC computation;
#' [ggplotify::as.ggplot()] for base-to-ggplot conversion.
#'
#' @export
tox_roc <- function(...,
                    ge_matrix,
                    metadata,
                    probes = NULL,
                    dose = NULL,
                    time = NULL,
                    nfold = 5,
                    nrep = 10,
                    model = "lr",
                    error_call = rlang::caller_env()) {
  comps_gr <- test_group(...)
  roc_res <- tox_eval(comps_gr,
                      ge_matrix = ge_matrix,
                      metadata = metadata,
                      probes = probes,
                      dose = dose,
                      time = time,
                      nfold = nfold,
                      nrep = nrep,
                      model = model,
                      error_call = error_call)

  perf <- roc_res[[1]]
  auc  <- roc_res[[2]]
  auc_mean <- mean(unlist(auc))
  auc_sd   <- stats::sd(unlist(auc))
  roc_plot <- ggplotify::as.ggplot(function() {
    graphics::par(cex.lab = 1.2)
    plot(perf, lty = 3, col = "#FB9A99", cex.axis = 1)
    plot(perf, avg = "threshold", add = TRUE, col = "#E31A1C", lwd = 3, main = "")
    graphics::legend(x = 0.54, y = 0.10,
      legend = bquote(bar(AUC)~"(SD)"==.(sprintf("%.2f", auc_mean))~
                        "(" * .(sprintf("%.2f", auc_sd)) * ")"),
      col = "#E31A1C", lwd = 3, bty = "n", cex = 1, y.intersp = 1.3
    )
  })
  return(roc_plot)
}

