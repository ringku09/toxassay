# #library(tidymodels)
#
#
# get_classifier <- function(...,
#                            ge_matrix,
#                            metadata,
#                            probes = NULL,
#                            dose = NULL,
#                            time = NULL,
#                            model = c("lr", "svm", "xgboost", "rf", "knn"),
#                            multicore = FALSE,
#                            store = FALSE,
#                            output_dir = missing_arg(),
#                            error_call = caller_env()) {
#   test_input(model, auto_input = TRUE)
#   comps_group <- test_group(...)
#   test_data(ge_matrix, metadata)
#   dt <- get_subset(comps_gr,
#                    ge_matrix = ge_matrix,
#                    metadata = metadata,
#                    probes = probes,
#                    dose = dose,
#                    time = time,
#                    multicore = multicore,
#                    store = store,
#                    output_dir = output_dir,
#                    error_call = error_call)
#   expr_tbl <- tibble::as_tibble(t(dt$expression)) %>%
#     dplyr::mutate(group = as.factor(dt$metadata$group))
#   if (identical(model, "svm")) {
#     mod <-
#       svm_rbf() %>%
#       set_mode("classification") %>%
#       set_engine("kernlab")
#   } else if (identical(model, "xgboost")) {
#     mod <-
#       boost_tree(learn_rate = 0.05,
#                  mtry = 0.5,
#                  min_n = 1,
#                  loss_reduction = 1,
#                  sample_size = 1,
#                  tree_depth = 4,
#                  trees = 200) %>%
#       set_engine("xgboost") %>%
#       set_mode("classification")
#   } else if (identical(model, "lr")) {
#     mod <-
#       logistic_reg() %>%
#       set_mode("classification") %>%
#       set_engine("glm")
#     # logistic_reg(penalty =  0.01684069, mixture = 0.1934516) %>%
#     # set_mode("classification") %>%
#     # set_engine("glmnet")
#   } else if (identical(model, "rf")) {
#     mod <-
#       rand_forest() %>%
#       set_mode("classification") %>%
#       set_engine("ranger", importance = "impurity")
#   } else if (identical(model, "knn")) {
#     mod <-
#       nearest_neighbor(neighbors = 4) %>%
#       set_mode("classification") %>%
#       set_engine("kknn")
#   }
#   mod_wf <- workflow() %>%
#     # add_variables(outcomes = group, predictors = everything()) %>%
#     add_formula(group ~ .) %>%
#     add_model(mod)
#   # parm_grid <- expand_grid(penalty = c(0,1,2),
#   #                               mixture = seq(0,1,by=0.2))
#   folds <-  vfold_cv(expr_tbl, v = nfold, repeats = nrep, strata = group)
#   mod_fit <- mod_wf %>%
#     fit_resamples(resamples = folds, control = control_resamples(save_pred = TRUE, verbose = TRUE)) %>%
#     suppressMessages()
#
#   labels <- lapply(mod_fit$.predictions, function(x)ifelse(x$group == levels(x$group)[1], 1, 0))
#   predictions <- lapply(mod_fit$.predictions, function(x) dplyr::pull(x[,1]))
#   pred <- prediction(predictions, labels)
#   auc_roc <- performance(pred, measure = "auc")
#   auc <- auc_roc@y.values
#   perf <- performance(pred, "tpr", "fpr")
#   return(list(perf, auc))
# }
#
# mod_perf <- function(...,
#                      expr_data,
#                      attr_data,
#                      probes = NULL,
#                      dose = NULL,
#                      time = NULL,
#                      nfold = 10,
#                      nrep = 1,
#                      multicore = FALSE,
#                      store = FALSE,
#                      output_dir = missing_arg(),
#                      error_call = caller_env()) {
#   dt <- get_space(...,
#                   expr_data = expr_data,
#                   attr_data = attr_data,
#                   probes = probes,
#                   dose = dose,
#                   time = time,
#                   multicore = multicore,
#                   store = store,
#                   output_dir = output_dir,
#                   error_call = error_call)
#
#   expr_tbl <- as_tibble(t(dt$expr_data)) %>%
#     mutate(group = as.factor(dt$attr_data$group))
#   folds <-  vfold_cv(expr_tbl, v = nfold, repeats = nrep, strata = group)
#
#   # Models
#   rf_mod <-
#     rand_forest() %>%
#     set_engine("ranger") %>%
#     set_mode("classification")
#   log_mod <- # your model specification
#     logistic_reg() %>%  # model type
#     set_engine(engine = "glm") %>%  # model engine
#     set_mode("classification") # model mode
#   knn_mod <-
#     nearest_neighbor( ) %>% # we can adjust the number of neighbors
#     set_engine("kknn") %>%
#     set_mode("classification")
#   svm_mod <-
#     svm_rbf() %>%
#     set_mode("classification") %>%
#     set_engine("kernlab")
#
#
#   log_wf <-
#     workflow() %>%
#     add_model(log_mod) %>%
#     add_formula(group ~ .)
#
#   knn_wf <-
#     workflow() %>%
#     add_model(knn_mod) %>%
#     add_formula(group ~ .)
#   rf_wf <-
#     workflow() %>%
#     add_model(rf_mod) %>%
#     add_formula(group ~ .)
#   svm_wf <-
#     workflow() %>%
#     add_model(svm_mod) %>%
#     add_formula(group ~ .)
#
#
#   rf_fit <-
#     rf_wf %>%
#     fit_resamples(resamples = folds,
#                   metrics = metric_set(recall, precision, f_meas, accuracy, kap,roc_auc, sens, spec),
#                   control = control_resamples(save_pred = TRUE))
#   log_fit <-
#     log_wf %>%
#     fit_resamples(resamples = folds,
#                   metrics = metric_set(recall, precision, f_meas, accuracy, kap,roc_auc, sens, spec),
#                   control = control_resamples(save_pred = TRUE)) %>%
#     suppressMessages()
#   knn_fit <-
#     knn_wf %>%
#     fit_resamples(resamples = folds,
#                   metrics = metric_set(recall, precision, f_meas, accuracy, kap,roc_auc, sens, spec),
#                   control = control_resamples(save_pred = TRUE))
#   svm_fit <-
#     svm_wf %>%
#     fit_resamples(resamples = folds,
#                   metrics = metric_set(recall, precision, f_meas, accuracy, kap,roc_auc, sens, spec),
#                   control = control_resamples(save_pred = TRUE))
#
#   # Evaluation
#   acc <- data.frame(
#   RF = glue("{round(collect_metrics(rf_fit)$mean,3)} ({round(collect_metrics(rf_fit)$std_err,3)})"),
#   KNN = glue("{round(collect_metrics(log_fit)$mean,3)} ({round(collect_metrics(log_fit)$std_err,3)})"),
#   SVM = glue("{round(collect_metrics(knn_fit)$mean,3)} ({round(collect_metrics(knn_fit)$std_err,3)})"),
#   LR = glue("{round(collect_metrics(svm_fit)$mean,3)} ({round(collect_metrics(svm_fit)$std_err,3)})")
#   )
#   rownames(acc) <- c("Recall", "Precision", "F_measure", "Accuracy",
#                      "Kappa", "AUC", "Sensitivity", "Specificity")
#   return(acc)
# }
#
# # mod_perf(
# #   comps_gr,
# #   expr_data = expr_data,
# #   attr_data = attr_data,
# #   probes = all_probe,
# #   dose = NULL,
# #   time = "24 hr",
# #   nrep = 10
# # )
#
#
#
# # install.packages("ranger")
# #---------------------------------------  ROCR    --------------------------
# library(ROCR)
# library(ggplotify)
# library(xgboost)
# library(ranger)
#
# roc_cv <- function(...,
#                    ge_matrix,
#                    metadata,
#                    probes = NULL,
#                    dose = NULL,
#                    time = NULL,
#                    nfold = 10,
#                    nrep = 2,
#                    model = "lr",
#                    multicore = FALSE,
#                    store = FALSE,
#                    output_dir = missing_arg(),
#                    error_call = caller_env()) {
#   test_input(model, c("lr", "svm", "xgboost", "rf", "knn"))
#   dt <- get_subset(comps_gr,
#                          ge_matrix = ge_matrix,
#                          metadata = metadata,
#                          probes = probes,
#                          dose = dose,
#                          time = time,
#                          multicore = multicore,
#                          store = store,
#                          output_dir = output_dir,
#                          error_call = error_call)
#   expr_tbl <- as_tibble(t(dt$expression)) %>%
#     mutate(group = as.factor(dt$metadata$group))
#   if (identical(model, "svm")) {
#     mod <-
#       svm_rbf() %>%
#       set_mode("classification") %>%
#       set_engine("kernlab")
#   } else if (identical(model, "xgboost")) {
#     mod <-
#       boost_tree(learn_rate = 0.05,
#                  mtry = 0.5,
#                  min_n = 1,
#                  loss_reduction = 1,
#                  sample_size = 1,
#                  tree_depth = 4,
#                  trees = 200) %>%
#       set_engine("xgboost") %>%
#       set_mode("classification")
#   } else if (identical(model, "lr")) {
#     mod <-
#       logistic_reg() %>%
#       set_mode("classification") %>%
#       set_engine("glm")
#       # logistic_reg(penalty =  0.01684069, mixture = 0.1934516) %>%
#       # set_mode("classification") %>%
#       # set_engine("glmnet")
#   } else if (identical(model, "rf")) {
#     mod <-
#       rand_forest() %>%
#       set_mode("classification") %>%
#       set_engine("ranger", importance = "impurity")
#   } else if (identical(model, "knn")) {
#     mod <-
#       nearest_neighbor(neighbors = 4) %>%
#       set_mode("classification") %>%
#       set_engine("kknn")
#   }
#   mod_wf <- workflow() %>%
#    # add_variables(outcomes = group, predictors = everything()) %>%
#     add_formula(group ~ .) %>%
#     add_model(mod)
#   # parm_grid <- expand_grid(penalty = c(0,1,2),
#   #                               mixture = seq(0,1,by=0.2))
#   folds <-  vfold_cv(expr_tbl, v = nfold, repeats = nrep, strata = group)
#   mod_fit <- mod_wf %>%
#     fit_resamples(resamples = folds, control = control_resamples(save_pred = TRUE, verbose = TRUE)) %>%
#     suppressMessages()
#
#   labels <- lapply(mod_fit$.predictions, function(x)ifelse(x$group == levels(x$group)[1], 1, 0))
#   predictions <- lapply(mod_fit$.predictions, function(x) dplyr::pull(x[,1]))
#   pred <- prediction(predictions, labels)
#   auc_roc <- performance(pred, measure = "auc")
#   auc <- auc_roc@y.values
#   perf <- performance(pred, "tpr", "fpr")
#   return(list(perf, auc))
# }
#
#
# compare_roc <- function(...,
#                         ge_matrix,
#                         metadata,
#                         probe_list = NULL,
#                         dose = NULL,
#                         time = NULL,
#                         nfold = 5,
#                         nrep = 10,
#                         model = "lr",
#                         multicore = FALSE,
#                         store = FALSE,
#                         output_dir = missing_arg(),
#                         error_call = caller_env()) {
#   n_compare <- length(probe_list)
#   roc_perf <- vector(mode = "list", length = n_compare)
#   for (i in 1:n_compare) {
#     roc_perf[[i]] <- roc_cv(comps_gr,
#                             ge_matrix = ge_matrix,
#                             metadata = metadata,
#                             probes = probe_list[[i]],
#                             dose = dose,
#                             time = time,
#                             nfold = nfold,
#                             nrep = nrep,
#                             model = model,
#                             multicore = multicore,
#                             store = store,
#                             output_dir = output_dir,
#                             error_call = error_call)
#   }
#   roc_compare <- as.ggplot(function() {
#     par(cex.lab=1.2)
#     plot(roc_perf[[6]][[1]], lty=3,col="#CD853F", cex.axis=5)
#     plot(roc_perf[[5]][[1]], lty=3, col="#CAB2D6", add=TRUE)
#     plot(roc_perf[[4]][[1]], lty=3, col="#87CEFA", add=TRUE)
#     plot(roc_perf[[3]][[1]], lty=3, col="#90EE90", add=TRUE)
#     plot(roc_perf[[2]][[1]], lty=3, col="#DDA0DD", add=TRUE)
#     plot(roc_perf[[1]][[1]], lty=3, col="#FB9A99", add=TRUE)
#     plot(roc_perf[[6]][[1]], avg = "threshold", add = TRUE, col = "#8B4513", lwd = 3, main = "")
#     plot(roc_perf[[5]][[1]], avg = "threshold", add = TRUE, col = "#6A3D9A", lwd = 3, main = "")
#     plot(roc_perf[[4]][[1]], avg = "threshold", add = TRUE, col = "#1F78B4", lwd = 3, main = "")
#     plot(roc_perf[[3]][[1]],avg= "threshold",add = TRUE,col= "#33A02C",lwd= 3,main= "")
#     plot(roc_perf[[2]][[1]],avg= "threshold",add = TRUE,col= "#FF00FF",lwd= 3,main= "")
#     plot(roc_perf[[1]][[1]],avg= "threshold",add = TRUE,col= "#E31A1C",lwd= 3,main= "")
#     legend(x = 0.54, y = 0.36, legend = c(bquote("           "~bar(AUC)~"(SD)") ,
#                                          bquote("Core  " ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[1]][[-1]]))))~(.(sprintf("%.2f", sd(unlist(roc_perf[[1]][[-1]])))))),
#                                          bquote("DEGs" ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[2]][[-1]]))))~(.(sprintf("%.2f", sd(unlist(roc_perf[[2]][[-1]])))))),
#                                          bquote("Set-A" ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[3]][[-1]]))))~(.(sprintf("%.2f", sd(unlist(roc_perf[[3]][[-1]])))))),
#                                          bquote("Set-B" ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[4]][[-1]]))))~(.(sprintf("%.2f", sd(unlist(roc_perf[[4]][[-1]])))))),
#                                          bquote("Set-C" ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[5]][[-1]]))))~(.(sprintf("%.2f", sd(unlist(roc_perf[[5]][[-1]])))))),
#                                          bquote("Set-D" ==
#                                                   .(sprintf("%.2f", mean(unlist(roc_perf[[6]][[-1]]))))~(.(round(sd(unlist(roc_perf[[6]][[-1]])),2))))),
#            col= c("white","#E31A1C", "#FF00FF", "#33A02C", "#1F78B4","#6A3D9A" , "#8B4513"),
#            lwd = 3, bty = "n", cex = 1, y.intersp= 1.3)})
#   dev.off()
#   return(roc_compare)
# }
