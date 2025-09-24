#
# gaogene <- read.table("E:/TGGATES_v0.2/GeneList/Gao2010all.txt")
# kswgene4 <- read.table("E:/TGGATES_v0.2/GeneList/Kiyosawa2006.txt")
# kswgene7 <- read.table("E:/TGGATES_v0.2/GeneList/Kiyosawa2007.txt")
# yamgene11 <- read.table("E:/TGGATES_v0.2/GeneList/Kiyosawa2011.txt")
#
# s1 <- probes2genes(unique(gaogene$V1))
# s2 <- probes2genes(unique(kswgene4$V1))
# s3 <- probes2genes(unique(kswgene7$V1))
# s4 <- probes2genes(unique(yamgene11$V1))
#
# g1 <- glue::glue("{s1$PROBEID} ({s1$SYMBOL})")
# g2 <- glue::glue("{s2$PROBEID} ({s2$SYMBOL})")
# g3 <- glue::glue("{s3$PROBEID} ({s3$SYMBOL})")
# g4 <- glue::glue("{s4$PROBEID} ({s4$SYMBOL})")
#
#
# max_length <- max(length(g1), length(g2), length(g3), length(g4))
# set_a <- c(g1, rep(NA, max_length - length(g1)))
# set_b <- c(g2, rep(NA, max_length - length(g2)))
# set_c <- c(g3, rep(NA, max_length - length(g3)))
# set_d <- c(g4, rep(NA, max_length - length(g4)))
#
# gene_list <- tibble::tibble(Set_A = set_b, Set_B = set_c, Set_C = set_a, Set_D = set_d)
#
# write.csv(gene_list, "D:/ToxAssay/tables/S1.csv")
#
# system.time(
#   res3 <- tgx_degs(comps_gr,
#                    ge_matrix = expr_data,
#                    metadata = attr_data,
#                    p_cutoff = 0.05,
#                    p_adjust = "bonferroni",
#                    gr_diff = TRUE,
#                    multicore = FALSE,
#                    log10p =  TRUE,
#                    output_dir = missing_arg(),
#                    error_call = caller_env())
# )
# table(res3$sig_type)
# res4 <- res3 %>%
#   filter(sig_type == "DE")
#
# write.csv(res4, "D:/ToxAssay/tables/S2.csv")
#
#
# dise_voc <- read_csv("D:/ToxAssay/ToxAssay_genes/liver_diseases.csv")
#
# ctd <- get_ctd(compounds = comps_gr$Positive)
# trans_data <- get_transaction(compound_gene = ctd$compound_gene,
#                               compound_disease = ctd$compound_disease,
#                               compounds = NULL,
#                               genes = res4$gene_symbol,
#                               diseases = dise_voc$DiseaseID)
# aop_data <- get_aops(transaction = trans_data[[1]],
#                      genes = trans_data[[2]],
#                      diseases = trans_data[[3]],
#                      ci_metric = "lift")
# aop_res <- aop_data %>%
#   filter(rhs %in% dise_voc$DiseaseID) %>%
#   mutate(Disease = dise_voc$DiseaseName[match(rhs, dise_voc$DiseaseID)],
#          Group = dise_voc$DiseaseGroup[match(rhs, dise_voc$DiseaseID)], .before = support) %>%
#   dplyr::mutate(lhs = block_fst(lhs),
#                 rhs = gsub("MESH:", "", rhs)) %>%
#   arrange(desc(lift))
#
#   write.csv(aop_res, "D:/ToxAssay/tables/S3.csv")
# #----------------------
# enrc <- get_enrichment(gene_df = res4,
#                             category = "WikiPathways",
#                             organism = "rat",
#                             path_n = 100,
#                             score_threshold = 200,
#                             version = "12")
#
# write.csv(enrc, "D:/ToxAssay/tables/S4d.csv")
#
# #-
# net_edgebet <- get_netdata(
#   res4,
#   organism = "rat",
#   score_threshold = 200,
#   cluster_method = "edge.betweenness",
#   version = "12")
#
#
# ppi_net <- net_edgebet$vertices %>%
#   dplyr::select(c("probe_id", "gene_symbol", "entrez_id","STRING_id", "gene_name",
#                   "degree", "betweenness", "closenes", "eigenes", "gene_class"))
#
# write.csv(ppi_net, "D:/ToxAssay/tables/S5b.csv")
#
# unique(net_edgebet$edges$from)
# net_edgebet$vertices$gene_symbol[!net_edgebet$vertices$gene_symbol %in% union(unique(net_edgebet$edges$to), unique(net_edgebet$edges$from))]
#
#
# #-----------------------
# library(tidymodels)
#
#
# class_pred <- function(...,
#                      ge_matrix,
#                      metadata,
#                      probe_list = NULL,
#                      dose = NULL,
#                      time = NULL,
#                      model = c("lr", "svm", "xgboost", "rf", "knn"),
#                      nfold = 10,
#                      nrep = 100,
#                      multicore = FALSE,
#                      store = FALSE,
#                      output_dir = missing_arg(),
#                      error_call = caller_env()) {
#   test_input(model, auto_input = TRUE)
#   comps_group <- test_group(...)
#   test_data(ge_matrix, metadata)
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
#   res_mat <- NULL
#   mod_fit <- vector(mode = "list", length = length(probe_list))
#   for (i in 1:length(probe_list)) {
#     dt <- get_subset(comps_gr,
#                      ge_matrix = ge_matrix,
#                      metadata = metadata,
#                      probes = probe_list[[i]],
#                      dose = dose,
#                      time = time,
#                      multicore = multicore,
#                      store = store,
#                      output_dir = output_dir,
#                      error_call = error_call)
#     expr_tbl <- tibble::as_tibble(t(dt$expression)) %>%
#       dplyr::mutate(group = as.factor(dt$metadata$group))
#     folds <-  vfold_cv(expr_tbl, v = nfold, repeats = nrep, strata = group)
#     mod_fit[[i]] <- mod_wf %>%
#       fit_resamples(resamples = folds,
#                     metrics = metric_set(recall, precision, f_meas, accuracy, kap,roc_auc, sens, spec),
#                     control = control_resamples(save_pred = TRUE, verbose = TRUE)) %>%
#       suppressMessages()
#     res <- glue::glue("{round(collect_metrics(mod_fit[[i]])$mean,3)} ({round(collect_metrics(mod_fit[[i]])$std_err,3)})")
#     res_mat <- cbind(res_mat, res)
#   }
#   #
#   # # Evaluation
#   # acc <- data.frame(
#   #   RF = glue("{round(collect_metrics(rf_fit)$mean,3)} ({round(collect_metrics(rf_fit)$std_err,3)})"),
#   #   KNN = glue("{round(collect_metrics(log_fit)$mean,3)} ({round(collect_metrics(log_fit)$std_err,3)})"),
#   #   SVM = glue("{round(collect_metrics(knn_fit)$mean,3)} ({round(collect_metrics(knn_fit)$std_err,3)})"),
#   #   LR = glue("{round(collect_metrics(svm_fit)$mean,3)} ({round(collect_metrics(svm_fit)$std_err,3)})")
#   # )
#   rownames(res_mat) <- c("Recall", "Precision", "F_measure", "Accuracy",
#                      "Kappa", "AUC", "Sensitivity", "Specificity")
#   colnames(res_mat) <- names(probe_list)
#   return(res_mat)
# }
#
#
# new_core <- pull(read.table("E:/TGGATES_v0.2/GeneList/new_core.txt", header = FALSE))
# ksw04_genes <- probes2genes(kswgene4$V1)$SYMBOL
# ksw07_genes <- probes2genes(kswgene7$V1)$SYMBOL
# gao10_genes <- probes2genes(gaogene$V1)$SYMBOL
# yam11_genes <- probes2genes(yamgene11$V1)$SYMBOL
# setA <- unique(res3$probe_id[res3$gene_symbol %in% ksw04_genes])
# setB <- unique(res3$probe_id[res3$gene_symbol %in% ksw07_genes])
# setC <- unique(res3$probe_id[res3$gene_symbol %in% gao10_genes])
# setD <- unique(res3$probe_id[res3$gene_symbol %in% yam11_genes])
# probe_list <- list(core = new_core, DEGs = res4$probe_id,
#                    setA = setA, setB = setB, setC = setC, setD = setD)
#
# bb <- class_pred(comps_gr, ge_matrix = expr_data, metadata = attr_data,
#                  probe_list = probe_list, nrep = 100)
#
# write.csv(bb, "D:/ToxAssay/tables/S10.csv")

