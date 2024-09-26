
#' Average gene expression in different levels
#'
#' @description
#'
#' The function `mean_expression()` is used to get average gene expression value in `group`,
#' `compound`, `dose` and `time` levels.
#'
#' @inheritParams update_data
#'
#' @return
#' A list of average data
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' avg_expr <- mean_expression(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
mean_expression <- function(...,
                            ge_matrix,
                            metadata,
                            probes = NULL,
                            multicore = FALSE,
                            store = FALSE,
                            output_dir = missing_arg(),
                            error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  comps_gr <- test_group(...)
  output_dir <- destination(output_dir)
  compounds <- unlist(comps_gr)
  ck_data <- update_data(
    compounds,
    ge_matrix = ge_matrix,
    metadata = metadata,
    probes = probes,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  ge_matrix <- ck_data$expression
  comp_dict <- ck_data$metadata
  lev_str <- data_str(comps_gr, metadata = comp_dict)
  design_mat <- get_matrix(lev_str)
  y <- ge_matrix # set apply() for multiple genes
  #  gene <- probes2genes(probe_id, organism = "rat")
  estYi <- y %*% design_mat$group_mat
  colnames(estYi) <- names(comps_gr)
  estYij <- y %*% design_mat$compound_mat
  colnames(estYij) <- set_names("compound", comp_dict, "compound_name")
  estYijk <- y %*% design_mat$dose_mat
  colnames(estYijk) <- set_names("compound-dose", comp_dict, "compound_name")
  estYijkl <- y %*% design_mat$time_mat
  colnames(estYijkl) <- set_names("compound-dose-time", comp_dict, "compound_name")
  # group_data <- estYi %>%
  #   tibble::as_tibble()
  # compound_data <- estYij %>%
  #   tibble::as_tibble()
  # dose_data <- estYijk %>%
  #   tibble::as_tibble()
  # time_data <- estYijkl %>%
  #   tibble::as_tibble()
  return(list(group_data = estYi,
              compound_data = estYij,
              dose_data = estYijk,
              time_data = estYijkl))
}

#' Average gene expression in different levels
#'
#' @description
#'
#' The function `get_avgFC()` is used to get average gene expression value in `group`,
#' `compound`, `dose` and `time` stages.
#'
#' @inheritParams update_data
#' @param probe_id The gene/probe ID.
#'
#' @return
#' A list of average data
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' avg_value <- get_avgFC(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata, probe_id = "DE1")
get_avgFC <- function(...,
                      ge_matrix,
                      metadata,
                      probe_id,
                      multicore = FALSE,
                      store = FALSE,
                      output_dir = missing_arg(),
                      error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  test_input(probe_id, rownames(ge_matrix))
  comps_gr <- test_group(...)
  output_dir <- destination(output_dir)
  compounds <- unlist(comps_gr)
  ck_data <- update_data(
    compounds,
    ge_matrix = ge_matrix,
    metadata = metadata,
    probes = NULL,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  ge_matrix <- ck_data$expression
  metadata <- ck_data$metadata
  avg_expr <- mean_expression(comps_gr,
                              ge_matrix = ge_matrix,
                              metadata = metadata,
                              probes = probe_id,
                              multicore = multicore,
                              store = store,
                              output_dir = output_dir,
                              error_call = error_call)

  gr_tbl <- tibble(group = factor(names(comps_gr), levels = names(comps_gr)),
                   avg_x = as.vector(avg_expr$group_data))
  idx <-lapply(comps_gr, function(x) which(compounds %in% x))
  grp <- vector("character", length = length(compounds))
  for (i in 1:length(idx)) grp[idx[[i]]] <- names(comps_gr[i])
  com_tbl <- tibble(group = factor(grp, levels = names(idx)),
                    compound = factor(compounds, levels = compounds),
                    avg_x = as.vector(avg_expr$compound_data))

  dose_df <- metadata %>%
    select(c(compound_name, dose_level)) %>%
    distinct(compound_name, dose_level, .keep_all = TRUE)

  idx1 <-lapply(comps_gr, function(x) which(dose_df$compound_name %in% x))
  grp1 <- vector("character", length = nrow(dose_df))
  for (i in 1:length(idx1)) grp1[idx1[[i]]] <- names(comps_gr[i])
  dose_tbl <- tibble(
    group = factor(grp1, levels = names(idx1)),
    compound = factor(dose_df$compound_name, levels = compounds),
    dose = factor(dose_df$dose_level, levels = unique(dose_df$dose_level)),
    avg_x = as.vector(avg_expr$dose_data)
  )

  time_df <- metadata %>%
    select(c(compound_name, dose_level, time_level)) %>%
    distinct(compound_name, dose_level, time_level, .keep_all = TRUE)

  idx2 <- lapply(comps_gr, function(x) which(time_df$compound_name %in% x))
  grp2 <- vector("character", length = nrow(time_df))
  for (i in 1:length(idx2)) grp2[idx2[[i]]] <- names(comps_gr[i])
  time_tbl <- tibble(
    group = factor(grp2, levels = names(idx2)),
    compound = factor(time_df$compound_name, levels = compounds),
    dose = factor(time_df$dose_level, levels = unique(dose_df$dose_level)),
    time = factor(time_df$time_level, levels = unique(time_df$time_level)),
    avg_x = as.vector(avg_expr$time_data)
  )

  idx3 <- lapply(comps_gr, function(x) which(metadata$compound_name %in% x))
  grp3 <- vector("character", length = nrow(metadata))
  for (i in 1:length(idx3)) grp3[idx3[[i]]] <- names(comps_gr[i])
  sample_tbl <- tibble(
    group = factor(grp3, levels = names(idx3)),
    compound = factor(metadata$compound_name, levels = compounds),
    dose = factor(metadata$dose_level, levels = unique(metadata$dose_level)),
    time = factor(metadata$time_level, levels = unique(metadata$time_level)),
    avg_x = as.vector(ge_matrix[probe_id,])
  )
  return(list(
    group_df = gr_tbl,
    compound_df = com_tbl,
    dose_df = dose_tbl,
    time_df = time_tbl,
    sample_df = sample_tbl
  ))
}



#' Get a subset of gene expression data based on selected parameters
#'
#' @description
#'
#' This function subsets gene expression data (`ge_matrix`) and metadata (`metadata`) based on selected
#' dose levels, time points, and other parameters. It performs error handling for incorrect dose
#' or time inputs and offers options for parallel processing and data storage.
#'
#' @inheritParams update_data
#' @param dose A character vector of dose levels to subset the metadata and expression matrix. If NULL, all dose levels will be used.
#' @param time A character vector of time points to subset the metadata and expression matrix. If NULL, all time points will be used.
#'
#' @details
#' The function first checks for valid dose and time levels in the `metadata` and throws an error if
#' any mismatch is found. It then subsets the gene expression data (`ge_matrix`) and metadata
#' (`metadata`) based on the specified dose and time levels. Parallel processing is optionally supported via
#' the `multicore` parameter, and the results can be stored to a specified directory if `store` is TRUE.
#'
#' @return A list with the following components:
#' \item{expression}{A matrix of gene expression data that has been subset based on the specified dose and time levels, and optionally the selected probes.}
#' \item{metadata}{A data frame containing metadata corresponding to the samples in the subsetted expression matrix, including dose, time, and compound information.}
#'
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata()
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sub_data <- get_subset(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata, probes = paste0("DE", 1:10), dose = "Dose1", time = "Time1")
get_subset <- function(...,
                       ge_matrix,
                       metadata,
                       probes = NULL,
                       dose = NULL,
                       time = NULL,
                       multicore = FALSE,
                       store = FALSE,
                       output_dir = missing_arg(),
                       error_call = caller_env()) {
  if (!is.null(dose) & !all(dose %in% unique(metadata$dose_level))) {
    cli_abort(c("The input dose level is incorrect.",
                "x" = "Dose level ({style_bold(col_red(backtick(dose)))}) not match with the data.",
                "i" = "Please choose dose level either {add_or(style_bold(col_green(backtick(unique(metadata$dose_level)))))} instead.")
              , call = error_call)
  }
  if (!is.null(time) & !any(time %in% unique(metadata$time_level))) {
    cli_abort(c("The input time level is incorrect.",
                "x" = "Time level ({style_bold(col_red(backtick(time)))}) not match with the data.",
                "i" = "Please choose time level either {add_or(style_bold(col_green(backtick(unique(metadata$time_level)))))} instead.")
              , call = error_call)
  }
  test_data(ge_matrix, metadata)
  comps_gr <- test_group(...)
  output_dir <- destination(output_dir)
  compounds <- unlist(comps_gr)
  ck_data <- update_data(
    compounds,
    ge_matrix = ge_matrix,
    metadata = metadata,
    probes = probes,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  ge_matrix <- ck_data$expression
  comp_dict <- ck_data$metadata
  #
  #   if (is.null(probes)) {
  #     expr_data <- expr_data
  #   } else{
  #     true_probes <- probes %in% rownames(expr_data)
  #     if (any(!true_probes)) {
  #       false_probes <- probes[!true_probes]
  #       n_false <- length(false_probes)
  #       cli_alert_warning(c("Given {.emph {n_false}} probe{?s} ",
  #                           "{style_bold(col_red(backtick(false_probes)))} ",
  #                           "{?is/are} not found in the gene expression data."),
  #                         wrap = TRUE)
  #     }
  #     expr_data <- expr_data[probes[true_probes],]
  #   }
  # lev_str <- data_str(comps_gr, metadata = comp_dict)
  # design_mat <- design_matrix(lev_str)

  if (is.null(time)) {
    expr_sp <- ge_matrix
  } else {
    comp_dict <- comp_dict %>%
      dplyr::filter(time_level %in% {{time}})
    expr_sp <- ge_matrix[, colnames(ge_matrix) %in% comp_dict$barcode]
  }
  if (is.null(dose)) {
    expr_sp <- expr_sp
  } else {
    comp_dict <- comp_dict %>%
      dplyr::filter(dose_level %in% {{dose}})
    expr_sp <- expr_sp[, colnames(expr_sp) %in% comp_dict$barcode]
  }

  #  For average data
  # if (identical(space, "dose")) {
  #   if (average) {
  #     expr_sp <- Y %*% design_mat$qGama
  #     #colnames(expr_sp) <- level_names(..., ..., ...)
  #   }
  #   dose_dict <- comp_dict %>%
  #     dplyr::filter(dose_level == "High")
  #   expr_sp <- Y[, colnames(Y) %in% dose_dict$barcode]
  #  dim(expr_sp)
  # }
  #
  # if (identical(space, "time")) {
  #   if (average) {
  #     expr_sp <- Y %*% design_mat$qDelta
  #     #colnames(expr_sp) <- level_names(..., ..., ...)
  #   }
  #   time_dict <- comp_dict %>%
  #     dplyr::filter(time_level == "24 hr")
  #   expr_sp <- Y[, colnames(Y) %in% time_dict$barcode]
  #   dim(expr_sp)
  # }
  tgx_class <- vector(length = nrow(comp_dict))
  for(i in 1:length(comps_gr)) {
    tgx_class[comp_dict$compound_name %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  comp_dict <- comp_dict %>%
    dplyr::mutate(group = tgx_class, .after = barcode)
  return(list(expression = expr_sp, metadata = comp_dict))
}


#' Compute Mean Expression of Subsets
#'
#' This function calculates the mean expression of specific subsets from a gene expression matrix based on metadata criteria such as dose and time.
#'
#' @inheritParams get_subset
#' @param names_format A string specifying the format for labeling samples in the output. Defaults to "compound-dose-time".
#'
#' @return A list containing two elements:
#' \item{expression}{A matrix of mean gene expression values for each selected subset.}
#' \item{metadata}{A dataframe of metadata for the resulting subsets, including sample IDs and group information.}
#'
#' @details This function first subsets the gene expression matrix based on the selected compounds, dose levels, and time points. It then computes the mean expression for each subset. The output includes both the averaged gene expression matrix and the corresponding metadata for each subset.
#'
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata()
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sub_avg <- mean_subset(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata, probes = paste0("DE", 1:10), dose = "Dose1", time = "Time1")
mean_subset <- function(...,
                        ge_matrix,
                        metadata,
                        probes = NULL,
                        dose = NULL,
                        time = NULL,
                        names_format = "compound-dose-time",
                        multicore = FALSE,
                        store = FALSE,
                        output_dir = missing_arg(),
                        error_call = caller_env()) {
  comps_gr <- test_group(...)
  space_data <- get_subset(comps_gr,
                           ge_matrix = ge_matrix,
                           metadata = metadata,
                           probes = probes,
                           dose = dose,
                           time = time,
                           multicore = multicore,
                           store = store,
                           output_dir = output_dir,
                           error_call = error_call)
  space_nest <- space_data$metadata %>%
    dplyr::group_by(compound_name, dose_level,time_level) %>%
    nest()
  space_barcd <- lapply(space_nest$data, function(x) x$barcode)
  barcd_expr <- lapply(space_barcd, function(x)
    rowMeans(space_data$expression[,match(x, colnames(space_data$expression))]))
  space_expr <- do.call(cbind, barcd_expr)
  lab <- set_names(names_format, space_data$metadata, "compound_name")
  space_attr <- space_nest %>%
    select(-data) %>%
    ungroup() %>%
    mutate(sample_id = lab)
  tgx_class <- vector(length = nrow(space_attr))
  for(i in 1:length(comps_gr)) {
    tgx_class[space_attr$compound_name %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  space_attr <- space_attr %>%
    dplyr::mutate(group = tgx_class, .after = sample_id)
  colnames(space_expr) <- lab
  return(list(expression = space_expr, metadata = space_attr))
}

