#' Compute mean expression at group, compound, dose, and time levels
#'
#' @description
#' `mean_expression()` calculates average gene (or probe) expression at multiple hierarchical
#' levels of the experimental design: group, compound, dose, and time. The function derives the
#' study structure from `metadata` and uses summary matrices (from `get_matrix()`) to compute
#' level-wise means by matrix multiplication.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`), defining
#'   the grouping of compounds used to compute group-level means.
#' @param ge_matrix Gene expression matrix with `probes` in rows and samples (barcodes) in columns.
#' @param metadata A data frame or tibble of sample metadata corresponding to columns of `ge_matrix`.
#'   Must include `barcode`, `compound_name`, `dose_level`, and `time_level`.
#' @param probes Optional vector of probe/gene identifiers to subset `ge_matrix` prior to computing
#'   means. If `NULL` (default), all rows are used.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @return A named list with four matrices:
#' \itemize{
#'   \item `group_data`: mean expression at the group level (columns correspond to groups in `...`)
#'   \item `compound_data`: mean expression at the compound level
#'   \item `dose_data`: mean expression at the compound-dose level
#'   \item `time_data`: mean expression at the compound-dose-time level
#' }
#' Rows correspond to probes/genes (rows of `ge_matrix`), and columns correspond to the relevant
#' level-specific labels.
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' avg_expr <- mean_expression(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
#'
#' @seealso
#' [data_str()] and [get_matrix()] for constructing the design summary matrices used internally.
#'
#' @export
mean_expression <- function(...,
                            ge_matrix,
                            metadata,
                            probes = NULL,
                            error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  comps_gr <- test_group(...)
  lev_str <- data_str(comps_gr, metadata = metadata)
  design_mat <- get_matrix(lev_str)

  if (is.null(probes)) {
    y <- ge_matrix
  } else{
    true_probes <- probes %in% rownames(ge_matrix)
    if (any(!true_probes)) {
      false_probes <- probes[!true_probes]
      n_false <- length(false_probes)
      cli::cli_alert_warning(c("Given {.emph {n_false}} probe{?s} ",
                               "{style_bold(col_red(backtick(false_probes)))} ",
                               "{?is/are} not found in the gene expression data."),
                             wrap = TRUE)
    }
    y <- ge_matrix[probes[true_probes],]
  }
  estYi <- y %*% design_mat$group_mat
  colnames(estYi) <- names(comps_gr)
  estYij <- y %*% design_mat$compound_mat
  colnames(estYij) <- set_names("compound", metadata, "compound_name")
  estYijk <- y %*% design_mat$dose_mat
  colnames(estYijk) <- set_names("compound-dose", metadata, "compound_name")
  estYijkl <- y %*% design_mat$time_mat
  colnames(estYijkl) <- set_names("compound-dose-time", metadata, "compound_name")
  return(list(group_data = estYi,
              compound_data = estYij,
              dose_data = estYijk,
              time_data = estYijkl))
}

#' Extract mean expression for a single probe across design levels
#'
#' @description
#' `get_avgFC()` extracts the mean expression profile of a single probe (or gene) across the
#' hierarchical levels of the experiment defined in `metadata`. It summarizes expression at the
#' group, compound, compound-dose, and compound-dose-time levels (via `mean_expression()`), and
#' also returns the per-sample expression values for the same probe.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`), defining
#'   the grouping of compounds used to compute group-level means.
#' @param ge_matrix Gene expression matrix with `probes` in rows and samples (barcodes) in columns.
#' @param metadata A data frame or tibble of sample metadata corresponding to columns of `ge_matrix`.
#'   Must include `barcode`, `compound_name`, `dose_level`, and `time_level`.
#' @param probe_id A single probe/gene identifier. Must match one entry in `rownames(ge_matrix)`.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @return A named list of tibbles:
#' \itemize{
#'   \item `group_df`: mean expression per group (columns: `group`, `avg_x`)
#'   \item `compound_df`: mean expression per compound (columns: `group`, `compound`, `avg_x`)
#'   \item `dose_df`: mean expression per compound-dose (columns: `group`, `compound`, `dose`, `avg_x`)
#'   \item `time_df`: mean expression per compound-dose-time (columns: `group`, `compound`, `dose`, `time`, `avg_x`)
#'   \item `sample_df`: per-sample expression values (columns: `group`, `compound`, `dose`, `time`, `avg_x`)
#' }
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' avg_value <- get_avgFC(gr,
#'                        ge_matrix = sim_data$expression,
#'                        metadata = sim_data$metadata,
#'                        probe_id = "DE1")
#'
#' @seealso
#' [mean_expression()] for computing level-wise mean expression matrices used internally.
#'
#' @export
get_avgFC <- function(...,
                      ge_matrix,
                      metadata,
                      probe_id,
                      error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  test_input(probe_id, rownames(ge_matrix))
  comps_gr <- test_group(...)
  compounds <- unlist(comps_gr)
  avg_expr <- mean_expression(comps_gr,
                              ge_matrix = ge_matrix,
                              metadata = metadata,
                              probes = probe_id,
                              error_call = error_call)

  gr_tbl <- tibble::tibble(group = factor(names(comps_gr), levels = names(comps_gr)),
                   avg_x = as.vector(avg_expr$group_data))
  idx <-lapply(comps_gr, function(x) which(compounds %in% x))
  grp <- vector("character", length = length(compounds))
  for (i in 1:length(idx)) grp[idx[[i]]] <- names(comps_gr[i])
  com_tbl <- tibble::tibble(group = factor(grp, levels = names(idx)),
                    compound = factor(compounds, levels = compounds),
                    avg_x = as.vector(avg_expr$compound_data))

  dose_df <- metadata %>%
    dplyr::select(c(compound_name, dose_level)) %>%
    dplyr::distinct(compound_name, dose_level, .keep_all = TRUE)

  idx1 <-lapply(comps_gr, function(x) which(dose_df$compound_name %in% x))
  grp1 <- vector("character", length = nrow(dose_df))
  for (i in 1:length(idx1)) grp1[idx1[[i]]] <- names(comps_gr[i])
  dose_tbl <- tibble::tibble(
    group = factor(grp1, levels = names(idx1)),
    compound = factor(dose_df$compound_name, levels = compounds),
    dose = factor(dose_df$dose_level, levels = unique(dose_df$dose_level)),
    avg_x = as.vector(avg_expr$dose_data)
  )

  time_df <- metadata %>%
    dplyr::select(c(compound_name, dose_level, time_level)) %>%
    dplyr::distinct(compound_name, dose_level, time_level, .keep_all = TRUE)

  idx2 <- lapply(comps_gr, function(x) which(time_df$compound_name %in% x))
  grp2 <- vector("character", length = nrow(time_df))
  for (i in 1:length(idx2)) grp2[idx2[[i]]] <- names(comps_gr[i])
  time_tbl <- tibble::tibble(
    group = factor(grp2, levels = names(idx2)),
    compound = factor(time_df$compound_name, levels = compounds),
    dose = factor(time_df$dose_level, levels = unique(dose_df$dose_level)),
    time = factor(time_df$time_level, levels = unique(time_df$time_level)),
    avg_x = as.vector(avg_expr$time_data)
  )

  idx3 <- lapply(comps_gr, function(x) which(metadata$compound_name %in% x))
  grp3 <- vector("character", length = nrow(metadata))
  for (i in 1:length(idx3)) grp3[idx3[[i]]] <- names(comps_gr[i])
  sample_tbl <- tibble::tibble(
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

#' Subset perturbed gene expression data and matching metadata
#'
#' @description
#' `get_subset()` filters a perturbed gene expression matrix (`ge_matrix`) and its corresponding
#' sample `metadata` using selected dose levels and/or time points. The function also assigns a
#' `group` label to each sample based on the compound groups provided in `...`.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`), defining
#'   which compounds belong to each group.
#' @param ge_matrix Gene expression matrix with `probes` in rows and samples (barcodes) in columns.
#' @param metadata A data frame or tibble of sample metadata corresponding to columns of `ge_matrix`.
#'   Must include `barcode`, `compound_name`, `dose_level`, and `time_level`.
#' @param probes Optional vector of probe/gene identifiers to subset rows of `ge_matrix`. If `NULL`,
#'   all rows are kept.
#' @param dose Optional character vector of dose levels used to subset `metadata` and `ge_matrix`.
#'   If `NULL`, all dose levels are retained.
#' @param time Optional character vector of time levels used to subset `metadata` and `ge_matrix`.
#'   If `NULL`, all time levels are retained.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @details
#' The function validates `dose` and `time` against the available levels in `metadata`. It then
#' subsets `metadata` and retains only the matching columns in `ge_matrix`. A `group` column is
#' added to the returned `metadata` based on the compound group definitions supplied in `...`.
#'
#' @return A list with two components:
#' \itemize{
#'   \item `expression`: a subset of `ge_matrix` containing the selected samples (and optionally probes)
#'   \item `metadata`: the corresponding subset of `metadata` with an added `group` column
#' }
#'
#' @examples
#' sim_data <- simulate_tgxdata()
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sub_data <- get_subset(gr,
#'                        ge_matrix = sim_data$expression,
#'                        metadata = sim_data$metadata,
#'                        probes = paste0("DE", 1:10),
#'                        dose = "Dose1",
#'                        time = "Time1")
#'
#' @export
get_subset <- function(...,
                       ge_matrix,
                       metadata,
                       probes = NULL,
                       dose = NULL,
                       time = NULL,
                       error_call = rlang::caller_env()) {
  if (!is.null(dose) & !all(dose %in% unique(metadata$dose_level))) {
    cli::cli_abort(c("The input dose level is incorrect.",
                "x" = "Dose level ({style_bold(col_red(backtick(dose)))}) not match with the data.",
                "i" = "Please choose dose level either {add_or(style_bold(col_green(backtick(unique(metadata$dose_level)))))} instead.")
              , call = error_call)
  }
  if (!is.null(time) & !any(time %in% unique(metadata$time_level))) {
    cli::cli_abort(c("The input time level is incorrect.",
                "x" = "Time level ({style_bold(col_red(backtick(time)))}) not match with the data.",
                "i" = "Please choose time level either {add_or(style_bold(col_green(backtick(unique(metadata$time_level)))))} instead.")
              , call = error_call)
  }
  test_data(ge_matrix, metadata)
  comps_gr <- test_group(...)
  compounds <- unlist(comps_gr)
  if (is.null(time)) {
    expr_sp <- ge_matrix
  } else {
    metadata <- metadata %>%
      dplyr::filter(time_level %in% {{time}})
    expr_sp <- ge_matrix[, colnames(ge_matrix) %in% metadata$barcode]
  }
  if (is.null(dose)) {
    expr_sp <- expr_sp
  } else {
    metadata <- metadata %>%
      dplyr::filter(dose_level %in% {{dose}})
    expr_sp <- expr_sp[, colnames(expr_sp) %in% metadata$barcode]
  }
  tgx_class <- vector(length = nrow(metadata))
  for(i in 1:length(comps_gr)) {
    tgx_class[metadata$compound_name %in% comps_gr[[i]]] <- names(comps_gr)[i]
  }
  tgx_class <- factor(tgx_class, levels = names(comps_gr))
  metadata <- metadata %>%
    dplyr::mutate(group = tgx_class, .after = barcode)
  return(list(expression = expr_sp, metadata = metadata))
}



#' Compute mean expression for compound–dose–time subsets
#'
#' @description
#' `mean_subset()` subsets perturbed gene expression data using `dose` and/or `time` (via
#' `get_subset()`), then averages expression across replicates for each unique
#' compound–dose–time combination. The resulting expression matrix contains one column per
#' subset, and the returned metadata records the corresponding design labels and group
#' membership.
#'
#' @inheritParams get_subset
#' @param names_format Character string specifying how subset columns are labeled in the output.
#'   Passed to `set_names()` to construct `sample_id` values (e.g., `"compound-dose-time"`).
#'   Default is `"compound-dose-time"`.
#'
#' @return A list with two components:
#' \itemize{
#'   \item `expression`: a matrix of mean expression values, with probes/genes in rows and
#'     averaged subsets in columns (one column per compound–dose–time combination)
#'   \item `metadata`: a data frame describing the averaged subsets, including `compound_name`,
#'     `dose_level`, `time_level`, the generated `sample_id`, and the assigned `group`
#' }
#'
#' @details
#' The function groups the subsetted metadata by `compound_name`, `dose_level`, and `time_level`,
#' collects the corresponding sample barcodes, and computes row means over those samples in
#' `ge_matrix`. Group labels are assigned from the compound groups provided in `...`.
#'
#' @examples
#' sim_data <- simulate_tgxdata()
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sub_avg <- mean_subset(gr,
#'                        ge_matrix = sim_data$expression,
#'                        metadata = sim_data$metadata,
#'                        probes = paste0("DE", 1:10),
#'                        dose = "Dose1",
#'                        time = "Time1")
#'
#' @seealso
#' [get_subset()] for subsetting samples by dose/time and assigning group labels.
#'
#' @export
mean_subset <- function(...,
                        ge_matrix,
                        metadata,
                        probes = NULL,
                        dose = NULL,
                        time = NULL,
                        names_format = "compound-dose-time",
                        error_call = caller_env()) {
  comps_gr <- test_group(...)
  space_data <- get_subset(comps_gr,
                           ge_matrix = ge_matrix,
                           metadata = metadata,
                           probes = probes,
                           dose = dose,
                           time = time,
                           error_call = error_call)
  space_nest <- space_data$metadata %>%
    dplyr::group_by(compound_name, dose_level, time_level) %>%
    tidyr::nest()
  space_barcd <- lapply(space_nest$data, function(x) x$barcode)
  barcd_expr <- lapply(space_barcd, function(x)
    rowMeans(space_data$expression[,match(x, colnames(space_data$expression))]))
  space_expr <- do.call(cbind, barcd_expr)
  lab <- set_names(names_format, space_data$metadata, "compound_name")
  space_attr <- space_nest %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_id = lab)
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
