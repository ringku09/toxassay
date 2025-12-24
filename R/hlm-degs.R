#' F-test for comparing group means in hierarchical expression data
#'
#' @description
#' `fh()` performs an F-test for differences in group means using a hierarchical linear model (HLM)
#' representation of the experimental design. Given an expression matrix `Y` and quadratic matrices
#' `A` and `B` derived from the study structure, the function computes an F-statistic for each row
#' (e.g., probe/gene) and returns the corresponding p-values.
#'
#' @param Y Gene expression data matrix with `probes` in rows and samples in columns.
#' @param A Quadratic matrix A (typically from `get_matrix()`).
#' @param B Quadratic matrix B (typically from `get_matrix()`).
#' @param a Integer specifying the number of compound groups. This value can be obtained
#'   from the output of `data_str()` (i.e., `data_str(...)[[1]]`).
#' @param p Integer specifying the number of clusters in the experiment (model degrees-of-freedom component). This value
#'   can be obtained from `data_str()` (e.g., `sum(data_str(...)[[4]])`).
#' @param n Integer specifying the total number of samples in the experiment. This value
#'   can be obtained from `data_str()` (e.g., `sum(data_str(...)[[5]])`).
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return A named numeric vector of p-values, one per row of `Y`.
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sim_str <- data_str(gr, metadata = sim_data$metadata)
#' summary_mat <- get_matrix(sim_str)
#' fh(Y = sim_data$expression,
#'    A = summary_mat$A,
#'    B = summary_mat$B,
#'    a = sim_str[[1]],
#'    p = sum(sim_str[[4]]),
#'    n = sum(sim_str[[5]]))
#'
#' @seealso
#' [data_str()] for deriving `a`, `p`, and `n`,
#' [get_matrix()] for constructing quadratic matrices `A` and `B`.
#'
#' @export
fh <- function(Y, A, B, a, p, n, error_call = rlang::caller_env()) {
  if (!(inherits(Y, "matrix"))) {
    cli::cli_abort(c("The {.arg Y} must be a matrix.",
             "i" = "Please make {.arg Y} as matrix."), call = error_call)
  }
  if (! ncol(Y) == n) {
    cli::cli_abort(c("The matrix {.arg Y} must have {style_bold(col_red(backtick(n)))} columns for given {.arg expr_str}.",
                "x" = "The number of columns in {.arg Y} ({ncol(Y)}) is not equal to the \\
                total sample size of the experiment ({n}).",
                "i" = "Please provide matrix {.arg Y} with number of columns {style_bold(col_red(backtick(n)))} \\
                or provide {.arg expr_str} appropriately."), call = error_call)
  }
  f_pi <- get_pi(Y, A)
  BB <- methods::as(B, 'dgCMatrix')
  f_psi <- get_psi(Y, BB)
  eta <- (n - p)/(a - 1)
  F_h <- eta* f_pi / f_psi
  pval <- stats::pf(F_h, a - 1, n - p, lower.tail = FALSE)
  names(pval) <- rownames(Y)
  return(pval)
}

#' Identify differentially expressed genes using an HLM F-test
#'
#' @description
#' `de_genes()` detects differentially expressed (DE) genes between two or more compound groups
#' using an F-test derived from a hierarchical linear model (HLM). The experimental design is
#' inferred from `metadata`, and the required quadratic and summary matrices are constructed from
#' the design structure. Genes are classified as `"DE"` when their (optionally adjusted) p-values
#' are less than or equal to `p_cutoff`; otherwise they are classified as `"EE"`.
#'
#' @param ... Compound groups supplied as vectors or as a single list (see `test_group()`), defining
#'   the groups to be compared.
#' @param ge_matrix Gene expression matrix with `probes` in rows and samples (barcodes) in columns.
#' @param metadata A data frame or tibble of sample metadata corresponding to columns of `ge_matrix`.
#'   Must include `barcode`, `compound_name`, `dose_level`, and `time_level` (see `test_data()`).
#' @param p_cutoff Numeric p-value cutoff used to classify probes/genes as `"DE"` versus `"EE"`.
#'   Default is `0.05`.
#' @param p_adjust Character string specifying the p-value adjustment method passed to
#'   [stats::p.adjust()]. Default is `"none"`.
#' @param gr_diff Logical indicating whether to include group-level mean expression columns in the
#'   output. Default is `TRUE`.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @return A tibble with at least the columns:
#' \itemize{
#'   \item `probe_id`: probe/gene identifier (row names of `ge_matrix`)
#'   \item `p_value`: (adjusted) p-value from the HLM F-test
#'   \item `sig_type`: `"DE"` if `p_value <= p_cutoff`, otherwise `"EE"`
#' }
#' If `gr_diff = TRUE`, additional columns are included giving the estimated group-level mean
#' expression for each group.
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' de_genes(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
#'
#' @seealso
#' [stats::p.adjust()] for p-value adjustment methods.
#'
#' @export
de_genes <- function(...,
                     ge_matrix,
                     metadata,
                     p_cutoff = 0.05,
                     p_adjust = "none",
                     gr_diff = TRUE,
                     error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  com_group <- test_group(...)
  lev_str <- data_str(com_group, metadata = metadata, error_call = error_call)
  a <- lev_str[[1]]
  p <- sum(lev_str[[4]])
  n <-  sum(lev_str[[5]])
  quadmat <- get_matrix(lev_str)
  pval <- fh(Y = ge_matrix, A = quadmat$A, B = quadmat$B, a = a, p = p, n = n)
  pvalues <- stats::p.adjust(pval, method = p_adjust)
  sig_probe <- pvalues <= p_cutoff
  sig_exp <- signif(pvalues, digits = 3)
  sig_df <- tibble::tibble(probe_id = names(sig_exp), p_value = sig_exp) %>%
    dplyr::mutate(sig_type = ifelse(sig_probe,"DE", "EE"))
  if (gr_diff) {
    avg_gfc <- ge_matrix %*% quadmat$group_mat
    colnames(avg_gfc) <- names(com_group)
    rownames(avg_gfc) <- rownames(ge_matrix)
    sig_df <- col_diff(avg_gfc) %>%
      dplyr::left_join(sig_df, by = "probe_id")
    sig_df <- tibble::as_tibble(avg_gfc, rownames = "probe_id") %>%
      dplyr::left_join(sig_df, by = "probe_id")
  }
  return(sig_df)
}

#' Identify outlier-free reduced differentially expressed genes
#'
#' @description
#' `tox_degs()` identifies an outlier-robust set of differentially expressed genes (DEGs)
#' by applying a leave-\eqn{(m-1)}-out cross-validation framework across compounds. Let \eqn{C} be
#' the set of all compounds, partitioned into \eqn{a} groups, where group \eqn{i} contains \eqn{m_i}
#' compounds, i.e., \eqn{C = \cup_{i=1}^{a} C_i}. First, DEGs are identified from the full dataset
#' using the proposed HLM F-test (see `de_genes()`), yielding an initial set \eqn{G}. Because \eqn{G}
#' may include genes driven by chemical-specific outlying expression, the function constructs
#' multiple reduced DEG sets by leaving out \eqn{(m_i-1)} compounds at a time within the relevant
#' group-wise structure.
#'
#' Specifically, the function generates \eqn{M} DEG sets \eqn{G_k} by recomputing DEGs using a reduced
#' compound set \eqn{C_k = C \setminus \{c_{-k}\}} (for \eqn{k = 1, \dots, M}), where each reduced set
#' excludes a different subset of compounds (equivalently, retains a different held-in compound).
#' The final outlier-free reduced DEG set is defined as the intersection:
#' \deqn{G^\* = \cap_{k=1}^{M} G_k.}
#' Genes significant in the full model but not retained in \eqn{G^\*} are flagged as compound-driven.
#'
#' @inheritParams de_genes
#' @param log10p Logical indicating whether to add a `log10p` column computed as `-log10(p_value)`.
#'   Default is `TRUE`.
#'
#' @return
#' A tibble (classed as `"toxassay"`) containing results from the full HLM analysis together with
#' cross-validation filtering. The output includes the following key columns:
#' \itemize{
#'   \item `probe_id` and, when available, gene annotations (e.g., `gene_symbol`, `entrez_id`, `gene_name`)
#'   \item `p_value` (optionally adjusted using `p_adjust`)
#'   \item `sig_type`, where:
#'     \itemize{
#'       \item `"DE"` indicates genes retained in the outlier-free reduced set \eqn{G^\*}
#'       \item `"EE"` indicates non-significant genes
#'       \item `"CE"` indicates genes significant in the full model but not retained in \eqn{G^\*}
#'     }
#'   \item `log10p` if `log10p = TRUE`
#' }
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' tox_degs(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
#'
#' @seealso
#' [de_genes()] for the underlying HLM F-test used to identify DEGs.
#'
#' @export
tox_degs <- function(...,
                      ge_matrix,
                      metadata,
                      p_cutoff = 0.05,
                      p_adjust = "none",
                      gr_diff = TRUE,
                      log10p = TRUE,
                      error_call = rlang::caller_env()) {
  comps_group <- test_group(...)
  test_data(ge_matrix, metadata)
  compounds <- as.vector(unlist(comps_group))
  full_model <- de_genes(
    comps_group,
    ge_matrix = ge_matrix,
    metadata = metadata,
    p_cutoff = p_cutoff,
    p_adjust = p_adjust,
    gr_diff = gr_diff,
    error_call = error_call
  )
  if ("arr_design" %in% names(metadata)) {
    chip <- unique(metadata$arr_design)
    if(identical(chip, "Rat230_2")) {
      organism = "rat"
    }else if(identical(chip, "HG-U133_Plus_2")) {
      organism = "human"
    }
    all_genes <- probes2genes(full_model$probe_id, organism)
    hlm_tab <- full_model[match(all_genes$probe_id, full_model$probe_id), ]
    hlm_df <- hlm_tab %>%
      dplyr::mutate(entrez_id = all_genes$entrez_id,
                    gene_symbol = all_genes$gene_symbol,
                    gene_name = all_genes$gene_name,
                    .after = 1)
  } else {
    hlm_df <- full_model %>%
      dplyr::mutate(gene_symbol = probe_id, .after = 1)
  }
  hlm_df <- hlm_df %>%
    dplyr::distinct(gene_symbol, .keep_all= TRUE) %>%
    #distinct(probe_id, .keep_all= TRUE) %>% # calculate average of rep probes
    dplyr::arrange(p_value)
  sig_probes <- hlm_df$probe_id[hlm_df$sig_type == "DE"]
  expr_loodata <- ge_matrix[sig_probes, ]
  loo_probes <- vector(mode = "list", length = length(compounds))
  for (i in seq_len(length(compounds))) {
    comp_grnew <- lapply(comps_group, function(x) {
      if (compounds[i] %in% x) x[x == compounds[i]] else x })
    comp_loo <- metadata$compound_name %in% unlist(comp_grnew)
    attr_new <- metadata[comp_loo, ]
    expr_new <- expr_loodata[,match(metadata$barcode[comp_loo], colnames(expr_loodata))]
    loo_model <- de_genes(  # .parallel setup run each lapply , that take more time,
      comp_grnew,    #   need to setup once
      ge_matrix = expr_new,
      metadata = attr_new,
      p_cutoff = p_cutoff,
      p_adjust = p_adjust,
      gr_diff = FALSE,
      error_call = error_call
    )
    loo_probes[[i]] <- loo_model$probe_id[loo_model$sig_type == "DE"]
  }
  comon_probes <- Reduce(intersect, loo_probes)
  sig_tab <- hlm_df %>%
    dplyr::mutate(sig_type = ifelse(!(hlm_df$probe_id %in% comon_probes) & (hlm_df$sig_type == "DE"),
                                    "CE", sig_type))
  if (log10p) {
    sig_tab <- sig_tab %>%
       dplyr::mutate(log10p = -log10(p_value),.after = p_value)
  }
  n_probe <- nrow(ge_matrix)
  full_n <- sum(sig_tab$sig_type == "DE" | sig_tab$sig_type == "CE")
  loo_n <- sum(sig_tab$sig_type == "DE")
  unique_genes <- nrow(sig_tab)
  fil1 <- n_probe - full_n
  fil2 <- full_n - loo_n
    cli::cli_ul(c(
      cli::bg_magenta("Total number of unique genes = {unique_genes} ({n_probe} probes)"),
      cli::bg_red("Number of significant genes = {loo_n} (at {'\u03B1'}  = {p_cutoff})")
    ))
    return(structure(sig_tab, class = c("toxassay", "tbl_df", "tbl", "data.frame")))
}
