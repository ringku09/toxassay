#' Hierarchical linear model
#'
#' The `hlm()` function models gene expression data utilizing the hierarchical linear
#' model (HLM).
#'
#' @param Y Gene expression data matrix, probes in row and samples in column.
#' @param A quadratic matrix A.
#' @param B quadratic matrix B.
#' @param a Number of compound groups.
#' @param p Number of possible cluster in the data.
#' @param n Total number of samples in the data.
#'
#' @return
#' A vector of p-values with a size equal to the number of rows in **Y**.
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sim_str <- data_str(gr, metadata = sim_data$metadata)
#' summary_mat <- get_matrix(sim_str)
#' hlm(Y = sim_data$expression, A = summary_mat$A, B =summary_mat$B, a = sim_str[[1]], p = sum(sim_str[[4]]), n = sum(sim_str[[5]]))
hlm <- function(Y, A, B, a, p, n, error_call = rlang::caller_env()) {
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

#' Identification of initial differentially expressed genes
#'
#' The `de_genes()` function identify initial differentially expressed genes using the HLM.
#'
#' `de_genes()` identifies genes that exhibit statistically significant differential
#' expression across two or more query compound groups. A gene is deemed statistically
#' significant if the *p-value* is significant in HLM.
#'
#' @param gr_diff Show the column(s) for group mean differences. The default is `FALSE`.
#' @param p_adjust Adjusted `p-values` using Bonferroni methods. default is `TRUE`.
#' @param p_cutoff P-value cutoff threshold.
#' @inheritParams expand_data
#'
#' @return
#' An object of `data.frame` or `tibble` about genes table.
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' de_genes(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
de_genes <- function(...,
                     ge_matrix,
                     metadata,
                     p_cutoff = 0.05,
                     p_adjust = "none",
                     gr_diff = FALSE,
                     multicore = FALSE,
                     store = FALSE,
                     output_dir = missing_arg(),
                     error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  com_group <- test_group(...)
  output_dir <- destination(output_dir)
  compounds <- as.vector(unlist(com_group))
    ck_data <- update_data(
      compounds,
      ge_matrix = ge_matrix,
      metadata = metadata,
      multicore = multicore,
      store = store,
      output_dir = output_dir,
      error_call = error_call
    )
    ge_matrix <- ck_data$expression
    metadata <- ck_data$metadata
# ekhane eddit korte hobe, kintu ki seta mone nai
  lev_str <- data_str(com_group, metadata = metadata, error_call = error_call)
  a <- lev_str[[1]]
  p <- sum(lev_str[[4]])
  n <-  sum(lev_str[[5]])
  quadmat <- get_matrix(lev_str)
  pval <- hlm(Y = ge_matrix, A = quadmat$A, B = quadmat$B, a = a, p = p, n = n)
  pvalues <- p.adjust(pval, method = p_adjust)
  sig_probe <- pvalues <= p_cutoff
  sig_exp <- signif(pvalues, digits = 3)
  sig_df <- tibble::tibble(probe_id = names(sig_exp), p_value = sig_exp) %>%
    dplyr::mutate(sig_type = ifelse(sig_probe,"DE", "EE"))
  if (gr_diff) {
    avg_gfc <- ge_matrix %**% quadmat$group_mat
    colnames(avg_gfc) <- names(com_group)
    rownames(avg_gfc) <- rownames(ge_matrix)
    sig_df <- col_diff(avg_gfc) %>%
      dplyr::left_join(sig_df, by = "probe_id")
    sig_df <- tibble::as_tibble(avg_gfc, rownames = "probe_id") %>%
      dplyr::left_join(sig_df, by = "probe_id")
  }
  return(sig_df)
}

#' Identification of causal differentially expressed genes
#'
#' The function `tgx_genes()` is used to find final gene set using Leave One Out Compound (LOOC) method.
#' Genes misrepresented by a compound were removed using Leave One Out Compound (LOOC) method.
#' In LOOC process, each compound is removed from the experiment and find differentially expressed (DE)
#' genes for every LOO. The final gene set is therefore the intersection of the LOOC genes. Genes which
#' are co-regulated by a specific compound then filter out in LOOC process.
#'
#' @inheritParams de_genes
#'
#' @return
#' A data frame of gene identification result.
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' tgx_degs(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
tgx_degs <- function(...,
                      ge_matrix,
                      metadata,
                      p_cutoff = 0.05,
                      p_adjust = "none",
                      gr_diff = FALSE,
                      log10p = FALSE,
                      multicore = FALSE,
                      store = FALSE,
                      output_dir = rlang::missing_arg(),
                      error_call = rlang::caller_env()) {
  comps_group <- test_group(...)
  test_data(ge_matrix, metadata)
  output_dir <- destination(output_dir)
  compounds <- as.vector(unlist(comps_group))
  check_data <- update_data(
    compounds,
    ge_matrix = ge_matrix,
    metadata = metadata,
    multicore = multicore,
    store = store,
    output_dir = output_dir,
    error_call = error_call
  )
  ge_matrix <- check_data$expression
  metadata <- check_data$metadata
  full_model <- de_genes(
    comps_group,
    ge_matrix = ge_matrix,
    metadata = metadata,
    p_cutoff = p_cutoff,
    p_adjust = p_adjust,
    gr_diff = gr_diff,
    multicore = multicore,
    store = FALSE,
    output_dir = output_dir,
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
    hlm_tab <- full_model[match(all_genes$PROBEID, full_model$probe_id), ]
    hlm_df <- hlm_tab %>%
      dplyr::mutate(entrez_id = all_genes$ENTREZID,
                    gene_symbol = all_genes$SYMBOL,
                    gene_name = all_genes$GENENAME,
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
      multicore = FALSE,
      store = FALSE,
      output_dir = output_dir,
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
    ), call = error_call)
    return(structure(sig_tab, class = c("ToxAssay", "tbl_df", "tbl", "data.frame")))
}
