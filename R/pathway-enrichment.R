#' Retrieve functional enrichment results from STRING
#'
#' @description
#' `get_enrichment()` performs functional enrichment analysis for a set of genes using the
#' STRING database. Gene symbols are mapped to STRING identifiers, enrichment is computed for a
#' selected annotation category (e.g., KEGG, GO, Reactome, WikiPathways), and the top enriched
#' terms are returned as a tidy table.
#'
#' @param gene_df A data frame containing gene information. Must include a `gene_symbol` column.
#' @param category Character string specifying the enrichment category to query in STRING.
#'   One of `"KEGG"`, `"Process"` (GO Biological Process), `"Component"` (GO Cellular Component),
#'   `"Function"` (GO Molecular Function), `"RCTM"` (Reactome), or `"WikiPathways"`.
#'   Default is `"KEGG"` (selected via `test_input(auto_input = TRUE)`).
#' @param path_n Integer specifying the number of top enriched terms to return. If `NULL`, all
#'   enriched terms are returned. Default is `10`.
#' @param organism Character string specifying the organism used for STRING mapping. One of
#'   `"rat"` or `"human"`.
#' @param score_threshold Integer specifying the minimum STRING interaction score passed to
#'   `setup_stringdb()`. Default is `200`.
#' @param version Character string specifying the STRING database version passed to
#'   `setup_stringdb()` (e.g., `"12"`). Default is `"12"`.
#'
#' @details
#' The function maps `gene_symbol` values to STRING IDs using `setup_stringdb()`. Enrichment
#' results are retrieved using `STRINGdb::get_enrichment()` for the selected `category`. If `path_n`
#' is larger than the number of available enriched terms, all terms are returned and a warning is
#' shown. The returned table includes the term identifiers, descriptions, gene counts, and
#' enrichment statistics. A comma-separated `genes` column is constructed from STRING's
#' `preferredNames` field.
#'
#' @return A tibble with one row per enriched term and the columns:
#' \itemize{
#'   \item `ID`: term identifier (e.g., KEGG/Reactome/WikiPathways/GO term ID)
#'   \item `pathway`: term description
#'   \item `n_genes`: number of mapped genes in the term
#'   \item `p_value`: enrichment p-value
#'   \item `fdr`: false discovery rate
#'   \item `genes`: comma-separated gene symbols mapped to the term
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Get probes for "Insulin signaling pathway" (Rat KEGG ID: 04910)
#' insulin_probes <- AnnotationDbi::select(rat2302.db::rat2302.db,
#'                                         keys = "04910",
#'                                         keytype = "PATH",
#'                                         columns = c("PROBEID", "SYMBOL"))
#' # Simulate gene expression data
#' sim_data <- simulate_data(n_gene = nrow(insulin_probes), n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' gene_data <- tox_degs(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
#' gene_data$probe_id <- insulin_probes$PROBEID
#' gene_data$gene_symbol <- insulin_probes$SYMBOL
#' enriched_genes <- get_enrichment(gene_data)
#' }
#'
#' @seealso
#' [setup_stringdb()] for initializing the STRING interface used internally.
#'
#' @export
get_enrichment <- function(gene_df,
                           category = c("KEGG", "Process", "Component", "Function", "RCTM", "WikiPathways"),
                           path_n = 10,
                           organism = c("rat", "human"),
                           score_threshold = 200,
                           version = "12") {
  test_column("gene_symbol", gene_df)
  test_input(category, auto_input = TRUE)
  test_input(organism, auto_input = TRUE)
  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold, version = version)
  gene_map <- string_db$map(as.data.frame(gene_df), "gene_symbol", removeUnmappedRows = TRUE)
  enrichment <- string_db$get_enrichment(gene_map, category = category)
  total_path <- nrow(enrichment)
  if (is.null(path_n)) {
    enrich_res <- enrichment
  } else if (path_n <= total_path) {
    enrich_res <- enrichment[1 : path_n, ]
  } else if (path_n > total_path) {
    cli::cli_alert_warning(c("You have selected {path_n} pathway{?s}, but only enriched {total_path} pathway{?s}."))
    enrich_res <- enrichment[1 : total_path, ]
  }
  path_tab <- enrich_res %>%
    dplyr::select(c(category, term, description, number_of_genes, p_value, fdr)) %>%
    dplyr::rename(ID = term, pathway = description, n_genes = number_of_genes) %>%
    dplyr::as_tibble()
  if (total_path > 0) {
    gene_list <- strsplit(enrich_res$preferredNames, ",")
    genes_path <- unlist(lapply(gene_list, function(x) {paste(x, collapse = ", ")}))
    path_tab <- path_tab %>%
      dplyr::mutate(genes = genes_path)
  } else {
    path_tab <- path_tab %>%
      dplyr::mutate(genes = enrich_res$preferredNames)
  }
  return(path_tab)
}

#' Create enrichment network data for pathway–gene visualization
#'
#' @description
#' `enrichment_netdata()` converts STRING enrichment results into a bipartite network-like
#' representation linking enriched terms (e.g., pathways or GO categories) to member genes.
#' The output is designed for network visualization, with vertex attributes describing
#' pathways and genes, and edge weights reflecting term/gene connectivity.
#'
#' @param gene_df A data frame containing gene information. Must include a `gene_symbol` column.
#' @param category Character string specifying the enrichment category queried in STRING.
#'   One of `"KEGG"`, `"Process"`, `"Component"`, `"Function"`, `"RCTM"`, or `"WikiPathways"`.
#' @param organism Character string specifying the organism for STRING mapping. One of `"rat"` or
#'   `"human"`.
#' @param path_n Integer specifying the number of top enriched terms to include. If `NULL`, all
#'   enriched terms are used. Default is `10`.
#' @param gene_wtcol Column name in `gene_df` used as a gene-level weight when assigning gene node
#'   sizes (e.g., `"p_value"`). If `NULL`, a constant size is used for all genes. Default is
#'   `"p_value"`.
#' @param score_threshold Integer specifying the minimum STRING interaction score passed to
#'   `setup_stringdb()`. Default is `200`.
#' @param version Character string specifying the STRING database version passed to
#'   `setup_stringdb()` (e.g., `"12"`). Default is `"12"`.
#'
#' @details
#' The function first calls `get_enrichment()` to obtain enriched terms and their member genes.
#' It then constructs:
#' \itemize{
#'   \item a pathway vertex table with term size scaled from enrichment p-values
#'   \item a gene vertex table restricted to genes present in the enrichment results, with size
#'     optionally scaled from `gene_wtcol`
#'   \item an edge table connecting each term to its genes
#' }
#' Edge weights are computed as `N / (n1 * n2)`, where `N` is the number of gene nodes, `n1` is the
#' term gene count, and `n2` is the number of terms connected to the gene (gene degree in the
#' bipartite graph).
#'
#' @return A list with two tibbles:
#' \itemize{
#'   \item `vertices`: vertex attributes with columns including `name`, `count`, `p_value` (for genes),
#'     `nodes` (`"pathway"` or `"gene"`), and `size`
#'   \item `edges`: edge list with columns `from` (term name), `to` (gene symbol), and `weight`
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Get probes for "Insulin signaling pathway" (Rat KEGG ID: 04910)
#' insulin_probes <- AnnotationDbi::select(rat2302.db::rat2302.db,
#'                                         keys = "04910",
#'                                         keytype = "PATH",
#'                                         columns = c("PROBEID", "SYMBOL"))
#' # Simulate gene expression data
#' sim_data <- simulate_data(n_gene = nrow(insulin_probes), n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' gene_data <- tox_degs(gr, ge_matrix = sim_data$expression, metadata = sim_data$metadata)
#' gene_data$probe_id <- insulin_probes$PROBEID
#' gene_data$gene_symbol <- insulin_probes$SYMBOL
#' enriched_net <- enrichment_netdata(gene_data)
#' }
#'
#' @seealso
#' [get_enrichment()] for retrieving STRING enrichment results used to build the network.
#'
#' @export
enrichment_netdata <- function(gene_df,
                               category = c("KEGG", "Process", "Component", "Function", "RCTM", "WikiPathways"),
                               organism = c("rat", "human"),
                               path_n = 10,
                               gene_wtcol = "p_value",
                               score_threshold = 200,
                               version = "12") {
  enrc_tab <- get_enrichment(gene_df = gene_df,
                             category = category,
                             path_n = path_n,
                             organism = organism,
                             score_threshold = score_threshold,
                             version = version)
  gene_list <- strsplit(enrc_tab$genes, split = "\\, ")
  gene_n <- as.list(enrc_tab$n_genes)
  rep_path <- mapply(function(x, y) {rep(y, x)},x = gene_n, y = enrc_tab$pathway)
  gene_name <- gene_df %>%
    dplyr::select(gene_symbol) %>%
    dplyr::pull()
  edge_data <- tibble::tibble(from = unlist(rep_path), to = unlist(gene_list), weight = NA) %>%
    dplyr::filter(to %in% gene_name)
  enrc_genes <- unique(edge_data$to)
  if(any(grepl(".*e-", enrc_tab$p_value))){
    weight <- format(enrc_tab$p_value, scientific = TRUE)
    path_size <- ceiling(scales::rescale(as.numeric(gsub(".*e-","", weight)), to = c(10, 15)))
  } else {
    path_size <- ceiling(scales::rescale(enrc_tab$p_value, to = c(10, 15)))
  }
  path_data <- enrc_tab %>%
    dplyr::select(c(pathway, n_genes, p_value)) %>%
    dplyr::mutate(nodes = "pathway", size = path_size) %>%
    dplyr::rename(name = pathway, count = n_genes)
  gene_df2 <- gene_df %>%
    dplyr::filter(gene_symbol %in% enrc_genes)
  if (is.null(gene_wtcol)) {
    gene_size <- 5
  } else {
    weight <- gene_df2 %>%
      dplyr::select({{gene_wtcol}}) %>%
      dplyr::pull()
    if(any(grepl(".*e-", weight))){
      weight <- format(weight, scientific = TRUE)
      gene_size <- ceiling(scales::rescale(as.numeric(gsub(".*e-","", weight)), to = c(3, 7)))
    } else {
      gene_size <- ceiling(scales::rescale(weight, to = c(3, 7)))
    }
  }
  gene_pval <- gene_df2 %>%
    dplyr::select(c(gene_symbol,{{gene_wtcol}})) %>%
    dplyr::mutate(nodes = "gene", size = gene_size) %>%
    dplyr::rename(name = gene_symbol, p_value = {{gene_wtcol}})
  gene_count <- edge_data %>%
    dplyr::rename(name = to) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(count = n())

  gene_data <- dplyr::left_join(gene_count, gene_pval, by = "name")
  vertx_data <- dplyr::bind_rows(path_data, gene_data)

  net_data <- list(vertices = vertx_data, edges = edge_data)
  N <- sum(net_data$vertices$nodes == "gene")
  for(i in 1:nrow(net_data$edges)) {
    n1 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$from]
    n2 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$to]
    net_data$edges$weight[i] <- N / (n1 * n2)
  }
  return(net_data)
}
