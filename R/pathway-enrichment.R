#' Get Gene Enrichment Results
#'
#' @description
#' The `get_enrichment` function enriches a set of genes based on specified categories such as KEGG pathways,
#' Gene Ontology (GO) biological processes, molecular functions, cellular components, Reactome, and WikiPathways.
#' It maps the genes to the STRING database, retrieves enrichment results, and processes the results to return a data frame of enriched pathways.
#'
#' @param gene_df A data frame containing the gene information with a column named `gene_symbol`.
#' @param category A character vector specifying the categories for enrichment. Options are "KEGG" (the default), "Process", "Component", "Function", "RCTM", and "WikiPathways".
#' @param path_n An integer specifying the number of top pathways to return (default is 10).
#' @param organism A character string specifying the organism. Options are "rat" or "human".
#' @param score_threshold A numeric value specifying the score threshold for STRING database interactions (default is 200).
#' @param version A character string specifying the version of the STRING database (default is "12").
#'
#' @return A tibble containing enriched pathways with columns:
#' \item{ID}{The term ID of the pathway.}
#' \item{pathway}{The name of the pathway.}
#' \item{n_genes}{The number of genes in the pathway.}
#' \item{p_value}{The p-value of the enrichment.}
#' \item{fdr}{The false discovery rate of the enrichment.}
#' \item{genes}{The genes involved in the pathway.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Sample gene data frame
#' gene_df <- data.frame(gene_symbol = c("Gene1", "Gene2", "Gene3"))
#'
#' # Get enrichment results using KEGG pathways for humans
#' enriched_genes <- get_enrichment(gene_df, category = "KEGG", path_n = 5, organism = "human")
#'
#' # View the enriched pathways
#' print(enriched_genes)
#' }
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
    cli_alert_warning(c("You have selected {path_n} pathway{?s}, but only enriched {total_path} pathway{?s}."))
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



#' Generate Enrichment Network Data
#'
#' @description
#' The `enrichment_netdata` function generates network like data for gene enrichment analysis based on specified categories such as KEGG pathways.
#' It maps the genes to pathways, calculates weights, and prepares the data for network visualization.
#'
#' @param gene_df A data frame containing the gene information with a column named `gene_symbol`.
#' @param category A character string specifying the category for enrichment (default is "KEGG").
#' @param organism A character string specifying the organism (default is "rat").
#' @param path_n An integer specifying the number of top pathways to return (default is 10).
#' @param gene_wtcol A character string specifying the column in `gene_df` to use for gene weights (default is "p_value").
#' @param score_threshold A numeric value specifying the score threshold for STRING database interactions (default is 200).
#' @param version A character string specifying the version of the STRING database (default is "12").
#'
#' @return A list containing two elements:
#' \item{vertices}{A tibble containing the vertices of the network with columns: `name`, `count`, `p_value`, `nodes`, and `size`.}
#' \item{edges}{A tibble containing the edges of the network with columns: `from`, `to`, and `weight`.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Sample gene data frame
#' gene_df <- data.frame(gene_symbol = c("Gene1", "Gene2", "Gene3"))
#'
#' # Generate enrichment network data for KEGG pathways in rats
#' net_data <- enrichment_netdata(gene_df, category = "KEGG", organism = "rat")
#'
#' # View the vertices and edges of the network
#' print(net_data$vertices)
#' print(net_data$edges)
#' }
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
  edge_data <- tibble(from = unlist(rep_path), to = unlist(gene_list), weight = NA) %>%
    filter(to %in% gene_name)
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

  gene_data <- left_join(gene_count, gene_pval, by = "name")
  vertx_data <- bind_rows(path_data, gene_data)

  net_data <- list(vertices = vertx_data, edges = edge_data)
  N <- sum(net_data$vertices$nodes == "gene")
  for(i in 1:nrow(net_data$edges)) {
    n1 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$from]
    n2 <- net_data$vertices$count[net_data$vertices$name %in% net_data$edges[i, ]$to]
    net_data$edges$weight[i] <- N / (n1 * n2)
  }
  return(net_data)
}
