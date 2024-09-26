#' Identify Core Differentially Expressed Genes (DEGs)
#'
#' This function identifies core differentially expressed genes (DEGs) using
#' multiple clustering methods in PPI network and pathway enrichment analysis.
#'
#' @param gene_df A data frame containing the gene information with a column named `gene_symbol`.
#' @param path_n Number of pathways to consider for enrichment analysis. Default is `NULL`.
#' @param enrich_cutoff Cutoff value for enrichment analysis. Default is `0.5`.
#' @param category Vector of categories for pathway enrichment analysis. Options are
#'  `kegg`, `go`, `reactome`, and `wikipathways`. Default is `kegg` if not supplied.
#' @param organism Organism of interest. Options are `"rat"` and `"human"`. Default is `c("rat", "human")`.
#' @param score_threshold Score threshold for filtering genes. Default is `200`.
#' @param version Version of the database to use for pathway enrichment analysis. Default is `"12"`.
#'
#' @return A vector of pivot genes that are core DEGs across different clustering methods and enrichment analyses.
#'
#' @details This function performs the following steps:
#'   \itemize{
#'     \item Validates the input categories and organisms.
#'     \item Retrieves network data for the specified gene expression data using three clustering methods: "fastgreedy", "walktrap", and "edge.betweenness".
#'     \item Performs pathway enrichment analysis for the specified categories and organism.
#'     \item Identifies clusters in the network data and performs enrichment analysis on each cluster.
#'     \item Calculates the proportion of enriched pathways shared with the overall enrichment analysis.
#'     \item Filters clusters based on the enrichment cutoff and extracts unique genes from the filtered clusters.
#'     \item Identifies and returns pivot genes that are common across the filtered clusters from different clustering methods.
#'   }
#'
#' @examples
#' \dontrun{
#'   core_genes <- core_degs(gene_df = res4, path_n = 10, enrich_cutoff = 0.5, category = "kegg", organism = "human", score_threshold = 200, version = "12")
#' }
#'
#' @export
core_degs <- function(gene_df = res4,
                        path_n = NULL,
                        enrich_cutoff = 0.5,
                        category = c("kegg", "go", "reactome", "wikipathways"),
                        organism = c("rat", "human"),
                        score_threshold = 200,
                        version = "12") {
  test_input(category, auto_input = TRUE)
  test_input(organism, auto_input = TRUE)
  input_db <- c(kegg = "KEGG",
                go = "Process",
                reactome = "RCTM",
                wikipathways = "WikiPathways")
  category <- input_db[category]
  net_fgreedy <- get_netdata(gene_df = gene_df,
                             organism = organism,
                             cluster_method = "fastgreedy",
                             score_threshold = score_threshold,
                             version = version)
  net_walktrap <- get_netdata(gene_df = gene_df,
                              organism = organism,
                              cluster_method = "walktrap",
                              score_threshold = score_threshold,
                              version = version)
  net_edgebet <- get_netdata(gene_df = gene_df,
                             organism = organism,
                             cluster_method = "edge.betweenness",
                             score_threshold = score_threshold,
                             version = version)

  enrc_path <- get_enrichment(gene_df = gene_df,
                              category = category,
                              organism = organism,
                              path_n = path_n,
                              score_threshold = score_threshold,
                              version = version)
  freddy_class <- unique(net_fgreedy$vertices$community)
  pect_fgreedy <- vector(length = length(freddy_class))
  enrc_fgreedy <- vector(mode = "list", length = length(freddy_class))
  for(i in 1:length(freddy_class)) {
    cluster_fgreedy <- net_fgreedy$vertices %>%
      dplyr::filter(community == i)
    gcluster_fgreedy <- gene_df %>%
      dplyr::filter(gene_symbol %in% cluster_fgreedy$gene_symbol)
    enrc_fgreedy[[i]] <- get_enrichment(gene_df = gcluster_fgreedy,
                                        category = category,
                                        organism = organism,
                                        path_n = path_n,
                                        score_threshold = score_threshold,
                                        version = version)
    pect_fgreedy[i] <- round(sum(enrc_fgreedy[[i]]$ID %in% enrc_path$ID)/length(enrc_path$ID),2)
  }
  walk_class <- unique(net_walktrap$vertices$community)
  pect_walk <- vector(length = length(walk_class))
  enrc_walk <- vector(mode = "list", length = length(walk_class))
  for(i in 1:length(walk_class)) {
    cluster_walk <- net_walktrap$vertices %>%
      dplyr::filter(community == i)
    gcluster_walk <- gene_df %>%
      dplyr::filter(gene_symbol %in% cluster_walk$gene_symbol)
    enrc_walk[[i]] <- get_enrichment(gene_df = gcluster_walk,
                                   category = category,
                                   organism = organism,
                                   path_n = path_n,
                                   score_threshold = score_threshold,
                                   version = version)
    pect_walk[i] <- round(sum(enrc_walk[[i]]$ID %in% enrc_path$ID)/length(enrc_path$ID),2)
  }

  edge_class <- unique(net_edgebet$vertices$community)
  pect_edge <- vector(length = length(edge_class))
  enrc_edge <- vector(mode = "list", length = length(edge_class))
  for(i in 1:length(edge_class)) {
    cluster_edge <- net_edgebet$vertices %>%
      dplyr::filter(community == i)
    gcluster_edge <- gene_df %>%
      dplyr::filter(gene_symbol %in% cluster_edge$gene_symbol)
    enrc_edge[[i]] <- get_enrichment(gene_df = gcluster_edge,
                                   category = category,
                                   organism = organism,
                                   path_n = path_n,
                                   score_threshold = score_threshold,
                                   version = version)
    pect_edge[i] <- round(sum(enrc_edge[[i]]$ID %in% enrc_path$ID)/length(enrc_path$ID),2)
  }
  subnet_fgreedy <- enrc_fgreedy[pect_fgreedy >= enrich_cutoff]
  if (length(subnet_fgreedy) > 0) {
    genes_fgreedy <- unlist(strsplit(unlist(lapply(subnet_fgreedy, function(x) x$genes)), ", "))
    unique_fgreedy <- unique(genes_fgreedy)
  } else unique_fgreedy <- NULL

  subnet_walk <- enrc_walk[pect_walk >= enrich_cutoff]
  if (length(subnet_walk) > 0) {
  genes_walk <- unlist(strsplit(unlist(lapply(subnet_walk, function(x) x$genes)), ", "))
  unique_walk <- unique(genes_walk)
  } else unique_walk <- NULL

  subnet_edge <- enrc_edge[pect_edge >= enrich_cutoff]
  if (length(subnet_edge) > 0) {
  genes_edge <- unlist(strsplit(unlist(lapply(subnet_edge, function(x) x$genes)), ", "))
  unique_edge <- unique(genes_edge)
  } else unique_edge <- NULL

  gene_list <- list(unique_fgreedy, unique_walk, unique_edge)
  # gene_list <- Filter(Negate(is.null), gene_list)
  pivot_genes <- Reduce(dplyr::intersect, gene_list)
  core_genes <- pivot_genes[pivot_genes %in% gene_df$gene_symbol]
 return(core_genes)
}


#
# xx3 <- core_degs(res4,
#                 path_n = NULL,
#                 enrich_cutoff = 0.25,
#                 category = "go",
#                 organism = "rat",
#                 score_threshold = 200,
#                 version = "12")
