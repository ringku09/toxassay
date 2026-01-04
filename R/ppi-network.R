#' Initialize a STRING database object
#'
#' @description
#' `setup_stringdb()` creates and configures a `STRINGdb` object for retrieving
#' protein–protein interaction networks from the STRING database for a specified
#' organism. The function allows control over the STRING version, interaction
#' score threshold, and local storage location for downloaded data.
#'
#' @param organism Character string specifying the organism of interest.
#'   Must be one of `"human"` or `"rat"`. Default is `"rat"`.
#' @param score_threshold Numeric value specifying the minimum combined interaction
#'   score required to include an edge in the network. Default is `200`.
#' @param version Character string specifying the STRING database version to use.
#'   Default is `"12"`.
#' @param file_path Character string specifying the directory in which STRING data
#'   files will be stored. If the required files already exist at this location,
#'   downloading is skipped. If empty or `NULL`, a temporary directory is used.
#'
#' @details
#' STRING downloads can be large and may fail on slow connections when the global R
#' download timeout is low. To reduce download failures, the function temporarily sets
#' `options(timeout = 300)` when the current timeout is smaller, and restores the
#' original timeout value on exit.
#'
#' @return An object of class `STRINGdb` configured for the specified organism and
#'   interaction score threshold.
#'
#' @examples
#' \dontrun{
#' # Initialize the STRING database for rat with a score threshold of 200
#' string_db <- setup_stringdb(
#'   organism = "rat",
#'   score_threshold = 200,
#'   version = "12")
#' }
#'
#' @seealso
#' [STRINGdb::STRINGdb] for details on the STRING database interface.
#'
#' @export
setup_stringdb <- function(organism = c("human", "rat"),
                           score_threshold = 200,
                           version = "12",
                           file_path = "") {
  test_input(organism, auto_input = TRUE)
  current_timeout <- getOption("timeout")
  on.exit(options(timeout = current_timeout), add = TRUE)
  if (is.null(current_timeout) || current_timeout < 300) {
    options(timeout = 300)
    cli::cli_alert_info(
      "Temporarily increased download `timeout` to {style_bold(col_red(300))} seconds for STRING data retrieval."
    )
  }
  taxa <- if (identical(organism, "human")) {
    9606
  } else if (identical(organism, "rat")) {
    10116
  }
  string_db <- STRINGdb::STRINGdb$new(
    version = version,
    species = taxa,
    network_type = "full",
    score_threshold = score_threshold,
    input_directory = file_path
  )
  return(string_db)
}


#' Build a protein–protein interaction network from STRING
#'
#' @description
#' `get_ppinet()` maps input gene symbols to STRING identifiers and retrieves the corresponding
#' protein–protein interaction (PPI) subnetwork from the STRING database. The function returns
#' edge and vertex tables suitable for downstream network analysis and visualization. It also
#' computes common node centrality measures and assigns community labels using a selected
#' clustering algorithm.
#'
#' @param gene_df A data frame containing gene information. Must include a `gene_symbol` column.
#' @param organism Character string specifying the organism for STRING mapping. One of
#'   `"human"` or `"rat"`.
#' @param cluster_method Character string specifying the community detection algorithm used to
#'   cluster the network. One of `"edge.betweenness"`, `"fastgreedy"`, `"walktrap"`, or
#'   `"spinglass"`. Default is `"edge.betweenness"` (selected via `test_input(auto_input = TRUE)`).
#' @param n_percent Numeric value (percentage) used to define the minimum community size. Communities
#'   with fewer than `ceiling(n_genes * n_percent / 100)` members are merged into `"SN0"`.
#' @param score_threshold Integer specifying the minimum STRING interaction score passed to
#'   `setup_stringdb()`.
#' @param version Character string specifying the STRING database version passed to
#'   `setup_stringdb()` (e.g., `"12"`).
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Initializes a STRINGdb object via `setup_stringdb()` using `organism`, `score_threshold`,
#'     and `version`.
#'   \item Maps `gene_symbol` to STRING identifiers and removes duplicated STRING IDs.
#'   \item Retrieves the STRING subnetwork induced by the mapped genes.
#'   \item Returns network tables (`edges`, `vertices`) and augments vertices with degree,
#'     betweenness, closeness, and eigenvector centralities.
#'   \item Detects communities using `cluster_method`; small communities are merged into `"SN0"`.
#' }
#'
#' @return A list with two data frames:
#' \itemize{
#'   \item `vertices`: node-level information including `gene_symbol`, `STRING_id`, centrality measures,
#'     and community annotations (`community`, `gene_class`)
#'   \item `edges`: edge-level information with mapped `from` and `to` gene symbols
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
#' net_data <- get_ppinet(
#'   gene_data,
#'   organism = "human",
#'   score_threshold = 400,
#'   cluster_method = "edge.betweenness"
#' )
#' }
#'
#' @seealso
#' [setup_stringdb()] for creating the STRINGdb interface used internally.
#'
#' @export
get_ppinet <- function(gene_df,
                        organism = c("human", "rat"),
                        cluster_method = c("edge.betweenness", "fastgreedy", "walktrap", "spinglass"),
                        n_percent = 5,
                        score_threshold = 200,
                        version = "12") {
  test_column("gene_symbol", gene_df)
  string_db <- setup_stringdb(organism = organism, score_threshold = score_threshold, version = version)
  gene_map <- string_db$map(as.data.frame(gene_df), "gene_symbol", removeUnmappedRows = TRUE)
  gene_map <- gene_map[!duplicated(gene_map$STRING_id),]
  net <- string_db$get_subnetwork(gene_map$STRING_id)
  net_df <- suppressWarnings(igraph::as_data_frame(net, what = "both"))
  net_df$edges <- tibble::as_tibble(net_df$edges) %>%
    dplyr::mutate(
      from = gene_map$gene_symbol[match(net_df$edges$from, gene_map$STRING_id)],
      to = gene_map$gene_symbol[match(net_df$edges$to, gene_map$STRING_id)]
      ) %>%
    dplyr::distinct()
  net_df$vertices <- gene_map  #[match(gene_map$STRING_id, net_df$vertices$name),]
  network <- igraph::graph_from_data_frame(d = net_df$edges, vertices = net_df$vertices, directed = TRUE)
  deg <- igraph::degree(network)            # Degree centrality
  clo <- igraph::closeness(network)         # Closeness centrality
  bet <- igraph::betweenness(network)       # Betweenness centrality
  eig <- igraph::eigen_centrality(network)$vector     # Eigenvector centrality
  net_df$vertices <- net_df$vertices %>%
    dplyr::mutate(
      degree = deg,
      betweenness = bet,
      closenes = clo,
      eigenes = eig) %>%
    tibble::as_tibble()
  test_input(cluster_method, auto_input = TRUE)
  cl_list <- string_db$get_clusters(net_df$vertices$STRING_id, algorithm = cluster_method)
  min_mem <- ceiling(sum(length(gene_map$STRING_id)) * n_percent/100)
  act_gr <- lapply(cl_list, function(x) {length(x) >= min_mem})
  others <- list(unlist(cl_list[!unlist(act_gr)]))
  names(others) <- "SN0"
  act_comunity <- cl_list[unlist(act_gr)]
  names(act_comunity) <- paste0("SN", 1:length(act_comunity))
  comunity <- append(act_comunity, others)
  if (length(comunity) > 10) {
    comunity9 <- comunity[1 : 9]
    others <- list(as.vector(unlist(comunity[10 : length(comunity)])))
    names(others) <- "SN0"
    comunity <- append(comunity9, others)
  }
  comunity <- lapply(comunity, function(x) {gene_map$gene_symbol[match(x, gene_map$STRING_id)]})
  cl_name <- 1 : length(comunity)
  community <- unlist(lapply(lapply(net_df$vertices$gene_symbol, function(x) {grep(x, comunity)}), "[",1))
  gene_class <- names(comunity)[unlist(lapply(lapply(net_df$vertices$gene_symbol, function(x) {grep(x, comunity)}), "[",1))]
  net_df$vertices <- net_df$vertices %>%
    dplyr::mutate(community = community,
                  gene_class = factor(gene_class, levels = names(comunity)),
                  STRING_id = sub(".*\\.", "", STRING_id))
  return(net_df)
}

#' Extract a hub subnetwork from STRING-derived network data
#'
#' @description
#' `get_hubdata()` filters a full network object (as returned by `get_ppinet()`) to retain a
#' hub subnetwork defined by a user-specified condition on node centrality metrics. Nodes that
#' satisfy the condition are kept, and the edge table is restricted to interactions where both
#' endpoints are among the retained hub nodes.
#'
#' @param net_data A network object containing `vertices` and `edges` components, as returned by
#'   `get_ppinet()`. `net_data$vertices` must contain a `gene_symbol` column and centrality metrics.
#' @param condition Character string specifying a filtering rule applied to `net_data$vertices`.
#'   The rule must be of the form `<metric><operator><value>` (spaces are allowed but ignored),
#'   for example `"degree >= 10"`. Valid metrics are `degree`, `betweenness`, `closenes`,
#'   and `eigenes`. Valid operators are `<`, `>`, `<=`, `>=`, and `==`. The condition must end
#'   with a numeric threshold. Default is `"degree >= 10"`.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @details
#' The function parses `condition` to extract:
#' \itemize{
#'   \item a metric name (one of `degree`, `betweenness`, `closenes`, `eigenes`)
#'   \item a comparison operator (`<`, `>`, `<=`, `>=`, `==`)
#'   \item a numeric cutoff value
#' }
#' The condition is then evaluated on `net_data$vertices`, and only rows meeting the condition
#' are retained. The edge list is subsequently filtered to keep only edges connecting two retained
#' hub nodes.
#'
#' @return A network object with the same structure as `net_data` (components `vertices` and `edges`),
#' containing only the hub nodes and their induced subnetwork.
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
#' net_data <- get_ppinet(
#'   gene_data,
#'   organism = "human",
#'   score_threshold = 400,
#'   cluster_method = "edge.betweenness"
#' )
#' hub_net <- get_hubdata(net_data = net_data, condition = "degree >= 10")
#' }
#'
#' @seealso
#' [get_ppinet()] for generating the full network with centrality metrics used by `condition`.
#'
#' @export
get_hubdata <- function(net_data,
                        condition = "degree >= 10",
                        error_call = caller_env()) {
  test_element(names(net_data), c("vertices", "edges"))
  test_column("gene_symbol", net_data$vertice)
  condition <- gsub(" ", "", condition)
  attrb <- gsub("[^(A-Za-z)]", "", condition)
  value <- as.numeric(gsub("[^0-9.-]", "", condition))
  cond <- gsub("[a-zA-Z0-9.]", "", condition)

  if (!any(sapply(list("degree", "betweenness", "closenes", "eigenes"), FUN = identical, attrb))) {
    cli::cli_abort(c("Metric used in condition must be a valid metric.",
                "x" = "Input {style_bold(col_red(backtick(attrb)))} is not a valid metric.",
                "i" = "Please use either `degree`, `betweenness`, `closenes`, or `eigenes` instread.")
              , call = error_call)
  }
  if (!any(sapply(list(">", "<", ">=", "<=", "=="), FUN = identical, cond))) {
    cli::cli_abort(c("Comparison operator used in condition must be a valid comparison operator.",
                "x" = "Input {style_bold(col_red(backtick(cond)))} is not a valid metric.",
                "i" = "Please use either `<`, `>`, `<=`, `>=`, or `==` instread.")
              , call = error_call)
  }
  if (!is.numeric(value) | is.na(value)) {
    cli::cli_abort(c("The condition must contain numeric value after symbol.",
                "x" = "The numeric value in condition are missing.",
                "i" = "Please provide appropriate numeric value in the condition.")
              , call = error_call)
  }
  valid_cond <- net_data$vertices %>%
    dplyr::select(!!rlang::parse_expr(attrb)) %>%
    unlist()
  if (all(value > valid_cond)) {
    max_val <- max(valid_cond)
    cli::cli_abort(c("The numeric value must be valid with condition.",
                "x" = "The numeric value {value} in condition exceed the real values.",
                "i" = "Please provide value smaller than {max_val} in the condition.")
              , call = error_call)
  }
  net_data$vertices <- net_data$vertices %>%
    dplyr::filter(!!rlang::parse_expr(condition)) %>%
    droplevels()
  hub_from <- (1 : nrow(net_data$edges))[net_data$edges$from %in% net_data$vertices$gene_symbol]
  hub_to <- (1 : nrow(net_data$edges))[net_data$edges$to %in% net_data$vertices$gene_symbol]
  hub_node <- intersect(hub_from, hub_to)
  net_data$edges <- net_data$edges[hub_node, ]
  return(net_data)
}
