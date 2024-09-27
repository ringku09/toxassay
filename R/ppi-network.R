#' Setup STRING Database
#'
#' The `setup_stringdb()` function initializes the STRING database for a specified organism and score threshold.
#'
#' @param organism A character string specifying the sample organism used for the experiment. Options are "human" or "rat". Default is "rat".
#' @param score_threshold A numeric value specifying the threshold for the combined scores of the interactions. Default is 200.
#' @param version A character string specifying the version of the STRING database. Default is "11.5".
#'
#' @return An object of class `STRINGdb` representing the STRING database for the specified organism.
#' @export
#'
#' @examples
#' \dontrun{
#' # Setup STRING database for rat with a score threshold of 200
#' string_db <- setup_stringdb(organism = "rat", score_threshold = 200)
#' }
setup_stringdb <- function(organism = c("human","rat"),
                           score_threshold = 200,
                           version = "12") {
  test_input(organism, auto_input = TRUE)
  if (identical(organism, "human")) {
    taxa = 9606
  } else if (identical(organism, "rat")) {
    taxa = 10116
  }
  string_db <- STRINGdb::STRINGdb$new(
    version = version,
    species = taxa,
    network_type = "full",
    score_threshold = score_threshold,
    input_directory = ""
    )
  return(string_db)
}

#' Get network data
#'
#' The function `get_netdata()` is used to extract network data from gene table.
#'
#' @param gene_df A data frame of gene information or a `tgxtool` class of object.
#' @param gene_col The name of the column for gene name in `gene_df`.
#' @param p.value_col The name of the column for p-values in `gene_df`.
#' @param cluster_method Cluster algorithm to find community. You can choose between "fastgreedy",
#' "walktrap", "spinglass" and "edge.betweenness" (the default is `NULL`).
#' @inheritParams setup_stringdb
#'
#' @return A list contain two data frames
#'  1) `vertices`: A data frame of information about nodes.
#'  2) `edges`: A data frame of information about edges.
#' @export
#'
#' @examples
#' \dontrun{
#' net_data <- get_ppinet(
#'   gene_data,
#'   gene_col = "gene_name",
#'   p.value_col = "group",
#'   organism = "rat",
#'   score_threshold = 400,
#'   cluster_method = "edge.betweenness"
#' )
#' }
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
  eig <- igraph::evcent(network)$vector     # Eigenvector centrality
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

#' Extract hub network data
#'
#' The function `get_hubdata()` is used to extract hub network data from full network data.
#'
#' @param net_data A network data of class `tgxtool`.
#' @param condition A character string specify the condition of hub network.
#' (default is `degree >= 20`). There are four metrics `degree`, `betweenness`, `closenes`,
#'  and `eigenes` can be used in the condition with the five comparison operators `<`, `>`, `<=`,
#'  `>=`, and `==`. At the end of the condition you must provide appropriate numeric value of metric
#'  used in the condition.
#'
#' @return A network data
#' @export
#'
#' @examples
#' \dontrun{
#' net_data <- get_netdata(
#'   gene_data,
#'   gene_col = "gene_name",
#'   p.value_col = "group",
#'   organism = "rat",
#'   score_threshold = 400,
#'   cluster_method = "edge.betweenness"
#' )
#' get_hubdata(net_data = net_dt, condition = "degree >= 20")
#' }
get_hubdata <- function(net_data, condition = "degree >= 10", error_call = caller_env()) {
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
