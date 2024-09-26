#' Get Transaction Data for Compounds, Genes, and Diseases
#'
#' @description
#' The `get_transaction()` function processes and retrieves transaction data for specified compounds,
#' genes, and diseases. It takes two data frames as input and filters the data based on the provided
#' compounds, gene IDs, and disease IDs. The function merges the gene and disease data and returns
#' a binary matrix indicating the presence of interactions between chemicals, genes, and diseases.
#'
#' @param compound_gene A data frame containing chemical and gene information.
#'  It should have columns `ChemicalName`, `GeneSymbol`, and `OrganismID`, similar to the CTD database.
#'  One can download this data using the function `get_ctd()`.
#' @param compound_disease A data frame containing chemical and disease information.
#'  It should have columns `ChemicalID` and `DiseaseID`, similar to the CTD database.
#'  One can download this data using the function `get_ctd()`.
#' @param compounds An optional character vector of compound names to filter.
#'  If `NULL` (default), all compounds from `compound_gene` are used.
#' @param gene_id An optional character vector of gene IDs to filter.
#'  If `NULL` (default), all gene IDs from `compound_gene` are used.
#' @param disease_id An optional character vector of disease IDs to filter.
#'  If `NULL` (default), all disease IDs from `compound_disease` are used.
#'
#' @return A list with three components:
#' \item{bin_data}{A binary matrix indicating the presence of interactions between chemicals, genes, and diseases.}
#' \item{genes}{A vector of gene symbols included in the binary matrix.}
#' \item{diseases}{A vector of disease IDs included in the binary matrix.}
#'
#' @examples
#' chemical_name <- paste0("Compound", 1:10)
#' gene_symbol <- paste0("GENE", 1:30)
#' disease_id <- paste0("Disease", 1:30)
#' compound_gene <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                             gene_symbol = sample(gene_symbol, 100, replace = TRUE),
#'                             organism_id = sample(c(9606, 10116), 100, replace = TRUE)) %>%
#'                             dplyr::distinct()
#' compound_disease <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                                disease_id = sample(disease_id, 100, replace = TRUE)) %>%
#'                                dplyr::distinct()
#' get_transaction(compound_gene, compound_disease)
#' @import arules
#' @importFrom reshape2 dcast
#' @export
get_transaction <- function(compound_gene,
                            compound_disease,
                            compounds = NULL,
                            genes = NULL,
                            diseases = NULL) {

  test_column(c("chemical_name", "gene_symbol", "organism_id"), compound_gene)
  test_column(c("chemical_name", "disease_id"), compound_disease)
  test_element(unique(na.omit(compound_gene$chemical_name)), unique(na.omit(compound_disease$chemical_name)))
  if (is.null(compounds)) {
    compounds <- compound_gene %>%
      dplyr::select(chemical_name) %>%
      unique() %>% dplyr::pull()
  }
  if (!is.null(genes)) {
    gene_idx <- toupper(genes) %in% unique(compound_gene$gene_symbol)
    if (all(!gene_idx)) {
      test_element(genes, unique(compound_gene$gene_symbol))
    } else {
      genes <- toupper(genes)[gene_idx]
    }
  }
  if (!is.null(diseases)) {
    disease_idx <- diseases %in% unique(compound_disease$disease_id)
    if (all(!disease_idx)) {
      test_element(diseases, unique(compound_disease$disease_id))
    } else {
      diseases <- diseases[disease_idx]
    }
  }
  if (is.null(diseases)) {
    diseases <- compound_disease %>% select(disease_id) %>% dplyr::pull()
  }
  if (is.null(genes)) {
    genes <- compound_gene %>% select(gene_symbol) %>% dplyr::pull()
  }
  gene_df <- compound_gene %>%
    dplyr::filter(chemical_name %in% compounds) %>%
    dplyr::filter(gene_symbol %in% toupper(genes)) %>%
    dplyr::filter(organism_id %in% c(9606, 10116)) %>%
    dplyr::select(c(chemical_name, gene_symbol)) %>%
    dplyr::distinct()
  dise_df <- compound_disease %>%
    dplyr::filter(chemical_name %in% compounds) %>%
    dplyr::filter(disease_id %in% diseases) %>%
    dplyr::select(c(chemical_name, disease_id)) %>%
    dplyr::distinct()
  gene_mat <- reshape2::dcast(gene_df, chemical_name ~ gene_symbol,
                              fun.aggregate = length, value.var = "gene_symbol")
  gene_mat[is.na(gene_mat)] <- 0
  dise_mat <- reshape2::dcast(dise_df, chemical_name ~ disease_id,
                              fun.aggregate = length, value.var = "disease_id")
  dise_mat[is.na(dise_mat)] <- 0
  genes <- colnames(gene_mat)[-1]
  diseases <- colnames(dise_mat)[-1]

  marg_df <- gene_mat %>%
    dplyr::inner_join(dise_mat, by = "chemical_name") %>%
    dplyr::distinct(chemical_name, .keep_all = TRUE)
  marg_df <- marg_df %>%
    dplyr::mutate(dplyr::across(-chemical_name, ~ ifelse(!is.na(.) & . > 0, 1, 0)))
  non_zero_cols <- colnames(marg_df)[colSums(marg_df[-1], na.rm = TRUE) > 1]
  non_zero_rows <- rowSums(marg_df[-1], na.rm = TRUE) > 1
  marg_df <- marg_df %>%
    dplyr::filter(non_zero_rows) %>%
    dplyr::select(chemical_name, all_of(non_zero_cols))
  marg_mat <- df2matrix(dplyr::select(marg_df,-chemical_name), dplyr::pull(marg_df, chemical_name))

  genes_idx <- which(genes %in% colnames(marg_mat))
  diseases_idx <- which(diseases %in% colnames(marg_mat))
  temp <- matrix(0, nrow = 4, ncol = ncol(marg_mat))
  temp[1, (length(genes_idx) + 1):ncol(marg_mat)] <- 1
  temp[2, 1:length(genes_idx)] <- 1
  temp[4, ] <- 1
  colnames(temp) <- colnames(marg_mat)
  rownames(temp) <- paste("temp", 1:4, sep = "")
  bin_data <- rbind(marg_mat, temp)
  return(list(bin_data, genes, diseases))
}


#' Get Adverse Outcome Pathways (AOPs) from Transaction Data
#'
#' @description
#' The `get_aops()` function processes transaction data to identify and retrieve
#' AOPs based on specified genes and diseases. The function uses the Apriori algorithm
#' to generate association rules, calculates confidence intervals,
#' and returns a data frame of the resulting rules with additional quality measures.
#'
#' @param transaction A binary matrix indicating the presence of interactions between chemicals, genes, and diseases.
#' @param genes A vector of gene symbols included in the binary matrix.
#' @param diseases A vector of disease IDs included in the binary matrix.
#' @param ci_metric A character string specifying the metric to use for confidence intervals. Default is `"lift"`.
#'
#' @return A tibble containing the resulting association rules with columns for the left-hand side (LHS) and right-hand side (RHS) of each rule,
#' along with various quality measures including support, confidence, lift, odds ratio, and confidence intervals.
#'
#' @examples
#' chemical_name <- paste0("Compound", 1:10)
#' gene_symbol <- paste0("GENE", 1:30)
#' disease_id <- paste0("Disease", 1:30)
#' compound_gene <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                             gene_symbol = sample(gene_symbol, 100, replace = TRUE),
#'                             organism_id = sample(c(9606, 10116), 100, replace = TRUE)) %>%
#'                             dplyr::distinct()
#' compound_disease <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                                disease_id = sample(disease_id, 100, replace = TRUE)) %>%
#'                                dplyr::distinct()
#' trans_data <- get_transaction(compound_gene, compound_disease)
#' get_aops(transaction = trans_data[[1]], genes = trans_data[[2]], diseases = trans_data[[3]])
#' @importFrom stats confint
#' @export
get_aops <- function(transaction,
                     genes,
                     diseases,
                     min_support = 0.3,
                     min_confidence = 0.5,
                     ci_metric = "lift") {
  if (!all(requireNamespace("arules", quietly = TRUE))) {
    cli::cli_abort(c("Packages `arules` required for AOP!",
                "i" = "Please install `arules`."),
              call = error_call)
  }
  transaction[transaction == 0] <- NA
  transaction <- as.data.frame(transaction) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), factor))
  bin_trans <- methods::as(transaction, "transactions")
  gene_items <- arules::intersect(colnames(bin_trans), paste(genes, 1, sep='='))
  dise_items <- arules::intersect(colnames(bin_trans), paste(diseases, 1, sep='='))
  rules <- suppressWarnings(arules::apriori(bin_trans,
                           parameter = list(support = min_support,
                                            confidence = min_confidence,
                                            minlen = 2,
                                            maxlen = 2)))
  rules_ap <- arules::subset(rules, subset = lhs %in% gene_items &
                               lift > 1 &
                               rhs %in% dise_items)
  ci <- stats::confint(rules_ap, ci_metric,  smoothCounts = 0.5, transactions = bin_trans)
  arules::quality(rules_ap) <- cbind(
    arules::quality(rules_ap),
    oddsRatio = arules::interestMeasure(rules_ap, "oddsRatio", bin_trans),
    CI = ci)
  df_ap <- methods::as(rules_ap, "data.frame")
  df_ap$rules <- as.character(df_ap$rules)
  temp <- sapply(strsplit(df_ap$rules, split="[{}]"), unlist)
  df_ap <- df_ap %>%
    dplyr::mutate(lhs = gsub("=1", "",temp[2,]), rhs = gsub("=1", "",temp[4,]), .before = support) %>%
    dplyr::select(-rules) %>%
    tibble::tibble()
  return(df_ap)
}


#' Calculate Disease Similarity Based on Gene Overlap
#'
#' This function calculates the similarity between diseases based on the overlap of associated genes.
#'
#' @param disease_data A data frame containing disease information.
#' @param gene_column The column name in `disease_data` that contains gene symbols.
#' @param target_column The column name in `disease_data` that contains the target or disease information.
#' @param genes An optional vector of genes to filter the data. If NULL, all genes in the data are used.
#' @param condition An optional condition to filter the data. Default is TRUE.
#'
#' @return A similarity matrix where each entry (i, j) represents the similarity between disease i and disease j.
#'
#' @examples
#' chemical_name <- paste0("Compound", 1:10)
#' gene_symbol <- paste0("GENE", 1:30)
#' disease_id <- paste0("Disease", 1:30)
#' compound_gene <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                             gene_symbol = sample(gene_symbol, 100, replace = TRUE),
#'                             organism_id = sample(c(9606, 10116), 100, replace = TRUE)) %>%
#'                             dplyr::distinct()
#' compound_disease <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                                disease_id = sample(disease_id, 100, replace = TRUE)) %>%
#'                                dplyr::distinct()
#' trans_data <- get_transaction(compound_gene, compound_disease)
#' aop_data <- get_aops(transaction = trans_data[[1]], genes = trans_data[[2]], diseases = trans_data[[3]])
#' disease_similarity(aop_data, gene_column = lhs, target_column = rhs)
#' @export
disease_similarity <- function(disease_data,
                               gene_column,
                               target_column,
                               genes = NULL,
                               condition = TRUE) {
  if (is.null(genes)) {
    genes <- unique(disease_data %>%
                      dplyr::select({{gene_column}}) %>%
                      dplyr::pull())
  }
  # gene_names <- disease_data %>%
  #   dplyr::select({{gene_column}}) %>%
  #   pull() %>%
  #   toupper()

  df_unique <- tibble::as_tibble(disease_data) %>%
    dplyr::filter(eval(parse(text = condition))) %>%
    dplyr::filter({{ gene_column }} %in%  toupper(genes)) %>%
    dplyr::distinct({{gene_column}}, {{target_column}}) %>%
    dplyr::mutate(Value = 1)
  df_wide <- df_unique %>%
    tidyr::pivot_wider(names_from = {{gene_column}}, values_from = Value, values_fill = list(Value = 0))
  disease_names <- dplyr::pull(df_wide[, 1])
  binary_matrix <- as.matrix(df_wide[, -1])
  rownames(binary_matrix) <- disease_names
  gene_counts <- rowSums(binary_matrix)
  similarity_matrix <- outer(1:nrow(binary_matrix), 1:nrow(binary_matrix), Vectorize(function(i, j) {
    if (i == j) {
      return(1)  # The similarity of a disease with itself is 1
    } else {
      # Calculate number of overlapped genes
      overlap <- sum(binary_matrix[i, ] & binary_matrix[j, ])
      # Calculate similarity
      similarity <- overlap / sqrt(gene_counts[i] * gene_counts[j])
      return(similarity)
    }
  }))
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- rownames(binary_matrix)
  return(similarity_matrix)
}




#' Calculate gene Similarity Based on Gene Overlap
#'
#' This function calculates the similarity between genes based on the overlap of associated diseases.
#'
#' @param link_data A data frame containing gene and disease information.
#' @param gene_column The column name in `disease_data` that contains gene symbols.
#' @param target_column The column name in `disease_data` that contains the target or disease information.
#' @param genes An optional vector of genes to filter the data. If NULL, all genes in the data are used.
#' @param condition An optional condition to filter the data. Default is TRUE.
#'
#' @return A similarity matrix where each entry (i, j) represents the similarity between disease i and disease j.
#'
#' @examples
#' chemical_name <- paste0("Compound", 1:10)
#' gene_symbol <- paste0("GENE", 1:30)
#' disease_id <- paste0("Disease", 1:30)
#' compound_gene <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                             gene_symbol = sample(gene_symbol, 100, replace = TRUE),
#'                             organism_id = sample(c(9606, 10116), 100, replace = TRUE)) %>%
#'                             dplyr::distinct()
#' compound_disease <- data.frame(chemical_name = sample(chemical_name, 100, replace = TRUE),
#'                                disease_id = sample(disease_id, 100, replace = TRUE)) %>%
#'                                dplyr::distinct()
#' trans_data <- get_transaction(compound_gene, compound_disease)
# aop_data <- get_aops(transaction = trans_data[[1]], genes = trans_data[[2]], diseases = trans_data[[3]])
# gene_similarity(aop_data, disease_column = rhs, target_column = lhs)
#' @export
gene_similarity <- function(gene_data,
                            disease_column,
                            target_column,
                            diseases = NULL,
                            condition = TRUE) {
  if (is.null(diseases)) {
    diseases <- unique(gene_data %>%
                         dplyr::select({{disease_column}}) %>%
                         dplyr::pull())
  }

  df_unique <- tibble::as_tibble(gene_data) %>%
    dplyr::filter(eval(parse(text = condition))) %>%
    dplyr::filter({{ disease_column }} %in% diseases) %>%
    dplyr::distinct({{disease_column}}, {{target_column}}) %>%
    dplyr::mutate(Value = 1)
  df_wide <- df_unique %>%
    tidyr::pivot_wider(names_from = {{disease_column}}, values_from = Value, values_fill = list(Value = 0))
  gene_names <- dplyr::pull(df_wide[, 1])
  binary_matrix <- as.matrix(df_wide[, -1])
  rownames(binary_matrix) <- gene_names
  disease_counts <- rowSums(binary_matrix)
  similarity_matrix <- outer(1:nrow(binary_matrix), 1:nrow(binary_matrix), Vectorize(function(i, j) {
    if (i == j) {
      return(1)  # The similarity of a gene with itself is 1
    } else {
      # Calculate number of overlapped diseases
      overlap <- sum(binary_matrix[i, ] & binary_matrix[j, ])
      # Calculate similarity
      similarity <- overlap / sqrt(disease_counts[i] * disease_counts[j])
      return(similarity)
    }
  }))
  rownames(similarity_matrix) <- colnames(similarity_matrix) <- rownames(binary_matrix)
  return(similarity_matrix)
}








#' Create Adverse Outcome Pathway (AOP) Network data
#'
#' @description
#' The `aop_network()` function constructs a network of Adverse Outcome Pathways (AOPs)
#' from the provided AOP data frame and disease vocabulary. It generates a network graph
#' with nodes representing diseases and genes, and edges representing the connections
#' between them, weighted by a specified column.
#'
#' @param aop_df A data frame containing AOP data with columns for the right-hand side (rhs)
#'   and left-hand side (lhs) of the associations.
#' @param dise_voc A data frame containing the disease vocabulary with columns `DiseaseID`, `DiseaseName`,
#'   and `DiseaseGroup`.
#' @param weight_column A string specifying the column in `aop_df` to use as the weight for the edges.
#'
#' @return A list with two components:
#' \item{edges}{A data frame representing the edges of the network, with columns `from`, `to`, and `weight`.}
#' \item{vertices}{A data frame representing the vertices of the network, with columns `name`, `nodes`, and `size`.}
#'
#' @examples
#' \dontrun{
#' # Sample data frames
#' aop_df <- data.frame(lhs = c("Gene1", "Gene2"), rhs = c("Disease1", "Disease2"), support = c(0.8, 0.6))
#' dise_voc <- data.frame(DiseaseID = c("Disease1", "Disease2"), DiseaseName = c("Disease A", "Disease B"), DiseaseGroup = c("Group1", "Group2"))
#'
#' # Create AOP network
#' aop_network_data <- aop_network(aop_df, dise_voc, "support")
#' edges <- aop_network_data$edges
#' vertices <- aop_network_data$vertices
#' }
#' @importFrom dplyr filter mutate select rename
#' @importFrom rlang enquo
#' @export
aop_network <- function(aop_df, dise_voc, weight_column) {
  dese_net <- aop_df %>%
    filter(rhs %in% dise_voc$DiseaseID) %>%
    mutate(Disease = dise_voc$DiseaseName[match(rhs, dise_voc$DiseaseID)],
           Group = dise_voc$DiseaseGroup[match(rhs, dise_voc$DiseaseID)], .before = support)
#    arrange(desc(oddsRatio))

  vert_df <- data.frame(name = c(unique(dese_net$Disease), unique(dese_net$lhs)),
                        nodes = rep(c("disease", "gene"),
                                    times = c(length(unique(dese_net$Disease)),
                                              length(unique(dese_net$lhs)))),
                        size= 10)
  vert_df$nodes[vert_df$nodes == "disease"] <- dise_voc$DiseaseGroup[match(vert_df$name[
    vert_df$nodes == "disease"], dise_voc$DiseaseName)]
  vert_df$name[vert_df$nodes == "gene"] <- block_fst(vert_df$name[vert_df$nodes == "gene"])
  dese <- vert_df$name[vert_df$nodes != "gene"]
  mac_desc <- match(dise_voc$DiseaseName,dese)
  idx <- mac_desc[!is.na(mac_desc)]
  vert_df[1:length(idx),] <- vert_df[idx,]

  edge_df <- dese_net %>%
    dplyr::select(Disease, lhs, {{weight_column}} ) %>%
    dplyr::rename(c("from" = Disease,"to" = lhs, "weight" = {{weight_column}} )) %>%
    dplyr::mutate(to = block_fst(to))
  aop_data <- list(edges = edge_df, vertices = vert_df)
  return(aop_data)
}


