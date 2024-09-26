#' Expand existing data by adding data for new compounds
#'
#' @description
#'
#' The `expand_data()` function adds data for extra compounds to the existing data.
#'
#' The function `expand_data()` is used to add data for one or more compounds to the
#' existing data. It processes new compound data and updates both `expression` and
#' `metadata` based on the database used. You can add compounds from both the
#' [TG-GATEs](https://toxico.nibiohn.go.jp/english/) and
#' [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) databases. The available
#' compounds in TG-GATEs and DrugMatrix can be found by calling the functions
#' [tggates_compounds()] and [dm_compounds()], respectively.
#'
#'
#' @param ge_matrix A matrix of gene expression data for given compounds.
#' @param metadata Metadata for given compounds, usually in the form of a data frame or
#'   tibble object.
#' @param database The database to be used, either `tggates` or `drugmatrix`.
#' @inheritParams get_tggates
#'
#' @family Data update helpers
#'
#' @returns A list.
#'   * A matrix of expanded gene expression data of compounds, probes in row and sample
#'   (CEL-ID) in column.
#'   * A data frame of expanded metadata about sample and phenotypic information of compounds.
#' @export
#'
#' @examples
#' \dontrun{
#' otg_data <- get_tggates(compounds = c("2-nitrofluorene", "N-nitrosomorpholine"))
#' otg_expndata <- expand_data(compounds = c("2-nitrofluorene", "N-nitrosomorpholine", "N-methyl-N-nitrosourea"),
#'                             ge_matrix = otg_data$expression,
#'                             metadata = otg_data$metadata)
#' head(otg_expndata$expression)
#' head(otg_expndata$metadata)
#' }
expand_data <- function(compounds,
                        ge_matrix,
                        metadata,
                        multicore = FALSE,
                        store = FALSE,
                        method = "rma",
                        output_dir = missing_arg(),
                        error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  is_metadata(metadata)
  is_compound(compounds)
  output_dir <- destination(output_dir)
  ext_comp <- unique(metadata$compound_name)
  ext_idx <- compounds %in% ext_comp
  if (any(ext_idx)) {
      cli::cli_warn(c(
        "Some of the compound already exist in the data",
        "i" = "{style_bold(col_red(compounds[ext_idx]))} {?is/are} already in the data,
      we ignore {?it/them} to add again."
      ))
  }
  new_dl <- compounds[!ext_idx]
  FC <- unique(metadata$fc)
  # if (length(FC) > 1) {
  #   cli_abort(
  #     c("Data must be unique type.",
  #       "x" = "You have used both `FC` and `normal` data.",
  #       "i" = "Please use only one type either `FC` or `normal` data.")
  #   )
  # }
  database <- unique(metadata$database)
  if (identical(database, "tggates")) {
    species = unique(metadata$species)
    data_type = gsub(" ", "_", unique(metadata$test_type))
    tissue =  unique(metadata$organ_id)
    dose_type =  unique(metadata$sin_rep_type)
    new_data <- get_tggates(compounds = new_dl,
                            species = species,
                            data_type = data_type,
                            tissue = tissue,
                            dose_type = dose_type,
                            fc = FC,
                            method = method,
                            multicore = multicore,
                            store = store,
                            output_dir = output_dir)
  }
  if (identical(database, "drugmatrix")) {
    tissue <- unique(metadata$organ_id)
    new_data <- get_drugmatrix(compounds = new_dl,
                               tissue = tissue,
                               fc = FC,
                               method = method,
                               multicore = multicore,
                               store = store,
                               output_dir = output_dir)
  }
  expr_mat <- cbind(ge_matrix, new_data$expression)
  dict_df <- dplyr::bind_rows(metadata, new_data$metadata)
  if (store) {
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "metadata.csv", sep = "/"))
  }
  tgx_data <- list(expr_mat, dict_df)
  names(tgx_data) <- c("expression", "metadata")
  return(tgx_data)
}


#' Shorten data by removing data for compounds
#'
#' @description
#'
#' The `shirink_data()` function is used to remove data for one or more compounds from the
#' existing data.
#'
#' The function `shrink_data()` is used to remove data for one or more compounds from the
#' existing data. It deletes corresponding samples of compounds from both the expression
#' and metadata tables, based on the database used.
#'
#' @inheritParams expand_data
#'
#' @returns A list.
#'   * A matrix of shortened gene expression values, with probes in rows and samples
#'   (CEL IDs) in columns.
#'   * A data frame of shortened metadata containing sample information for compounds.
#'
#' @export
#'
#' @examples
#' compounds <-  c("2NF", "NMOR", "MNU")
#' \dontrun{
#' otg_data <- get_tggates(compounds = c("2-nitrofluorene", "N-nitrosomorpholine", "N-methyl-N-nitrosourea"))
#' otg_srnkdata <- expand_data(compounds = c("2-nitrofluorene", "N-nitrosomorpholine"),
#'                            expr_data = otg_data$expression,
#'                            attr_data = otg_data$metadata,
#'                            database = "tggates")
#'
#' head(otg_srnkdata$expression)
#' head(otg_srnkdata$metadata)
#' }
shrink_data <- function(compounds,
                        ge_matrix,
                        metadata,
                        store = FALSE,
                        output_dir = missing_arg()) {
  test_data(ge_matrix, metadata)
  output_dir <- destination(output_dir)
  # test_input(database, c("tggates", "drugmatrix"))
  # is_compound(compounds)
  # if (identical(database, "tggates")) {
  #   rm_comp <- suppressMessages(abbr2name(compounds))
  # } else if (identical(database, "drugmatrix")) {
  #   rm_comp <- compounds
  # }
  ret_idx <- metadata$compound_name %in% compounds
  # can add cli_abort("data not found for compounds {compounds[retidx]})
  dict_df <- metadata[!ret_idx, ]
  act_bar <- match(dict_df$barcode, colnames(ge_matrix))
  expr_mat <- ge_matrix[, act_bar]
  if (store) {
    expr_mat %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    dict_df %>%
      readr::write_csv(file = paste(output_dir, "attribute.csv", sep = "/"))
  }
  srnk_data <- list(expr_mat, dict_df)
  names(srnk_data) <- c("expression", "metadata")
  return(srnk_data)
}


#' Clean and refresh the data
#'
#' @description
#'
#' The `update_data()` function updates the given gene expression and metadata by adding or
#' removing data for compounds based on queried compounds.
#'
#' The function `update_data()` is used to add or remove data for compounds in the existing
#' dataset. The compounds you wish to add must be available in the specified database.
#' The available compounds in [TG-GATEs](https://toxico.nibiohn.go.jp/english/) and
#' [DrugMatrix](https://ntp.niehs.nih.gov/data/drugmatrix) databases can be found by
#' calling the functions [tggates_compounds()] and [dm_compounds()], respectively.
#' Furthermore, the compounds you wish to remove must have both expression and attribute data.
#'
#' @param probes A vector of Affymetrix probe IDs to be included. The default is `NULL`,
#'   which includes all probes..
#' @inheritParams expand_data
#'
#' @seealso
#'   [expand_data()] for expanding the data by adding data for query compounds,
#'   [shrink_data()] for reducing the data by removing data for unwanted compounds.
#'
#' @return A list.
#'   * A matrix of gene update expression values, probes in row and sample (CEL ID) in column.
#'   * A data frame of update metadata about sample information of given compounds.
#' @export
#'
#' @examples
#' \dontrun{
#' otg_data <- get_tggates(compounds = c("2-nitrofluorene", "N-nitrosomorpholine"))
#' otg_cnrdata <- update_data(compounds = c("2-nitrofluorene", "N-methyl-N-nitrosourea"),
#'                         ge_matrix = otg_data$expression,
#'                         metadata = otg_data$metadata)
#' head(otg_cnrdata$expression)
#' head(otg_cnrdata$metadata)
#' }
update_data <- function(compounds,
                     ge_matrix,
                     metadata,
                     probes = NULL,
                     multicore = FALSE,
                     store = FALSE,
                     method = "rma",
                     output_dir = missing_arg(),
                     error_call = caller_env()) {
  test_data(ge_matrix, metadata)
  output_dir <- destination(output_dir)
  mis_comps <- as.vector(compounds[!(compounds %in% metadata$compound_name)])
  attr_comp <- unique(metadata$compound_name)
  ext_comps <- as.vector(attr_comp[!(attr_comp %in% compounds)])
  if (length(mis_comps) == 0 & length(ext_comps) == 0) {
    ge_matrix <- ge_matrix
    metadata <- metadata
  } else {
    # database  <- unique(metadata$database)
    if (!rlang::is_empty(ext_comps)) {
      cli::cli_alert_warning(c("The input data for the compound{?s}
                        {style_bold(col_red(backtick(ext_comps)))}
                        has been removed because {?this/these}
                        compound{?s} {?is/are} not listed in the query list."),wrap = TRUE)
      ext_data <- shrink_data(
        compounds = ext_comps,
        ge_matrix = ge_matrix,
        metadata = metadata,
        store = FALSE
      )
      ge_matrix <- ext_data$expression
      metadata <- ext_data$metadata %>%
        dplyr::mutate(dplyr::across(where(is.factor), droplevels))
    }
    if (!rlang::is_empty(mis_comps)) {
      exp_data <- expand_data(
        compounds = mis_comps,
        ge_matrix = ge_matrix,
        metadata = metadata,
        store = FALSE,
        multicore = multicore,
        output_dir = output_dir
      )
      ge_matrix <- exp_data$expression
      metadata <- exp_data$metadata
    }
  }
  attr_new <- metadata %>%
    dplyr::mutate(compound_name = factor(compound_name, levels = compounds)) %>%
    dplyr::arrange(compound_name)
  act_bar <- match(attr_new$barcode, colnames(ge_matrix))
  expr_new <- ge_matrix[, act_bar]
  if (is.null(probes)) {
    expr_new <- expr_new
  } else{
    true_probes <- probes %in% rownames(expr_new)
    if (any(!true_probes)) {
      false_probes <- probes[!true_probes]
      n_false <- length(false_probes)
      cli::cli_alert_warning(c("Given {.emph {n_false}} probe{?s} ",
                          "{style_bold(col_red(backtick(false_probes)))} ",
                          "{?is/are} not found in the gene expression data."),
                        wrap = TRUE)
    }
    expr_new <- expr_new[probes[true_probes],]
  }
  if (store) {
    expr_new %>%
      tibble::as_tibble(rownames = "probes") %>%
      readr::write_csv(file = paste(output_dir, "expression.csv", sep = "/"))
    attr_new %>%
      readr::write_csv(file = paste(output_dir, "metadata.csv", sep = "/"))
  }
  return(list(expression = expr_new, metadata = attr_new))
}

