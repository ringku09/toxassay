#' Compounds available in the Open TG-GATEs database
#'
#' @description
#' `tggates_compounds()` returns the compound names available in the Open TG-GATEs resource for a
#' specified combination of `species`, `data_type`, `tissue`, and (when applicable) `dose_type`.
#'
#' @section Open TG-GATEs:
#' Open TG-GATEs is a toxicogenomics resource containing transcriptomics profiles for a large set of
#' compounds (primarily pharmaceuticals) across *in vivo* and *in vitro* studies. The dataset includes
#' rat studies under single-dose and repeated-dose protocols and *in vitro* hepatocyte studies from
#' rat and human donors.
#'
#' @param species A character string specifying the species. Supported values are `"Rat"` and `"Human"`.
#' @param data_type A character string specifying the study type. Supported values are `"in_vivo"` and
#'   `"in_vitro"`.
#' @param tissue A character string specifying the tissue. Supported values are `"Liver"` and `"Kidney"`.
#' @param dose_type A character string specifying the dosing protocol for `data_type = "in_vivo"`.
#'   Supported values are `"Single"` and `"Repeat"`. Ignored for `data_type = "in_vitro"`.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @return A character vector of compound names available for the selected criteria.
#'
#' @examples
#' # Rat in vivo liver: single-dose
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#'
#' # Rat in vivo liver: repeated-dose
#' tggates_compounds(species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Repeat")
#'
#' # Rat in vitro liver
#' tggates_compounds(species = "Rat", data_type = "in_vitro", tissue = "Liver")
#'
#' # Human in vitro liver
#' tggates_compounds(species = "Human", data_type = "in_vitro", tissue = "Liver")
#'
#' @export
tggates_compounds <- function(species = c("Rat", "Human"),
                              data_type = c("in_vivo", "in_vitro"),
                              tissue = c("Liver", "Kidney"),
                              dose_type = c("Single", "Repeat"),
                              error_call = caller_env()) {
  species <- block_fst(species)
  data_type <- gsub("\\W", "_", data_type)
  tissue <- block_fst(tissue)
  dose_type <- block_fst(dose_type)
  test_input(species, auto_input = TRUE)
  test_input(data_type, auto_input = TRUE)
  test_input(tissue, auto_input = TRUE)
  test_input(dose_type, auto_input = TRUE)
  tggates_data <- compounds_tggates %>%
    dplyr::mutate(dplyr::across(- c("compound_name", "compound_abbr"), ~ ifelse(is.na(.), 0, 1))) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ . == 1))
  if (identical(data_type, "in_vivo")) {
    query_column <- glue::glue("{species}_{data_type}_{tissue}_{dose_type}")
  } else {
    query_column <- glue::glue("{species}_{data_type}_{tissue}")
  }
  if (! query_column %in% names(tggates_data)[-(1:2)]) {
    cli::cli_abort(c("Invalid data type parameter selection",
                     "x" = "Data type {style_bold(col_red(backtick(query_column)))} is \\
                not available in the open TG-GATEs database.",
                     "i" = " Please use one of the following valid parameter combinations: \\
                {style_italic(col_blue(backtick(names(tggates_data)[-(1:2)])))}.")
                   , call = error_call)
  }
  tgp_com <- tggates_data %>%
    dplyr::filter(!!rlang::sym(query_column)) %>%
    dplyr::select(compound_name) %>%
    dplyr::pull() %>%
    sort()
  #cli_alert(glue("Data available for the data type {style_bold(col_green(backtick(query_column)))} of compounds:"))
  return(tgp_com)
}

#' Compounds available in the DrugMatrix database
#'
#' @description
#' `drugmatrix_compounds()` returns compound names available in the DrugMatrix resource for a given
#' `tissue`.
#'
#' @section DrugMatrix tissues:
#' DrugMatrix contains toxicogenomics profiles for a large set of compounds measured in rat studies.
#' Available tissues include:
#' \itemize{
#'   \item `Liver`
#'   \item `Kidney`
#'   \item `Heart`
#'   \item `Hepatocytes` (rat primary hepatocytes; in vitro)
#' }
#'
#' @param tissue A character string specifying the tissue. Supported values are `"Liver"`, `"Kidney"`,
#'   `"Heart"`, and `"Hepatocytes"`.
#'
#' @return A character vector of compound names with available data for the selected `tissue`.
#'
#' @examples
#' # Compounds with data available for liver
#' drugmatrix_compounds(tissue = "Liver")
#'
#' # Compounds with data available for kidney
#' drugmatrix_compounds(tissue = "Kidney")
#'
#' # Compounds with data available for heart
#' drugmatrix_compounds(tissue = "Heart")
#'
#' # Compounds with data available for hepatocytes (in vitro)
#' drugmatrix_compounds(tissue = "Hepatocytes")
#'
#' @export
drugmatrix_compounds <- function(tissue = c("Liver", "Kidney", "Heart", "Hepatocytes")) {
  dm_comps <- dm_metadata %>%
    dplyr::filter(Compound != "Control") %>%
    dplyr::select(c(Compound, Tissue)) %>%
    dplyr::distinct()
  dm_comps$values <- 1
  dm_data <- dm_comps %>%
    tidyr::pivot_wider(names_from = Tissue, values_from = values, values_fill = 0) %>%
    dplyr::rename("Hepatocytes (in vitro)" = Hepatocytes) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ . == 1))
  tissue <- block_fst(tissue)
  test_input(tissue, auto_input = TRUE)
  if (identical(tissue, "Liver")) {
    dm_com <- dm_data %>%
      dplyr::filter(Liver) %>%
      dplyr::select(Compound) %>%
      dplyr::pull( ) %>%
      sort()
  } else if (identical(tissue, "Kidney")) {
    dm_com <- dm_data %>%
      dplyr::filter(Kidney) %>%
      dplyr::select(Compound) %>%
      dplyr::pull() %>%
      sort()
  } else if (identical(tissue, "Heart")) {
    dm_com <- dm_data %>%
      dplyr::filter(Heart) %>%
      dplyr::select(Compound) %>%
      dplyr::pull() %>%
      sort()
  } else if (identical(tissue, "Hepatocytes")) {
    dm_com <- dm_data %>%
      dplyr::filter(`Hepatocytes (in vitro)`) %>%
      dplyr::select(Compound) %>%
      dplyr::pull() %>%
      sort()
  }
  #cli_alert(glue("Data available for the {style_bold(col_green(backtick(tissue)))} tissue of compounds:"))
  return(dm_com)
}

#' Validate compound names against TG-GATEs and DrugMatrix databases
#'
#' @description
#' `validate_compounds()` checks whether the provided `compounds` are available in the selected database.
#' When `database` is missing, the function validates compound names against both Open TG-GATEs and
#' DrugMatrix. When `database` is provided, the function also validates the corresponding selection
#' parameters (e.g., `species`, `data_type`, `tissue`, and `dose_type`) and checks compound
#' availability within that subset.
#'
#' @param compounds A character vector of compound names to validate.
#' @param database Database to validate against. Must be `"tggates"` or `"drugmatrix"`. If missing,
#'   compounds are checked against both databases.
#' @param species Species used to subset Open TG-GATEs data when `database = "tggates"`. Must be
#'   `"Rat"` or `"Human"`.
#' @param data_type Data type used to subset Open TG-GATEs when `database = "tggates"`. Must be
#'   `"in_vivo"` or `"in_vitro"`.
#' @param tissue Tissue used to subset the selected database. For `database = "tggates"`, valid values
#'   are `"Liver"` and `"Kidney"`. For `database = "drugmatrix"`, valid values are `"Liver"`, `"Kidney"`,
#'   `"Heart"`, and `"Hepatocytes"`.
#' @param dose_type Dose type used to subset Open TG-GATEs when `database = "tggates"`. Must be
#'   `"Single"` or `"Repeat"`.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It returns `TRUE` invisibly if all compound
#' names are valid for the requested database/subset; otherwise it signals an error.
#'
#' @examples
#' # Check compounds against a specific Open TG-GATEs subset
#' validate_compounds(
#'   compounds = c("aspirin", "acetaminophen"),
#'   database = "tggates",
#'   species = "Rat",
#'   data_type = "in_vivo",
#'   tissue = "Liver",
#'   dose_type = "Single")
#'
#' \dontrun{
#' # Check compounds against both Open TG-GATEs and DrugMatrix (when database is missing)
#' validate_compounds(compounds = c("aspirin", "paracetamol"))
#' }
#'
#' @export
validate_compounds <- function(compounds,
                               database = rlang::missing_arg(),
                               species = rlang::missing_arg(),
                               data_type = rlang::missing_arg(),
                               tissue = rlang::missing_arg(),
                               dose_type = rlang::missing_arg(),
                               error_call = rlang::caller_env()) {

  if (rlang::is_empty(compounds)) {
    cli::cli_abort(c("{.var comp_name} must be non-empty.",
                     "i" = "You have supplied an empty vector, please provide compound(s) name instred.")
                   , call = error_call)
  }
  if (rlang::is_missing(database)) {
    comp_tg <- compounds_tggates
    comp_dm <- dm_metadata
    available_com <- unique(c(comp_tg$compound_name, comp_dm$Compound))
    comp_is <- compounds %in% comp_tg$compound_name
    if (any(!comp_is)) {
      avail_tgp <- rlang::englue("tggates_compounds()")
      avail_dm <- rlang::englue("dm_compounds()")
      cli::cli_abort(c("Invalid compound name.",
                       "x" = "{style_bold(col_red(backtick(compounds[!comp_is])))} compound{?s} \\
                {?is/are} not available in open TG-GATEs and DrugMatrix database.",
                       "i" = "Please find the name of available compound by calling the function \\
                {style_italic(col_blue(backtick(avail_tgp)))} for TG-GATEs and \\
                {style_italic(col_blue(backtick(avail_dm)))} for DrugMatrix database.")
                     , call = error_call)
    }
  } else {
    test_input(database, c("tggates", "drugmatrix"))
    if (identical(database, "tggates")) {
      test_input(species, c("Rat", "Human"))
      test_input(data_type, c("in_vivo", "in_vitro"))
      test_input(tissue, c("Liver", "Kidney"))
      test_input(dose_type, c("Single", "Repeat"))
      comp_tggates <- tggates_compounds(species = species,
                                        data_type = data_type,
                                        tissue = tissue,
                                        dose_type = dose_type)
      comp_is <- compounds %in% comp_tggates
      if (any(!comp_is)) {
        avail_com <- rlang::englue("tggates_compounds()")
        cli::cli_abort(c("Compound must available in open TG-GATEs database.",
                         "x" = "{style_bold(col_red(backtick(compounds[!comp_is])))} compound name{?s} \\
                {?is/are} not available in open TG-GATEs database.",
                         "i" = "Please find the name of available compound by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
                       , call = error_call)
      }
    } else if (identical(database, "drugmatrix")) {
      test_input(tissue, c("Liver", "Kidney", "Heart", "Hepatocytes"))
      comp_dm <- drugmatrix_compounds(tissue = tissue)
      comp_is <- compounds %in% comp_dm
      if (any(!comp_is)) {
        avail_com <- rlang::englue(" drugmatrix_compounds()")
        cli::cli_abort(c("Compound must available in DrugMatrix database.",
                         "x" = "{style_bold(col_red(backtick(compounds[!comp_is])))} compound name{?s} \\
                {?is/are} not available in DrugMatrix database.",
                         "i" = "Please find the name of available compound by calling the function {style_italic(col_blue(backtick(avail_com)))}.")
                       , call = error_call)
      }
    }
  }
}
