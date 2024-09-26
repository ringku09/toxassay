#' Validate and Auto-Correct Function Inputs
#'
#' This function validates the input provided by the user against a list of expected inputs.
#' It can automatically correct the input if `auto_input` is set to `TRUE`,
#' or it will prompt the user with an error message if the input is invalid.
#'
#' @param input The input provided by the user to be validated.
#' @param inputs A vector of valid inputs to compare against. If `inputs` is missing, the function will attempt to retrieve the default inputs from the parent environment.
#' @param auto_input Logical value indicating whether to automatically correct the input to the first valid entry in `inputs` if multiple or invalid inputs are provided. Default is `FALSE`.
#' @param error_call The environment in which to evaluate the error call. Default is `caller_env()`.
#'
#' @return This function does not return a value. It either validates the input or throws an error if the input is invalid.
#'
#' @examples
#' # Assume valid inputs are "A", "B", and "C"
#' valid_inputs <- c("A", "B", "C")
#'
#' # Example 1: Correct input
#' test_input(input = "A", inputs = valid_inputs)
#'
#' # Example 2: Auto-correct input (if auto_input is TRUE)
#' input_val = c("A", "D")
#' test_input(input = input_val, inputs = valid_inputs, auto_input = TRUE)
#'
#' # Example 3: Handle missing input
#' input_val = NULL
#' test_input(input = input_val, inputs = valid_inputs)
#'
#' @export
test_input <- function(input, inputs, auto_input = FALSE, error_call = caller_env()) {
  if (rlang::is_missing(inputs)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    inputs <- eval(formal.args[[as.character(substitute(input))]],
                   envir = sys.frame(sysP))
  }
  size <- length(inputs)
  if (size > 9L) {
    new_inputs <- c(inputs[1:9], "...", inputs[length(inputs)])
  } else {
    new_inputs <- inputs
  }
  if (rlang::is_missing(input) | rlang::is_empty(input)) {
    cli::cli_abort(c("{backtick(deparse(substitute(input)))} is missing.",
                "i" = "Please choose either {add_or(style_bold(col_green(backtick(new_inputs))))} instead."),
              wrap = TRUE, call = error_call)
  }
  if (length(input) > 1L) {
    if (auto_input) {
      eval.parent(substitute(input <- inputs[1]))
    } else {
      cli::cli_abort(c(
        "Multiple {backtick(deparse(substitute(input)))} are not allowed when `auto_input = FALSE`",
        "x" = "You have provided {backtick(deparse(substitute(input)))} of {style_bold(col_red(backtick(input)))}.",
        "i" = "Please choose either {add_or(style_bold(col_green(backtick(input))))} instead."
      ), wrap = TRUE, call = error_call)
    }
  }
  if (length(input) == 0L) {
    cli::cli_abort(c("Input {backtick(deparse(substitute(input)))} must have at least 1 location.",
                "i" = "Please choose either {add_or(style_bold(col_green(backtick(new_inputs))))} instead."),
              wrap = TRUE, call = error_call)
  }
  if (any(!input %in% inputs)) {
    cli::cli_abort(c(
      "Wrong {backtick(deparse(substitute(input)))}",
      "x" = "The {backtick(deparse(substitute(input)))} {style_bold(col_red(backtick(input)))} you provided is not recognized.",
      "i" = "Please choose either {add_or(style_bold(col_green(backtick(new_inputs))))} instead."
    ), wrap = TRUE, call = error_call)
  }
}


#' Validate Expression and Attribute Data
#'
#' This function checks the validity of the expression data and attribute data provided by the user.
#' It ensures that the expression data is a matrix and that the attribute data is a data frame.
#' Additionally, it verifies that the column names in the expression data match the barcodes specified in the attribute data.
#'
#' @param expr_data A matrix or array containing expression data. The matrix should have `probes` as rows and `barcode` as columns.
#' @param attr_data A data frame or tibble containing the attribute data, with required columns including `barcode`, `compound_name`, `dose_level`, and `time_level`.
#'
#' @return This function does not return a value. It is used for validation and will throw an error if the data is not valid.
#'
#' @details
#' - If `expr_data` is not a matrix or array, an error is thrown indicating that the expression data must be a matrix.
#' - If `attr_data` is not a data frame, an error is thrown indicating that the attribute data must be a data frame.
#' - The function also checks that all column names in `expr_data` are present in the `barcode` column of `attr_data`. If any are missing, an error is thrown.
#'
#' @examples
#' tgx_data <- simulate_tgxdata(n_de = 20, n_ee = 10, n_com = c(3, 3))
#' test_data(tgx_data$expression, tgx_data$metadata)
#' \dontrun{
#' # Example with invalid data
#' expr_data_invalid <- array(rnorm(100), dim = c(10, 10, 1))
#' attr_data_invalid <- list(barcode = paste0("sample", 1:10))
#' test_data(expr_data_invalid, attr_data_invalid) # This will throw an error
#' }
#'
#' @export
test_data <- function(expr_data, attr_data) {
  if (!inherits(expr_data, c("matrix", "array"))) {
    cli_abort(c("The expression data must be a matrix.",
                "x" = "The class {style_bold(col_cyan(backtick(class(expr_data))))} of \\
                expression data {?is/are} not supported.",
                "i" = "Please make expression data as matrix \\
             (`probes` in rows and `barcode` in columns)."))

  }
  if (!rlang::inherits_any(attr_data, c("tbl_df", "tbl", "data.frame"))) {
    cli::cli_abort(c("The attribute data must be a data frame.",
                "x" = "The class {style_bold(col_cyan(backtick(class(attr_data))))} of \\
                attribute data {?is/are} not supported.",
                "i" = "Please make attribute data as data frame."))

  }
  test_column(c("barcode", "compound_name", "dose_level", "time_level"), attr_data)
  match_bar <- colnames(expr_data) %in% attr_data$barcode
  if (any(!match_bar)) {
    cli::cli_abort(c("All columns/samples in expression data must given in metadata.",
                "x" = "{style_bold(col_cyan(backtick(colnames(expr_data)[match_bar])))} column{?s} \\
                {?is/are} not described in the metadata.",
                "i" = "Please remove unspecified columns from the expression data. After that you \\
               can update your data by calling the function `update_data()`"))
  }
}


#' Validate and Group Compounds
#'
#' This function checks and organizes compound groups provided as input arguments.
#' It ensures that the inputs are of appropriate types (character, factor, or list) and
#' provides informative error messages if the inputs are not valid.
#'
#' @param ... Compound groups supplied as vectors, lists, or factors containing compound names or abbreviations.
#' If a single list is provided, it will be used directly. Otherwise, multiple vectors or lists can be supplied.
#' @param error_call The environment in which the error is raised. Default is the caller environment (`caller_env()`).
#'
#' @return A named list of compound groups. If names are not provided, they will be automatically assigned as "Group_A", "Group_B", etc.
#'
#' @examples
#' # Example with valid input
#' group1 <- c("compound1", "compound2")
#' group2 <- c("compound3", "compound4")
#' test_group(group1, group2)
#'
#' # Example with a list
#' groups <- list(group1 = c("compound1", "compound2"), group2 = c("compound3", "compound4"))
#' test_group(groups)
#'
#' # Example with mixed input types (this will trigger an error)
#' \dontrun{
#' # test_group(group1, group2, list("compound5"))
#' }
#' @export
test_group <- function(..., error_call = rlang::caller_env()) {
  comps_gr <- rlang::list2(...)
  if (length(comps_gr) == 1 && rlang::is_bare_list(comps_gr[[1]])) {
    comps_gr <- comps_gr[[1]]
}
  # if (!length(comps_gr)>1) {
  #   cli_abort(c("The number of compound groups must be greater than one.",
  #               "i" = "Please use multiple compound groups for comparison."),
  #             wrap =TRUE, call = error_call)
  # }
  arg_checker <- unlist(lapply(comps_gr, function(x) inherits(x, c("character", "factor", "list"))))
  right_arg <- comps_gr[arg_checker]
  right_class <- unlist(lapply(right_arg, class))
  if (length(unique(right_class))>1) {
    right_name <- sapply(substitute(list(...))[-1], as.character)[arg_checker]
    cli::cli_abort(c("Invalid inpute of compound groups.",
                "x" = "You have supplied the object{?s} {paste(style_bold(col_red(right_name)), style_bold(col_blue(right_class)), sep = ' as a ')}
                class{?, respectively}.",
                "i" = "Please provide the compound groups as either a list or multiple vectors containing compound names or abbreviations."),
              wrap = TRUE, call = error_call)
  }
  # wrong_arg <- which(!arg_checker)
  wrong_arg <- sapply(substitute(list(...))[-1], as.character)[!arg_checker]
  if (!all(arg_checker)) {
    arg_nm <- paste("argument_name", LETTERS[1:length(wrong_arg)])
    cli::cli_abort(c("Invalid use of the argument in the function.",
                "x" = "You have supplied object{?s} {style_bold(col_red(wrong_arg))}  without assigning
                {?it/them} as an argument{?s}. However, {?it/they} {?does/may} not contain any compound names.",
                "i" = "Please provide the argument explicitly by supplying the compound names correctly.
                Alternatively, if applicable, you can assign the appropiate argument{?s} in the function as:
                {paste(style_bold(col_green('argument_name')), style_bold(col_red(wrong_arg)), sep = ' = ')}."),
              wrap = TRUE, call = error_call)
  }
  if (is.null(names(comps_gr))) {
    names(comps_gr) <- paste("Group", LETTERS[1:length(comps_gr)], sep = "_")
  }
  return(comps_gr)
}

#' Validate Input Against a List of Elements
#'
#' This function checks whether a given `input` is present in a specified list of `elements`.
#' If the input is missing, empty, or contains elements not found in the list, the function
#' throws an informative error message using `cli_abort`.
#'
#' @param input A vector representing the input to be validated against the `elements`.
#' @param elements A vector of valid elements that `input` should be compared against.
#' @param error_call The environment to be used in the `cli_abort` call for error reporting. Default is `caller_env()`.
#'
#' @return This function does not return a value. It is used for validation and will throw an error if validation fails.
#'
#' @examples
#' valid_elements <- c("apple", "banana", "cherry")
#' test_element("apple", valid_elements)  # No error
#' \dontrun{
#' test_element("orange", valid_elements) # Throws an error
#' }
#' @export
test_element <- function(input, elements, error_call = caller_env()) {
  if (length(elements) > 9L) {
    new_inputs <- c(elements[1:9], "...", elements[length(elements)])
  } else {
    new_inputs <- elements
  }
  if (rlang::is_missing(input) | rlang::is_empty(input)) {
    cli::cli_abort(c("{backtick(deparse(substitute(input)))} is missing.",
                "i" = "Please choose either {add_or(style_bold(col_green(backtick(new_inputs))))} instead."),
              wrap = TRUE, call = error_call)
  }
  # if (length(input) == size) {
  #   cli_abort(c("{backtick(deparse(substitute(input)))} is missing.",
  #               "i" = "Please choose either {add_or(style_bold(col_green(backtick(inputs))))} instead."),
  #             wrap = TRUE, call = error_call)
  # }

  if (length(input) == 0L) {
    cli::cli_abort(c("Input {backtick(deparse(substitute(input)))} must have at least 1 location.",
                "i" = "Please choose either {add_or(style_bold(col_green(backtick(new_inputs))))} instead."),
              wrap = TRUE, call = error_call)
  }

  exist_idx <- input %in% elements
  miss_element <- input[!exist_idx]
  if (length(miss_element) > 9L) {
    new_element <- c(miss_element[1:9], "...", miss_element[length(miss_element)])
  } else {
    new_element <- miss_element
  }
  if (any(!exist_idx)) {
    cli::cli_abort(c("Input {backtick(deparse(substitute(input)))} must available in the \\
                {backtick(deparse(substitute(elements)))}.",
                "x" = "{style_bold(col_red(backtick(new_element)))} element{?s} \\
                {?is\are} not available in {backtick(deparse(substitute(elements)))}.",
                "i" = "Please select element from {style_italic(col_blue(backtick(new_inputs)))}."))
  }
}

#' Validate Column Presence in Data Frame or List
#'
#' This function checks whether a specified column exists in a provided data frame, tibble, or list. It provides informative error messages if the input does not meet the expected criteria.
#'
#' @param column A character vector representing the name(s) of the column(s) to be checked.
#' @param df A data frame, tibble, or list in which the presence of the specified column(s) will be tested.
#'
#' @return The function does not return a value but throws an error if the column does not exist in the provided data structure.
#'
#' @details
#' - If `df` is not a data frame, tibble, or list, an error is thrown.
#' - If the `column` argument is empty or the specified column(s) do not exist in `df`, the function provides a detailed error message indicating the issue and suggesting possible corrections.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
#' test_column("a", df) # No error
#' \dontrun{
#' test_column("d", df) # Throws an error: "d" does not exist in `df`
#' }
#' @export
test_column <- function(column, df) {
  if (! (is.data.frame(df) | tibble::is_tibble(df) || is.list(df))) {
    cli::cli_abort(c("{.arg df} must be a `data.frame` or `tibble` or `list`.",
                "i" = "Please provide the correct {.arg df}."))
  }
  columns <- names(df)
  if (length(columns) > 9L) {
    new_columns <- c(columns[1:9], "...", columns[length(columns)])
  } else {
    new_columns <- columns
  }
  if (length(column) == 0) {
    cli::cli_abort(c("Input {gsub('_', ' ', deparse(substitute(column)))} must have at least 1 location.",
                "i" = "Please select column from {style_italic(col_blue(backtick(new_columns)))} \\
                ,or provide the correct {backtick(deparse(substitute(df)))}."))
  }
  # if (length(column) > 1) {
  #   cli_abort(c(
  #     "Multiple {backtick(gsub('_', ' ', deparse(substitute(columns))))} are not allowed",
  #     "x" = "You have provided {backtick(gsub('_', ' ', deparse(substitute(columns))))} {style_bold(col_red(backtick(column)))}.",
  #     "i" = "Please choose either {add_or(style_bold(col_green(backtick(column))))} instead."
  #   ))
  # }
  #
  exist_idx <- column %in% columns
  miss_col <- column[!exist_idx]
  if (length(miss_col) > 9L) {
    new_column <- c(miss_col[1:9], "...", miss_col[length(miss_col)])
  } else {
    new_column <- miss_col
  }
  if (any(!exist_idx)) {
    cli::cli_abort(c("Input {.arg column} must available in the {backtick(deparse(substitute(df)))}.",
                "x" = "{style_bold(col_red(backtick(new_column)))} column{?s} \\
                {?is\are} not available in {backtick(deparse(substitute(df)))}.",
                "i" = "Please select column from {style_italic(col_blue(backtick(new_columns)))} \\
                ,or provide the correct {backtick(deparse(substitute(df)))}."))
  }
}

#' Validate Data Structure for ToxAssay
#'
#' This function checks whether the provided data structure is a valid `ToxAssay` object or a list with the correct format.
#' It ensures that the input is either a `ToxAssay` object or a list of length 5, containing the elements `group`, `compound`, `dose`, `time`, and `replication`.
#'
#' @param expr_str An object to be validated. It should be either an object of class `ToxAssay` or a list.
#'
#' @return If the input is valid, the function returns `TRUE`. If the input is invalid, an error message is triggered.
#'
#' @examples
#' # Example of a valid ToxAssay object or list
#' valid_list <- list(group = "A", compound = "Methimazole", dose = "High", time = "24h", replication = 3)
#' test_datastr(valid_list)
#'
#' # Example of an invalid list (less than 5 elements)
#' invalid_list <- list(group = "A", compound = "Methimazole")
#' # This will trigger an error
#' \dontrun{
#' test_datastr(invalid_list)
#' }
#' @export
test_datastr <- function(expr_str) {
  if (!is.list(expr_str) & !inherits(expr_str, "ToxAssay")) {
    cli::cli_abort(c("{.var expr_str} must be object of class `ToxAssay` or `list`.",
                "x" = "The class {style_bold(col_cyan(backtick(class(expr_str))))} of \\
                {.var expr_str} is not supported.",
                "i" = "Please provide {.var expr_str} as a list of size 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."))
  } else if (length(expr_str)  !=  5) {
    cli::cli_abort(c("The length of {.var expr_str} must be 5.",
                "i" = "You have supplied a list {.var expr_str} of size {length(expr_str)}, please
                make sure {.var expr_str} has a length of 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."))
  }
}










#' Validate Compound Names Against Databases
#'
#' This function checks whether the provided compound names exist in the specified database
#' (TG-GATEs or DrugMatrix) and validates additional parameters such as species, data type,
#' tissue, and dose type. If no database is specified, it checks the compounds against both
#' the TG-GATEs and DrugMatrix databases.
#'
#' @param compounds A character vector of compound names to be checked.
#' @param database The name of the database to check against. Should be either `"tggates"` or `"drugmatrix"`. If not specified, both databases are checked.
#' @param species The species to filter for when using the `"tggates"` database. Should be `"Rat"` or `"Human"`.
#' @param data_type The type of data to filter for when using the `"tggates"` database. Should be `"in_vivo"` or `"in_vitro"`.
#' @param tissue The tissue type to filter for. For `"tggates"`, it should be `"Liver"` or `"Kidney"`. For `"drugmatrix"`, it should be `"Liver"`, `"Kidney"`, `"Heart"`, or `"Hepatocytes"`.
#' @param dose_type The dose type to filter for when using the `"tggates"` database. Should be `"Single"` or `"Repeat"`.
#' @param error_call Environment used to capture errors. Default is `caller_env()`.
#'
#' @return Returns `TRUE` if all compound names are valid, otherwise throws an error.
#'
#' @examples
#' # Check if a compound is in the TG-GATEs or DrugMatrix database
#' is_compound(compounds = c("aspirin", "acetaminophen"), database = "tggates", species = "Rat", data_type = "in_vivo", tissue = "Liver", dose_type = "Single")
#'\dontrun{
#' # Check if a compound exists in either the TG-GATEs or DrugMatrix database
#' is_compound(compounds = c("aspirin", "paracetamol"))
#' }
#' @export
is_compound <- function(compounds,
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

#' Validate Metadata for Toxicogenomics Datasets
#'
#' This function checks the integrity and consistency of metadata provided for toxicology datasets,
#' specifically for `tggates` and `drugmatrix` databases. It ensures that the metadata contains
#' the necessary columns, and that certain attributes (e.g., `database`, `fc`, `species`, `tissue`)
#' are consistent within the dataset.
#'
#' @param metadata A data frame containing metadata for the toxicogenomics dataset.
#'  The metadata must include specific columns such as `compound_name`, `dose_level`, `time_level`, `organ_id`, `database`, `arr_design`, and `fc`.
#'
#' @details The function performs the following checks:
#' - Ensures that the metadata contains required columns.
#' - Verifies that only one `database` (either `tggates` or `drugmatrix`) is present in the metadata.
#' - Checks that the `fc` (fold change) data type is consistent (i.e., only `FC` or `normal` data).
#' - For `tggates` database: Ensures metadata consistency for `species`, `test_type`, `organ_id`, and `sin_rep_type`.
#' - For `drugmatrix` database: Ensures metadata consistency for `organ_id`.
#'
#'
#' @return This function does not return a value. It is used for validation purposes and will raise an error if the metadata is inconsistent.
#'
#' @examples
#' # Example usage:
#' metadata <- data.frame(
#'   compound_name = c("Compound A", "Compound B"),
#'   dose_level = c("High", "Low"),
#'   time_level = c("24h", "48h"),
#'   organ_id = "Liver",
#'   database = "drugmatrix",
#'   arr_design = c("Design1", "Design2"),
#'   fc = "FC"
#' )
#' is_metadata(metadata)
#'
#' @export
is_metadata <- function(metadata) {
  test_column(c("compound_name", "dose_level", "time_level", "organ_id", "database", "arr_design","fc"), metadata)
  database <- unique(metadata$database)
  test_input(database, c("tggates", "drugmatrix"))
  if (length(database) > 1) {
    cli::cli_abort(c("More than one `database` not allowed.",
                "x" = "{style_bold(col_red(backtick(metadata)))} database {?is/are}  \\
                present in the metadata.",
                "i" = "Please use metadata with only one `database` for analysis.")
              , call = error_call)
  }
  FC <- unique(metadata$fc)
  if (length(FC) > 1) {
    cli::cli_abort(
      c("Data must be unique type.",
        "x" = "You have used both `FC` and `normal` data.",
        "i" = "Please use only one type either `FC` or `normal` data.")
    )
  }
  if (identical(database, "tggates")) {
    test_column(c("species", "test_type", "sin_rep_type"), metadata)
    species = unique(metadata$species)
    data_type = gsub(" ", "_", unique(metadata$test_type))
    tissue =  unique(metadata$organ_id)
    dose_type =  unique(metadata$sin_rep_type)
    if (length(species) > 1) {
      cli::cli_abort(c("More than one `species` not allowed.",
                  "x" = "{style_bold(col_red(backtick(species)))} species {?is/are}  \\
                present in the metadata.",
                  "i" = "Please use metadata with only one `species` for subsequent analysis.")
                , call = error_call)
    }
    if (length(data_type) > 1) {
      cli::cli_abort(c("More than one `data_type` not allowed.",
                  "x" = "{style_bold(col_red(backtick(data_type)))} test type{?s} {?is/are}  \\
                present in the metadata.",
                  "i" = "Please use metadata with only one `data_type` for subsequent analysis.")
                , call = error_call)
    }
    if (length(tissue) > 1) {
      cli::cli_abort(c("More than one `tissue` not allowed.",
                  "x" = "{style_bold(col_red(backtick(tissue)))} tissue{?s} {?is/are}  \\
                present in the metadata.",
                  "i" = "Please use metadata with only one `tissue` for subsequent analysis.")
                , call = error_call)
    }
    if (length(dose_type) > 1) {
      cli::cli_abort(c("More than one `dose_type` not allowed.",
                  "x" = "{style_bold(col_red(backtick(dose_type)))} experiment type{?s} {?is/are} \\
                present in the attribute data.",
                  "i" = "Please use attribute data with only one `dose_type` for subsequent analysis.")
                , call = error_call)
    }
  }
  if (identical(database, "drugmatrix")) {
    tissue =  unique(metadata$organ_id)
    if (length(tissue) > 1) {
      cli::cli_abort(c("More than one `tissue` not allowed.",
                  "x" = "{style_bold(col_red(backtick(tissue)))} tissue{?s} {?is/are}  \\
                present in the metadata.",
                  "i" = "Please use metadata with only one `tissue` for subsequent analysis.")
                , call = error_call)
    }
  }
}
