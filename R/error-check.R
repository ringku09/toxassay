#' Validate and optionally auto-correct function input
#'
#' @description
#' `test_input()` checks whether a user-supplied argument matches a set of allowed values. If the
#' input is missing, empty, invalid, or contains multiple values, the function raises an informative
#' error. When `auto_input = TRUE`, invalid or multiple inputs are automatically replaced with the
#' first valid value in `inputs`.
#'
#' @param input The value provided by the user to be validated.
#' @param inputs A character vector of allowed values. If missing, the function attempts to retrieve
#'   the default valid values from the parent function’s formal arguments.
#' @param auto_input Logical indicating whether to automatically replace invalid or multiple inputs
#'   with the first element of `inputs`. Default is `FALSE`.
#' @param error_call The environment used for error reporting. Default is `caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It either validates (and possibly modifies)
#' `input` or signals an error; no value is returned.
#'
#' @examples
#' # Validate input against a set of allowed values
#' valid_inputs <- c("A", "B", "C")
#' test_input(input = "A", inputs = valid_inputs)
#'
#' \dontrun{
#' # Automatically correct invalid input
#' test_input(input = "X", inputs = valid_inputs, auto_input = TRUE)
#' }
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

#' Validate required column names in a data object
#'
#' @description
#' `test_column()` checks whether one or more specified column names are present in a given data
#' object. The function supports data frames, tibbles, and named lists, and provides informative
#' error messages when required columns are missing or the input object is of an unsupported type.
#'
#' @param column A character vector specifying the name(s) of the column(s) to validate.
#' @param df A data frame, tibble, or named list in which the presence of `column` is checked.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It returns `NULL` invisibly if all specified
#' columns are present; otherwise it signals an error.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item `df` is a data frame, tibble, or list.
#'   \item `column` is non-empty.
#'   \item All elements of `column` are present in `names(df)`.
#' }
#' If columns are missing, the error message lists unavailable columns and suggests valid
#' alternatives when possible.
#'
#' @examples
#' df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
#' test_column("a", df)  # Valid column
#'
#' \dontrun{
#' # Missing column triggers an error
#' test_column("d", df)
#' }
#'
#' @export
test_column <- function(column, df, error_call = rlang::caller_env()) {
  if (! (is.data.frame(df) | tibble::is_tibble(df) || is.list(df))) {
    cli::cli_abort(c("{.arg df} must be a `data.frame` or `tibble` or `list`.",
                     "i" = "Please provide the correct {.arg df}."), call = error_call)
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
                ,or provide the correct {backtick(deparse(substitute(df)))}."), call = error_call)
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
                ,or provide the correct {backtick(deparse(substitute(df)))}."), call = error_call)
  }
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



#' Validate input values against a reference set
#'
#' @description
#' `test_element()` verifies that one or more input values are present in a predefined set of valid
#' elements. If the input is missing, empty, or contains values not included in `elements`, the
#' function raises an informative error message.
#'
#' @param input A vector of values to be validated.
#' @param elements A vector defining the set of valid values against which `input` is checked.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It returns `NULL` invisibly if validation
#' succeeds; otherwise it signals an error.
#'
#' @examples
#' valid_elements <- c("apple", "banana", "cherry")
#' test_element("apple", valid_elements)
#'
#' \dontrun{
#' # Invalid input triggers an error
#' test_element("orange", valid_elements)
#' }
#'
#' @export
test_element <- function(input, elements, error_call = rlang::caller_env()) {
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
                     "i" = "Please select element from {style_italic(col_blue(backtick(new_inputs)))}."), call = error_call)
  }
}

#' Validate perturbed gene expression data and corresponding metadata
#'
#' @description
#' `test_data()` validates perturbed gene expression data and its corresponding metadata prior to
#' downstream analysis. It ensures that the expression matrix and metadata table are compatible,
#' correctly formatted, and internally consistent with respect to sample identifiers.
#'
#' @param ge_matrix A matrix or array of perturbed gene expression values with `probes` in rows and
#'   `barcode` in columns.
#' @param metadata A data frame or tibble containing sample metadata. Must include the columns
#'   `barcode`, `compound_name`, `dose_level`, and `time_level`.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It returns `NULL` invisibly if validation
#' passes; otherwise it signals an error.
#'
#' @details
#' The function checks that:
#' \itemize{
#'   \item `ge_matrix` is a matrix or array.
#'   \item `metadata` is a data frame (including tibbles).
#'   \item Required metadata columns are present.
#'   \item All expression columns are described by `metadata$barcode`.
#' }
#'
#' @examples
#' tgx_data <- simulate_tgxdata(n_de = 20, n_ee = 10, n_com = c(3, 3))
#' test_data(tgx_data$expression, tgx_data$metadata)
#'
#' \dontrun{
#' # Invalid class of metadata
#' gemat_invalid <- array(rnorm(100), dim = c(10, 10, 1))
#' metadata_invalid <- list(barcode = paste0("sample", 1:10))
#' test_data(gemat_invalid, metadata_invalid)
#'
#' # Invalid barcode
#' tgx_data <- simulate_tgxdata(n_de = 20, n_ee = 10, n_com = c(3, 3))
#' tgx_data$metadata$barcode <- paste0("sample", 1:nrow(tgx_data$metadata))
#' test_data(tgx_data$expression, tgx_data$metadata)
#' }
#'
#' @export
test_data <- function(ge_matrix, metadata, error_call = rlang::caller_env()) {
  if (!inherits(ge_matrix, c("matrix", "array"))) {
    cli_abort(c("The expression data must be a matrix.",
                "x" = "The class {style_bold(col_cyan(backtick(class(ge_matrix))))} of \\
                expression data {?is/are} not supported.",
                "i" = "Please make expression data as matrix \\
             (`probes` in rows and `barcode` in columns)."), call = error_call)

  }
  if (!rlang::inherits_any(metadata, c("tbl_df", "tbl", "data.frame"))) {
    cli::cli_abort(c("The metadata must be a data frame.",
                "x" = "The class {style_bold(col_cyan(backtick(class(metadata))))} of \\
                attribute data {?is/are} not supported.",
                "i" = "Please make meatadata as data frame."), call = error_call)

  }
  test_column(c("barcode", "compound_name", "dose_level", "time_level"), metadata)
  match_bar <- colnames(ge_matrix) %in% metadata$barcode
  if (any(!match_bar)) {
    cli::cli_abort(c("All columns/samples in expression data must given in metadata.",
                "x" = "{style_bold(col_cyan(backtick(colnames(ge_matrix)[match_bar])))} column{?s} \\
                {?is/are} not described in the metadata.",
                "i" = "Please remove unspecified columns from the ge_matrix. After that you \\
               can update your data by calling the function `update_data()`"), call = error_call)
  }
}


#' Validate and standardize compound groups
#'
#' @description
#' `test_group()` validates compound groups supplied through `...` and returns them in a consistent
#' format. Inputs may be provided as multiple vectors/factors/lists, or as a single list containing
#' groups. The function checks that the supplied groups are of supported types and provides
#' informative errors when inputs are malformed or inconsistent.
#'
#' @param ... Compound groups supplied as vectors, factors, or lists containing compound names.
#'   If a single list is provided, it is treated as the set of groups; otherwise,
#'   each argument in `...` is treated as a separate group.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return A named list of compound groups. If group names are not provided, they are assigned as
#' `"Group_A"`, `"Group_B"`, etc.
#'
#' @examples
#' # Two groups provided as separate vectors
#' group1 <- c("compound1", "compound2")
#' group2 <- c("compound3", "compound4")
#' test_group(group1, group2)
#'
#' # Groups provided as a single list
#' groups <- list(group1 = c("compound1", "compound2"),
#'                group2 = c("compound3", "compound4"))
#' test_group(groups)
#'
#' \dontrun{
#' # Mixed input types may trigger an error depending on the supplied objects
#' test_group(group1, group2, list("compound5"))
#' }
#'
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
  if (!all(arg_checker )) {
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

#' Validate data structure for the `toxassay` class
#'
#' @description
#' `test_datastr()` validates that an input object conforms to the expected structure of a
#' `ToxAssay`. The input must be either an object of class `toxassay` or a list of length five,
#' representing the experimental hierarchy: group, compound, dose, time, and replication.
#'
#' @param expr_str An object to be validated. Must be either a `toxassay` object or a list with
#'   exactly five elements corresponding to `group`, `compound`, `dose`, `time`, and `replication`.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @return
#' This function is called for its side effects only. It returns `TRUE` invisibly if validation
#' succeeds; otherwise it signals an error.
#'
#' @examples
#' # Example of a valid list structure
#' valid_list <- list(
#'   group = "A",
#'   compound = "Methimazole",
#'   dose = "High",
#'   time = "24h",
#'   replication = 3)
#' test_datastr(valid_list)
#'
#' \dontrun{
#' # Invalid list (incorrect length)
#' invalid_list <- list(group = "A", compound = "Methimazole")
#' test_datastr(invalid_list)
#' }
#'
#' @export
test_datastr <- function(expr_str, error_call = rlang::caller_env()) {
  if (!is.list(expr_str) & !inherits(expr_str, "toxassay")) {
    cli::cli_abort(c("{.var expr_str} must be object of class `toxassay` or `list`.",
                "x" = "The class {style_bold(col_cyan(backtick(class(expr_str))))} of \\
                {.var expr_str} is not supported.",
                "i" = "Please provide {.var expr_str} as a list of size 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."), call = error_call)
  } else if (length(expr_str)  !=  5) {
    cli::cli_abort(c("The length of {.var expr_str} must be 5.",
                "i" = "You have supplied a list {.var expr_str} of size {length(expr_str)}, please
                make sure {.var expr_str} has a length of 5
                (for `group`, `compound`, `dose`, `time` and `replication`)."), call = error_call)
  }
}

#' Validate metadata for toxicogenomics datasets
#'
#' @description
#' `is_metadata()` checks that metadata used for toxicogenomics analyses are complete and internally
#' consistent. The function validates required columns and enforces that key study descriptors are
#' unique within a dataset (e.g., a single `database` and a single data type in `fc`). Additional
#' dataset-specific checks are applied for `tggates` and `drugmatrix`.
#'
#' @param metadata A data frame containing sample metadata. Must include the columns
#'   `compound_name`, `dose_level`, `time_level`, `organ_id`, `database`, `arr_design`, and `fc`.
#'   For `database = "tggates"`, the columns `species`, `test_type`, and `sin_rep_type` are also required.
#' @param error_call The environment used for error reporting. Default is `rlang::caller_env()`.
#'
#' @details
#' The function performs the following checks:
#' \itemize{
#'   \item Required columns are present in `metadata`.
#'   \item `database` is valid and unique (one of `"tggates"` or `"drugmatrix"`).
#'   \item `fc` is unique (e.g., only one of `"FC"` or `"normal"`).
#'   \item For `database = "tggates"`: `species`, `test_type`, `organ_id`, and `sin_rep_type` are each unique.
#'   \item For `database = "drugmatrix"`: `organ_id` is unique.
#' }
#'
#' @return
#' This function is called for its side effects only. It returns `NULL` invisibly if validation
#' succeeds; otherwise it signals an error.
#'
#' @examples
#' metadata <- data.frame(
#'   compound_name = c("Compound A", "Compound B"),
#'   dose_level = c("High", "Low"),
#'   time_level = c("24h", "48h"),
#'   organ_id = "Liver",
#'   database = "drugmatrix",
#'   arr_design = c("Design1", "Design2"),
#'   fc = "FC")
#' test_metadata(metadata)
#' \dontrun{
#' # Invalid: multiple databases present
#' metadata_bad_db <- data.frame(
#' compound_name = c("Compound A", "Compound B"),
#' dose_level = c("High", "Low"),
#' time_level = c("24h", "48h"),
#' organ_id = c("Liver", "Kidney"),
#' database = c("drugmatrix", "tggates"),
#' arr_design = c("Design1", "Design2"),
#' fc = "FC")
#' test_metadata(metadata_bad_db)
#'
#' # Invalid: mixed FC and normal data
#' metadata_bad_fc <- data.frame(
#' compound_name = c("Compound A", "Compound B"),
#' dose_level = c("High", "Low"),
#' time_level = c("24h", "48h"),
#' organ_id = "Liver",
#' database = "drugmatrix",
#' arr_design = c("Design1", "Design2"),
#' fc = c("FC", "normal"))
#' test_metadata(metadata_bad_fc)
#' }
#'
#' @export
test_metadata <- function(metadata, error_call = rlang::caller_env()) {
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
