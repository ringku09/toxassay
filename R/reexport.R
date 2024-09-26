#-------------------------------  Base packages ---------------------
#' @importFrom stats na.omit
#' @export
stats::na.omit

#' @importFrom stats confint
#' @export
stats::confint

#' @importFrom stats pf
#' @export
stats::pf

#' @importFrom methods as
#' @export
methods::as

#' @importFrom scales rescale
#' @export
scales::rescale

#' @importFrom utils unzip
#' @export
utils::unzip

#' @importFrom utils unzip
#' @export
utils::unzip

#' @importFrom utils read.table
#' @export
utils::read.table

#' @importFrom utils txtProgressBar
#' @export
utils::txtProgressBar

#' @importFrom utils setTxtProgressBar
#' @export
utils::setTxtProgressBar

#' @importFrom utils download.file
#' @export
utils::download.file

#' @importFrom utils unzip
#' @export
utils::unzip

#' @importFrom utils packageVersion
#' @export
utils::packageVersion

#' @importFrom tools file_ext
#' @export
tools::file_ext

#' @importFrom stats median
#' @export
stats::median

#-------------------------------  Biological database ---------------------
  #' @importFrom AnnotationDbi select
  #' @export
  AnnotationDbi::select

  #' @importFrom rat2302.db rat2302.db
  #' @export
  rat2302.db::rat2302.db


#' @importFrom magrittr %>%
#' @export
 magrittr::`%>%`

 #' @importFrom affy ReadAffy
 #' @export
 affy::ReadAffy

 #' @importFrom affy mas5
 #' @export
 affy::mas5

 #' @importFrom affy rma
 #' @export
 affy::rma

 #' @importFrom Biobase exprs
 #' @export
 Biobase::exprs

 #----------------------------- import for `dplyr`--------------------
 #' @importFrom dplyr across
 #' @export
 dplyr::across

 #' @importFrom dplyr pull
 #' @export
 dplyr::pull

 #' @importFrom dplyr distinct
 #' @export
 dplyr::distinct

 #' @importFrom dplyr ungroup
 #' @export
 dplyr::ungroup

 #' @importFrom dplyr arrange
 #' @export
 dplyr::arrange

 #' @importFrom dplyr group_by
 #' @export
 dplyr::group_by

 #' @importFrom dplyr select
 #' @export
 dplyr::select

 #' @importFrom dplyr select
 #' @export
 dplyr::arrange

 #' @importFrom dplyr filter
 #' @export
 dplyr::filter

 #' @importFrom dplyr rename
 #' @export
 dplyr::rename

 #' @importFrom dplyr rename_with
 #' @export
 dplyr::rename_with

 #' @importFrom dplyr left_join
 #' @export
 dplyr::left_join

 #' @importFrom  dplyr bind_rows
 #' @export
 dplyr::bind_rows



 #------------------------------------- Tidyselect ---------------

 #' @importFrom tidyselect everything
 #' @export
 tidyselect::everything

 #' @importFrom dplyr mutate_at
 #' @export
 dplyr::mutate_at

 #' @importFrom tibble as_tibble
 #' @export
 tibble::as_tibble

 #' @importFrom readr write_csv
 #' @export
 readr::write_csv

 #------------------------------------ parralel computing --------------------
 #' @importFrom parallel detectCores
 #' @export
 parallel::detectCores

 #' @importFrom parallel splitIndices
 #' @export
 parallel::splitIndices

 #' @importFrom parallel makeCluster
 #' @export
 parallel::makeCluster

 #' @importFrom parallel clusterExport
 #' @export
 parallel::clusterExport

 #' @importFrom parallel clusterCall
 #' @export
 parallel::clusterCall

 #' @importFrom parallel stopCluster
 #' @export
 parallel::stopCluster

 #' @importFrom doParallel registerDoParallel
 #' @export
 doParallel::registerDoParallel

 #' @importFrom foreach getDoParName
 #' @export
 foreach::getDoParName

 #' @importFrom foreach getDoParWorkers
 #' @export
 foreach::getDoParWorkers

 #' @importFrom foreach registerDoSEQ
 #' @export
 foreach::registerDoSEQ

 #' @importFrom foreach foreach
 #' @export
 foreach::foreach

 #' @importFrom foreach `%dopar%`
 #' @export
 foreach::`%dopar%`

 #' @importFrom foreach `%do%`
 #' @export
 foreach::`%do%`



 #----------------------------- import for `rlang`--------------------
 #' @importFrom rlang is_empty
 #' @export
 rlang::is_empty

 #' @importFrom rlang .data
 #' @export
 rlang::.data

 #' @importFrom rlang englue
 #' @export
 rlang::englue

 #' @importFrom rlang list2
 #' @export
 rlang::list2

 #' @importFrom rlang is_bare_list
 #' @export
 rlang::is_bare_list

 #' @importFrom rlang inherits_any
 #' @export
 rlang::inherits_any

 #' @importFrom rlang is_missing
 #' @export
 rlang::is_missing

 #' @importFrom rlang missing_arg
 #' @export
 rlang::missing_arg

 #' @importFrom rlang parse_expr
 #' @export
 rlang::parse_expr
 #----------------------------- import for `cli`--------------------
 #' @importFrom cli cli_alert_warning
 #' @export
 cli::cli_alert_warning

 cli_abort
 #' @importFrom cli cli_abort
 #' @export
 cli::cli_abort

 #' @importFrom cli style_bold
 #' @export
 cli::style_bold

 #' @importFrom cli cli_ul
 #' @export
 cli::cli_ul

 #' @importFrom cli bg_magenta
 #' @export
 cli::bg_magenta

 #' @importFrom cli bg_red
 #' @export
 cli::bg_red

 #' @importFrom cli bg_red
 #' @export
 cli::bg_red

 #' @importFrom cli col_red
 #' @export
 cli::col_red

 #' @importFrom cli col_green
 #' @export
 cli::col_green

 #' @importFrom cli cli_alert_info
 #' @export
 cli::cli_alert_info

 #' @importFrom cli col_br_red
 #' @export
 cli::col_br_red

 #' @importFrom cli cli_alert_success
 #' @export
 cli::cli_alert_success

 #----------------------------- import for `glue`--------------------
 #' @importFrom glue backtick
 #' @export
 glue::backtick

 #' @importFrom glue glue
 #' @export
 glue::glue

 #' @importFrom glue glue_collapse
 #' @export
 glue::glue_collapse

 #----------------------------- import for `tidyr`--------------------
 #' @importFrom tidyr nest
 #' @export
 tidyr::nest

 #' @importFrom tidyr pivot_wider
 #' @export
 tidyr::pivot_wider

 #----------------------------- import for `purrr`--------------------

 #' @importFrom  purrr map
 #' @export
 purrr::map

 #' @importFrom  purrr imap
 #' @export
 purrr::imap

 #----------------------------- import for `readr`--------------------

 #' @importFrom  readr write_tsv
 #' @export
 readr::write_tsv

 #' @importFrom  readr write_csv
 #' @export
 readr::write_csv

 #----------------------------- import for `tibble`--------------------
 #' @importFrom tibble tibble
 #' @export
 tibble::tibble

 #' @importFrom tibble as_tibble
 #' @export
 tibble::as_tibble

 #---------------------------------- httr -----------------------

 #' @importFrom httr GET
 #' @export
 httr::GET

 #' @importFrom httr write_disk
 #' @export
 httr::write_disk

 #' @importFrom httr progress
 #' @export
 httr::progress

 #' @importFrom  httr write_disk
 #' @export
 httr::write_disk


 #############################################

#
#  #' @import arules
#  #' @export

#' @import STRINGdb
#' @export


