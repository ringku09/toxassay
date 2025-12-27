#' Matrix class imports
#'
#' @importClassesFrom Matrix dgCMatrix
NULL

#' Re-export selected functions from stats
#'
#' @importFrom stats na.omit confint pf median p.adjust
#' @export
stats::na.omit
#' @export
stats::confint
#' @export
stats::pf
#' @export
stats::median
#' @export
stats::p.adjust

#' Re-export selected functions from methods
#'
#' @importFrom methods as
#' @export
methods::as

#' Re-export selected functions from scales
#'
#' @importFrom scales rescale
#' @export
scales::rescale

#' Re-export selected functions from utils
#'
#' @importFrom utils unzip read.table txtProgressBar setTxtProgressBar download.file packageVersion
#' @export
utils::unzip
#' @export
utils::read.table
#' @export
utils::txtProgressBar
#' @export
utils::setTxtProgressBar
#' @export
utils::download.file
#' @export
utils::packageVersion

#' Re-export selected functions from tools
#'
#' @importFrom tools file_ext
#' @export
tools::file_ext

#--------------------------- Biological databases -------------------

#' Re-export selected functions/objects from Bioconductor annotation packages
#'
#' @importFrom AnnotationDbi select
#' @importFrom rat2302.db rat2302.db
#' @export
AnnotationDbi::select
#' @export
rat2302.db::rat2302.db

#------------------------------- magrittr ---------------------------

#' Re-export the pipe operator
#'
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#------------------------------- affy / Biobase ---------------------

#' Re-export selected functions from affy
#'
#' @importFrom affy ReadAffy mas5 rma
#' @export
affy::ReadAffy
#' @export
affy::mas5
#' @export
affy::rma

#' Re-export selected functions from Biobase
#'
#' @importFrom Biobase exprs
#' @export
Biobase::exprs

#------------------------------- dplyr ------------------------------

#' Re-export selected functions from dplyr
#'
#' @importFrom dplyr across pull distinct ungroup arrange group_by select filter rename rename_with left_join bind_rows mutate_at
#' @export
dplyr::across
#' @export
dplyr::pull
#' @export
dplyr::distinct
#' @export
dplyr::ungroup
#' @export
dplyr::arrange
#' @export
dplyr::group_by
#' @export
dplyr::select
#' @export
dplyr::filter
#' @export
dplyr::rename
#' @export
dplyr::rename_with
#' @export
dplyr::left_join
#' @export
dplyr::bind_rows
#' @export
dplyr::mutate_at

#------------------------------ tidyselect --------------------------

#' Re-export selected functions from tidyselect
#'
#' @importFrom tidyselect everything
#' @export
tidyselect::everything

#------------------------------- tidyr ------------------------------

#' Re-export selected functions from tidyr
#'
#' @importFrom tidyr nest pivot_wider
#' @export
tidyr::nest
#' @export
tidyr::pivot_wider

#------------------------------- purrr ------------------------------

#' Re-export selected functions from purrr
#'
#' @importFrom purrr map imap
#' @export
purrr::map
#' @export
purrr::imap

#------------------------------- readr ------------------------------

#' Re-export selected functions from readr
#'
#' @importFrom readr write_csv write_tsv
#' @export
readr::write_csv
#' @export
readr::write_tsv

#------------------------------- tibble -----------------------------

#' Re-export selected functions from tibble
#'
#' @importFrom tibble tibble as_tibble
#' @export
tibble::tibble
#' @export
tibble::as_tibble

#------------------------------ parallel ----------------------------

#' Re-export selected functions from parallel
#'
#' @importFrom parallel detectCores splitIndices makeCluster clusterExport clusterCall stopCluster
#' @export
parallel::detectCores
#' @export
parallel::splitIndices
#' @export
parallel::makeCluster
#' @export
parallel::clusterExport
#' @export
parallel::clusterCall
#' @export
parallel::stopCluster

#----------------------------- doParallel ---------------------------

#' Re-export selected functions from doParallel
#'
#' @importFrom doParallel registerDoParallel
#' @export
doParallel::registerDoParallel

#------------------------------- foreach ----------------------------

#' Re-export selected functions from foreach
#'
#' @importFrom foreach getDoParName getDoParWorkers registerDoSEQ foreach %dopar% %do%
#' @export
foreach::getDoParName
#' @export
foreach::getDoParWorkers
#' @export
foreach::registerDoSEQ
#' @export
foreach::foreach
#' @export
foreach::`%dopar%`
#' @export
foreach::`%do%`

#------------------------------- rlang ------------------------------

#' Re-export selected functions from rlang
#'
#' @importFrom rlang is_empty is_missing caller_env .data englue list2 is_bare_list inherits_any missing_arg parse_expr
#' @export
rlang::is_empty
#' @export
rlang::is_missing
#' @export
rlang::caller_env
#' @export
rlang::.data
#' @export
rlang::englue
#' @export
rlang::list2
#' @export
rlang::is_bare_list
#' @export
rlang::inherits_any
#' @export
rlang::missing_arg
#' @export
rlang::parse_expr

#-------------------------------- cli -------------------------------

#' Re-export selected functions from cli
#'
#' @importFrom cli cli_abort cli_alert_warning cli_alert_info cli_alert_success style_bold cli_ul
#' @importFrom cli bg_magenta bg_red col_red col_green col_br_red
#' @export
cli::cli_abort
#' @export
cli::cli_alert_warning
#' @export
cli::cli_alert_info
#' @export
cli::cli_alert_success
#' @export
cli::style_bold
#' @export
cli::cli_ul
#' @export
cli::bg_magenta
#' @export
cli::bg_red
#' @export
cli::col_red
#' @export
cli::col_green
#' @export
cli::col_br_red

#-------------------------------- glue ------------------------------

#' Re-export selected functions from glue
#'
#' @importFrom glue backtick glue glue_collapse
#' @export
glue::backtick
#' @export
glue::glue
#' @export
glue::glue_collapse

#-------------------------------- httr ------------------------------

#' Re-export selected functions from httr
#'
#' @importFrom httr GET write_disk progress
#' @export
httr::GET
#' @export
httr::write_disk
#' @export
httr::progress

#-------------------------------- arules ------------------------------

#' Re-export selected functions from arules
#'
#' @importFrom arules intersect apriori subset %in% quality interestMeasure
#' @export
arules::intersect
#' @export
arules::apriori
#' @export
arules::subset
#' @export
arules::`%in%`
#' @export
arules::quality
#' @export
arules::interestMeasure

#------------------------------ Package imports ---------------------

#' Import STRINGdb
#'
#' @import STRINGdb
NULL
