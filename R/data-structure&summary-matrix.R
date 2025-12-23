#' Summarize dose, time, and replicate structure for a compound
#'
#' @description
#' `compound_str()` summarizes the experimental design for a single compound by reporting (1) the
#' number of dose levels, (2) the number of time points within each dose level, and (3) the number of
#' samples (replicates) for each dose–time combination.
#'
#' The returned values describe the hierarchical structure of the data. For example, if a compound has
#' 3 dose levels and each dose is measured at 4 time points, then `dose = 3` and `time = c(4, 4, 4)`.
#' If each dose–time combination has 3 replicates, then `replication` will contain 12 entries (3 doses ×
#' 4 times), each equal to 3.
#'
#' @param compound Name of the compound.
#' @inheritParams expand_data
#'
#' @return A list with three elements:
#' \item{dose}{Number of dose levels for `compound`.}
#' \item{time}{Integer vector giving the number of time points within each dose level.}
#' \item{replication}{Integer vector giving the number of replicates for each dose–time combination.}
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' compound_str(compound = "Compound1", metadata = sim_data$metadata)
#'
#' @export
compound_str <- function(compound, metadata) {
  if (length(compound) > 1) {
    cli_abort(c(
      "Multiple compounds are not allowed",
      "x" = "You have provided compounds {style_bold(col_red(backtick(compound)))}.",
      "i" = "Please choose either {add_or(style_bold(col_green(backtick(compound))))} instead."
    ))
  }
  test_column(c("compound_name", "dose_level", "time_level"), metadata)
  if (! compound %in% metadata$compound_name) {
    cli_abort(c(
      "Sample information of quired compound must be available in {.arg metadata}.",
      "x" = "Sample information of {style_bold(col_red(backtick(compound)))} are not available in {.arg metadata}.",
      "i" = "Please choose appropiate compound or provide {.arg metadata} correctly."
    ))
  }
  attr_data2 <- metadata %>% dplyr::filter(.data$compound_name == compound)
  dose_df <- attr_data2 %>%
    dplyr::group_by(.data$dose_level) %>%
    tidyr::nest()
  time_df <-lapply(dose_df$data, function(x) x %>%
                     dplyr::group_by(.data$time_level) %>%
                     tidyr::nest())
  dose_lev <- nrow(dose_df)
  time_lev <- sapply(time_df, nrow)
  rep_lev <- unlist(lapply(time_df, function(x) sapply(x$data, nrow)))
  lev_list <- list(dose_lev, time_lev, rep_lev)
  names(lev_list) <- c("dose", "time", "replication")
  return(lev_list)
}

#' Summarize experimental data structure for compound groups
#'
#' @description
#' `data_str()` summarizes the hierarchical structure of an experimental design involving multiple
#' compounds organized into one or more groups. For each compound, the function determines the
#' number of dose levels, time points, and replicates using `compound_str()`, and then aggregates
#' this information across compounds and groups to describe the overall data structure.
#'
#' @param ... One or more character vectors of compound names, or a list of character vectors
#'   defining compound groups.
#' @inheritParams expand_data
#'
#' @return
#' An object of class `toxassay`, represented as a list with the following elements:
#' \item{group}{Number of compound groups.}
#' \item{compound}{Integer vector giving the number of compounds in each group.}
#' \item{dose}{Integer vector giving the number of dose levels for each compound.}
#' \item{time}{Integer vector giving the number of time points for each dose level.}
#' \item{replication}{Integer vector giving the number of replicates for each dose–time combination.}
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' data_str(gr, metadata = sim_data$metadata)
#'
#' @seealso
#' [compound_str()] for summarizing the structure of a single compound.
#'
#' @export
data_str <- function(..., metadata, error_call = caller_env()) {
  comps_gr <- test_group(...)
  group <- length(comps_gr)
  comp_gr <- unlist(lapply(comps_gr, length))
  comps_list <- as.list(unlist(comps_gr))
  zz <- lapply(comps_list, compound_str, metadata = metadata)
  dose_lev <- do.call(c,lapply(zz, "[[", 1))
  time_lev <- do.call(c,lapply(zz, "[[", 2))
  rep_lev <- do.call(c,lapply(zz, "[[", 3))
  comps_lev <- list(group, comp_gr, dose_lev, time_lev, rep_lev)
  names(comps_lev) <- c("group", "compound", "dose", "time", "replication")
  structure(comps_lev, class = "toxassay")
  #return(comps_lev)
}

#' Quadratic and summary matrices for F-statistics of HLM
#'
#' @description
#' `get_matrix()` constructs quadratic matrices used for F-statistics and summary (averaging) matrices
#' for each level of the experimental hierarchy (group, compound, dose, and time). The input
#' `expr_str` describes the nesting structure of the experiment and can be created with [data_str()].
#'
#' @param expr_str The structure of the experiment, either a list describing the total number of
#'   labels at each level or a `ToxAssay` object representing the data structure (e.g., returned by
#'   [data_str()]).
#'
#' @return
#' An object of class `toxasaay` containing:
#' \item{A}{Quadratic matrix used in the F-statistic calculation.}
#' \item{B}{Quadratic matrix used in the F-statistic calculation.}
#' \item{group_mat}{Summary matrix used to compute average expression at the group level.}
#' \item{compound_mat}{Summary matrix used to compute average expression at the compound level.}
#' \item{dose_mat}{Summary matrix used to compute average expression at the dose level.}
#' \item{time_mat}{Summary matrix used to compute average expression at the time level.}
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5, 5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sim_str <- data_str(gr, metadata = sim_data$metadata)
#' summary_mat <- get_matrix(sim_str)
#'
#' @export
get_matrix <- function(expr_str) {
  test_datastr(expr_str)
  a <- expr_str[[1]]
  b_i <- expr_str[[2]]
  c_ij <- expr_str[[3]]
  d_ijk <- expr_str[[4]]
  n_ijkl <- expr_str[[5]]
  N <-  sum(expr_str[[5]])
  n_ijk. = n_ij.. = n_i... <- vector()
  Qa = Qd = qAlfa1 = qBeta1 = qGama1 = qDelta1 <-  matrix(nrow=0, ncol=0)
  for (i in 1 : a) {
    ni... <- 0
    for (j in 1 : b_i[i]) {
      nij.. <- 0
      for (k in 1 : c_ij[j]) {
        nijk. <- 0
        for (l in 1 : d_ijk[k]) {
          nijk. <- nijk. + n_ijkl[l]
          qDelta1 <- direct_sum(qDelta1, matrix(rep(1, n_ijkl[l]), ncol = 1) / n_ijkl[l])
          Qd <- direct_sum(Qd, mat_one(n_ijkl[l], n_ijkl[l]) / n_ijkl[l])
        }
        n_ijkl <- n_ijkl[-(1 : l)]
        n_ijk. <- c(n_ijk., nijk.)
        qGama1 <- direct_sum(qGama1, matrix(rep(1,  nijk.), ncol = 1) / nijk.)
        nij.. <- nij.. + nijk.
      }
      n_ij.. <- c(n_ij.., nij..)
      qBeta1 <- direct_sum(qBeta1, matrix(rep(1, nij..), ncol = 1) / nij..)
      d_ijk <- d_ijk[-(1 : k)]
      ni... <- ni... + nij..
    }
    n_i... <- c(n_i..., ni...)
    qAlfa1 <- direct_sum(qAlfa1, matrix(rep(1, ni...), ncol = 1) / ni...)
    Qa <- direct_sum(Qa, mat_one(ni..., ni...) / ni...)
    c_ij  <- c_ij[-(1 : j)]
  }
  QA <- Qa - mat_one(N, N) / N
  R <- diag(N) - Qd
  retn <- list(
    A = QA,
    B = R,
    group_mat = qAlfa1,
    compound_mat = qBeta1,
    dose_mat = qGama1,
    time_mat = qDelta1
  )
  structure(retn, class = "toxassay")
}

