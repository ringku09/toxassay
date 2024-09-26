
#' Structure of levels in the data for a compound
#'
#' @description
#'
#' The `compound_str()` function is used to find the total number of dose levels,
#' the time points within each dose, and the number of replicates used for each
#' dose-time combination.
#'
#' The `compound_str()` function calculates the total number of dose levels used for a
#' compound in the experiment, as well as the total number of time points for each dose
#' and the total number of replicates for each time point within each dose. For instance,
#' if a compound has 3 dose labels, there are 3 labels for the dose. If each dose label
#' is measured at 4 time points, then the number of labels for time levels would be
#' (4, 4, 4) corresponding to the 3 dose labels. Finally, if each dose-time combination
#' uses 3 replicates, then the number of labels for replicates would be 3 * 4 = 12
#' (3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3). Therefore, the total number of labels at each
#' level is a product of the labels from earlier levels.
#'
#' @param compound Compound name (or abbreviation for `TG-GATEs` database)
#' @inheritParams expand_data
#'
#' @returns
#' A list of dose, time and replication
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' compound_str(compound = "Compound1",  metadata = sim_data$metadata)
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

#' Structure of experimental data for analyzing targeted toxicity
#'
#' The `data_str()` function finds the number of levels in dose, time point, and
#' replication for a list of compounds in different groups. `data_str()` uses the
#' function `data_str()` to determine the structure of levels in compounds and
#' concatenates them to create the data structure for multiple compounds in a group.
#'
#' @param ... Either multiple vectors of individual compound names or multiple vectors of
#'   compound names grouped in a list.
#' @inheritParams expand_data
#'
#' @seealso
#'    [compound_str()] for getting structure of a compound
#'
#' @return
#' A list total levels in different stage or a class of `ToxAssay` object .
#'  * group
#'  * compound
#'  * dose
#'  * time
#'  * replication
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' data_str(gr, metadata = sim_data$metadata)
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
  structure(comps_lev, class = "ToxAssay")
  #return(comps_lev)
}

#' Quadratic matrices for F-statistics
#'
#' The `get_quadmat()` function calculates the quadratic matrices **A** and **B**
#'
#' @param expr_str The structure of the experiment, either a list detailing the total
#'   number of labels at each level or a `ToxAssay` object representing the data
#'   structure, can be obtained using the [data_str()] function.
#'
#' @return a list of design matrices and others information in the experiment
#'   * `Fstat_A`, quadratic matrix A
#'   * `BFstat_A`, quadratic matrix B
#'   * `group_mat`, summary matrix used to calculate average expression at group level
#'   * `compound_mat`, summary matrix used to calculate average expression at compound level
#'   * `dose_mat`, summary matrix used to calculate average expression at dose level
#'   * `time_mat`, summary matrix used to calculate average expression at time level
#'
#' @export
#'
#' @examples
#' sim_data <- simulate_tgxdata(n_de = 10, n_ee = 10, n_com = c(5,5))
#' gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
#' sim_str <- data_str(gr, metadata = sim_data$metadata)
#' summary_mat <- get_matrix(sim_str)
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
  structure(retn, class = "ToxAssay")
}

