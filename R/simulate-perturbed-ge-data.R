#' Simulate variance components from a target ICC
#'
#' @description
#' `sim_variance()` generates a set of variance components whose total proportion corresponds to a
#' specified intracluster correlation coefficient (ICC). The components can be generated with equal
#' expected contribution or with unequal contributions across components.
#'
#' @param icc A numeric value specifying the intracluster correlation coefficient (ICC).
#' @param n_component An integer specifying the number of variance components to simulate.
#' @param equal_var A logical value indicating whether variance components are generated with equal
#'   expected weights (`TRUE`) or unequal weights (`FALSE`). Default is `TRUE`.
#'
#' @return A numeric vector of length `n_component` containing simulated variance components whose
#' sum equals `icc`.
#'
#' @examples
#' # Simulate three variance components with equal expected weights
#' sim_variance(icc = 0.5, n_component = 3, equal_var = TRUE)
#'
#' # Simulate three variance components with unequal weights
#' sim_variance(icc = 0.5, n_component = 3, equal_var = FALSE)
#'
#' @export
sim_variance <- function(icc, n_component, equal_var = TRUE) {
  if (equal_var) {
    x <- rdirichlet(1, rep(100, n_component))
    x_norm <- x * icc
  } else {
    x <- runif(n_component, 0, 1)
    x_norm <- (x / sum(x)) * icc
  }
  return(x_norm)
}

#' Simulate perturbed gene expression data
#'
#' @description
#' `simulate_ge()` generates simulated expression values for a single gene under a multi-factor
#' experimental design defined by compounds, dose levels, time points, and replicates. The simulation
#' includes a group effect (controlled by `d`), correlation within experimental units (controlled by
#' `icc`), and residual (error) variation such that the total variance is standardized to 1 when the
#' requested settings are feasible.
#'
#' @param n_com Numeric vector indicating the number of compounds in each group (e.g., toxic and
#'   non-toxic). Typically length 2.
#' @param n_dose Integer specifying the number of dose levels. Default is `3`.
#' @param n_time Integer specifying the number of time points. Default is `4`.
#' @param n_rep Integer specifying the number of replicates for each dose–time condition. Default is `3`.
#' @param mu Numeric value specifying the overall mean expression level. Default is a random value in
#'   `[-3, 3]`.
#' @param d Numeric value specifying the effect size (Cohen's d) between the two groups. Default is `0.5`.
#' @param icc Numeric value specifying the intracluster correlation coefficient. Values must be in
#'   `[0, 1)`. Default is `0.5`.
#' @param error_call Environment used for error reporting. Default is `caller_env()`.
#'
#' @details
#' The function first determines the explained variation attributable to the requested `icc` and group
#' effect implied by `d`. If the resulting residual variance would be negative, an error is raised.
#' Variance components are then simulated (compound, dose, and time), and expression values are drawn
#' from a normal distribution with mean determined by `mu`, the group effect, and the simulated
#' hierarchical effects, and variance equal to the residual variance.
#'
#' @return A numeric vector of length `sum(n_com) * n_dose * n_time * n_rep` containing simulated gene
#' expression values.
#'
#' @examples
#' # Simulate gene expression data
#' simulated_data <- simulate_ge(n_com = c(3, 3))
#'
#' # Simulate with custom design parameters
#' simulated_data <- simulate_ge(n_com = c(5, 5), n_dose = 2, n_time = 3, n_rep = 2, d = 0.8, icc = 0.3)
#'
#' @export
simulate_ge <- function(n_com,
                        n_dose = 3,
                        n_time = 4,
                        n_rep = 3,
                        mu = runif(1, -3, 3),
                        d = 0.5,
                        icc = 0.5,
                        error_call = caller_env()) {
  icc <- abs(icc)
  r2 <- (d/sqrt(d^2+4))^2
  exp_var <- icc + r2
  v_err <- 1 - exp_var
  # if (length(n_com) > 2) {
  #   cli::cli_abort(c("Simulation only support two groups of compound.",
  #               "x" = "You have supplied {style_bold(col_red(n_com))} compound groups.",
  #               "i" = "Please provide number of compound groups of length 2."),
  #             wrap = TRUE, call = error_call)
  # }
  if(icc >= 1){
    cli::cli_abort(c("The data must be standerdize.",
                "x" = "You have supplied overall ICC of {style_bold(col_red(icc))}.",
                "i" = "Please provide ovarall ICC in between 0 and 1."),
              wrap = TRUE, call = error_call)
  }

  if(v_err < 0){
    cli::cli_abort(c("The data must be standerdize.",
                "x" = "You have supplied {.arg icc} and {.arg d} which result in sum of \\
                explained (icc) and unexplained (r2) variation is {style_bold(col_red(round(exp_var,2)))}.",
                "i" = "Please provide reduce values of {.arg icc} and {.arg d}."),
              wrap = TRUE, call = error_call)
  }
  n_coms <- sum(n_com)
  nsample <- n_coms * n_dose * n_time * n_rep
  s <- sim_variance(icc = icc, n_component = 3)
  s_com <- sqrt(s[1])
  s_dose <- sqrt(s[2])
  s_time <- sqrt(s[3])
  m_com <- rnorm(n_coms, 0, s_com)
  m_dose <- lapply(as.list(unlist(m_com)), function(x) rnorm(n_dose,x, s_dose))
  m_time <- lapply(as.list(unlist(m_dose)), function(x) rnorm(n_time,x,s_time))
 # y_temp1 <- unlist(lapply(as.list(unlist(m_time)), function(x) rnorm(n_rep,x, sqrt(v_err))))
  y_temp <- mu + 0.5*rep(c(-d,d), each = nsample/2) + rep(unlist(m_time), n_rep)
  y <- rnorm(nsample, y_temp , sqrt(v_err))
  return(y)
}

#' Simulate perturbed gene expression matrix
#'
#' @description
#' `simulate_data()` generates a simulated gene expression matrix for a multi-factor experimental
#' design with groups of compounds, dose levels, time points, and replicates. Expression values for
#' each gene are generated by repeatedly calling `simulate_ge()`, and sample-level metadata are
#' created to describe the simulated design.
#'
#' @param n_gene Integer specifying the number of genes to simulate.
#' @param n_com Integer vector specifying the number of compounds in each group.
#' @param n_dose Integer specifying the number of dose levels. Default is `3`.
#' @param n_time Integer specifying the number of time points. Default is `4`.
#' @param n_rep Integer specifying the number of replicates for each dose–time condition. Default is `3`.
#' @param mu Numeric value specifying the mean expression level used in simulation. Default is a random
#'   value in `[-3, 3]`.
#' @param d Numeric value specifying the effect size (Cohen's d) used for group differences. Default is `0.5`.
#' @param icc Numeric value specifying the intracluster correlation coefficient used by `simulate_ge()`.
#'   Default is `0.5`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{expression}{A numeric matrix of simulated gene expression values. Rows are genes and columns
#'   are samples (barcodes).}
#'   \item{metadata}{A tibble describing each sample, including `barcode`, `group`, `compound_name`,
#'   `dose_level`, and `time_level`.}
#' }
#'
#' @examples
#' # Simulate gene expression data for 10 genes and 5 compounds in two groups
#' sim_data <- simulate_data(n_gene = 10, n_com = c(5, 5))
#' head(sim_data$expression)
#' head(sim_data$metadata)
#'
#' @export
simulate_data <- function(n_gene,
                          n_com,
                          n_dose = 3,
                          n_time = 4,
                          n_rep = 3,
                          mu = runif(1, -3, 3),
                          d = 0.5,
                          icc = 0.5
) {
  dt_mat <- replicate(n_gene,simulate_ge(n_com = n_com,
                                         n_dose = n_dose,
                                         n_time = n_time,
                                         n_rep = n_rep,
                                         mu = mu,
                                         d = d,
                                         icc = icc))
 n_samples <- sum(n_com) * n_dose * n_time * n_rep
 group <- rep(paste0("Group", 1:length(n_com)), times = n_com*n_dose * n_time * n_rep)
 compound_name <-  rep(paste0("Compound", 1:sum(n_com)), each = n_dose * n_time * n_rep)
 dose_level <-  rep(paste0("Dose", 1 : n_dose), each = n_time * n_rep, times = sum(n_com))
 time_level <-  rep(paste0("Time", 1 : n_time), each = n_rep, times = n_dose*sum(n_com))
 barcode <- paste0(rep(sapply(1:sum(n_com), add_zero), each = n_dose * n_time * n_rep),
                  rep(sapply(1 : n_dose, add_zero), each = n_time * n_rep, times = sum(n_com)),
                  rep(sapply(1 : n_time, add_zero), each = n_rep, times = n_dose*sum(n_com)),
                  rep(sapply(1 : n_rep, add_zero), times = n_dose * n_time * sum(n_com)))
 metadata <- tibble::tibble(barcode ,
                        group = factor(group, levels = paste0("Group", 1:length(n_com))),
                        compound_name = factor(compound_name, levels = paste0("Compound", 1:sum(n_com))),
                        dose_level = factor(dose_level, levels = paste0("Dose", 1 : n_dose)),
                        time_level = factor(time_level, levels =paste0("Time", 1 : n_time)))
 expr_data <- t(dt_mat)
 colnames(expr_data) <- barcode #rep(set_names("compound-dose-time", metadata, "compound_name"), each = n_rep)
 rownames(expr_data) <- paste0("Gene", 1 : n_gene)
 return(list(expression = expr_data, metadata = metadata))
}

#' Simulate toxicogenomics expression data
#'
#' @description
#' `simulate_tgxdata()` generates simulated toxicogenomics (TGx) expression data containing both
#' differentially expressed (DE) genes and equivalently expressed (EE) genes between compound groups.
#' DE genes are simulated using the specified effect size (`d`) and intracluster correlation (`icc`),
#' while EE genes are simulated with no group effect and no intracluster correlation (`d = 0`, `icc = 0`).
#' The experimental design is controlled by the number of compounds, dose levels, time points, and
#' replicates.
#'
#' @param n_de Integer specifying the number of differentially expressed (DE) genes to simulate.
#' @param n_ee Integer specifying the number of equivalently expressed (EE) genes to simulate.
#' @inheritParams simulate_data
#'
#' @return A list with three elements:
#' \item{expression}{A numeric matrix of simulated expression values for DE and EE genes combined.}
#' \item{metadata}{A tibble describing the simulated experiment (copied from the DE simulation).}
#' \item{gene_type}{A character vector indicating the type of each gene (`"DE"` or `"EE"`).}
#'
#' @examples
#' # Simulate TGx expression data with 20 DE genes and 10 EE genes
#' tgx_data <- simulate_tgxdata(n_de = 20, n_ee = 10, n_com = c(3, 3))
#'
#' @export
simulate_tgxdata <- function(n_de = 10,
                             n_ee = 10,
                             n_com = c(5,5),
                             n_dose = 3,
                             n_time = 4,
                             n_rep = 3,
                             mu = runif(1, -3, 3),
                             d = 0.5,
                             icc = 0.5) {
  de_data <- simulate_data(n_gene = n_de,
                         n_com = n_com,
                         n_dose = n_dose,
                         n_time = n_time,
                         n_rep = n_rep,
                         mu = mu,
                         d = d,
                         icc = icc)
  rownames(de_data[[1]]) <- paste0("DE", 1 : n_de)
  ee_data <- simulate_data(n_gene = n_ee,
                         n_com = n_com,
                         n_dose = n_dose,
                         n_time = n_time,
                         n_rep = n_rep,
                         mu = mu,
                         d = 0,
                         icc = 0)
  rownames(ee_data[[1]]) <- paste0("EE", 1 : n_ee)
  exprdata <- rbind(de_data[[1]], ee_data[[1]])
  gene_cl <- rep(c("DE" ,"EE"), times = c(n_de,n_ee))
  return(list(expression = exprdata, metadata = de_data[[2]], gene_type = gene_cl))
}
