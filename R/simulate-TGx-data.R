#' Simulate Variance Components for ICC
#'
#' This function generates simulated variance components based on a specified intracluster
#' correlation coefficient (ICC). It can simulate either equal or unequal variance components
#' across a specified number of components.
#'
#' @param icc Numeric value representing the intracluster correlation coefficient (ICC).
#' @param n_component Integer value indicating the number of variance components to simulate. Default is 3.
#' @param equal_var Logical value indicating whether the variance components should be equal (`TRUE`) or unequal (`FALSE`). Default is `TRUE`.
#'
#' @return A numeric vector of length `n_component` representing the simulated variance components, normalized by the ICC.
#'
#' @examples
#' # Simulate 3 variance components with equal variance
#' sim_variance(icc = 0.5, n_component = 3, equal_var = TRUE)
#'
#' # Simulate 3 variance components with unequal variance
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

#' Simulate Gene Expression Data for Toxicogenomics Studies
#'
#' This function generates simulated gene expression data based on specified parameters,
#' such as the number of compounds, doses, time points, and replicates. It also incorporates
#' a given effect size, intracluster correlation coefficient (ICC), and error variance.
#'
#' @param n_com Numeric vector of length 2 indicating the number of compounds in each group
#'  (e.g., toxic and non-toxic). Default is `c(5, 5)`.
#' @param n_dose Integer specifying the number of dose levels. Default is `3`.
#' @param n_time Integer specifying the number of time points. Default is `4`.
#' @param n_rep Integer specifying the number of replicates for each condition. Default is `3`.
#' @param mu Numeric value representing the overall mean of the gene expression. Default is a random value between -3 and 3.
#' @param d Numeric value for the effect size (Cohen's d) between the two groups. Default is `0.5`.
#' @param icc Numeric value for the intracluster correlation coefficient. It should be between 0 and 1. Default is `0.5`.
#' @param error_call Environment used to capture and handle errors with informative messages. Default is `caller_env()`.
#'
#' @return A numeric vector of simulated gene expression values.
#'
#' @details
#' The function simulates gene expression data by first generating variance components based on the provided ICC.
#' It then generates compound, dose, and time point effects using normal distributions.
#'
#' @examples
#' # Simulate gene expression data with default settings
#' simulated_data <- simulate_ge(n_com = c(3, 3))
#'
#' # Simulate gene expression data with specific parameters
#' simulated_data <- simulate_ge(n_com = c(3, 3))
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
  if (length(n_com) > 2) {
    cli::cli_abort(c("Simulation only support two groups of compound.",
                "x" = "You have supplied {style_bold(col_red(n_com))} compound groups.",
                "i" = "Please provide number of compound groups of length 2."),
              wrap = TRUE, call = error_call)
  }
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

#' Simulate Gene Expression Data
#'
#' This function generates simulated gene expression data across multiple conditions,
#' including different compounds, dose levels, and time points.
#'
#' @param n_gene Integer. Number of genes to simulate.
#' @param n_com Integer vector. Number of compounds in each group.
#' @param n_dose Integer. Number of dose levels. Default is 3.
#' @param n_time Integer. Number of time points. Default is 4.
#' @param n_rep Integer. Number of replicates for each condition. Default is 3.
#' @param mu Numeric. Mean expression level for the simulated genes. Default is a random value between -3 and 3.
#' @param d Numeric. Effect size for differential expression. Default is 0.5.
#' @param icc Numeric. Intracluster correlation coefficient. Default is 0.5.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{expression}{A matrix of simulated gene expression data. Rows represent genes and columns represent samples.}
#'   \item{metadata}{A data frame containing metadata for each sample, including compound, dose level, time point, and replicate information.}
#' }
#'
#' @examples
#' # Simulate gene expression data for 10 genes for 5 compounds in two groups.
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


#' Simulate Toxicogenomics Expression Data
#'
#' This function generates simulated toxicogenomics (TGx) expression data, including differentially expressed (DE) genes
#' and equivalently expressed (EE) genes between groups. The function allows for customization of the number of DE and EE genes,
#' as well as various experimental conditions such as the number of compounds, doses, time points, and replicates.
#'
#' @param n_de Integer. The number of differentially expressed (DE) genes to simulate.
#' @param n_ee Integer. The number of equivalently expressed (EE) genes to simulate.
#' @inheritParams simulate_data
#'
#' @return A list containing:
#' \item{expression}{A matrix of simulated expression data for both DE and EE genes.}
#' \item{metadata}{A data frame of metadata associated with the simulated experiment.}
#' \item{gene_cl}{A vector indicating the class of each gene (DE or EE).}
#'
#' @examples
#' # Simulate TGx expression data with 20DEGs and 10 EEGs.
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

# xx <- simulate_data(n_gene = 1000,
#                        n_com = c(5, 5),
#                        n_dose = 3,
#                        n_time = 4,
#                        n_rep = 3,
#                        mu = runif(1, -3, 3),
#                        d = 0.2,
#                        icc = 0.1)
#
# gr <- list(A = paste0("Compound", 1:5), B = paste0("Compound", 6:10))
# vv <- de_genes(gr, ge_matrix = xx$expression, metadata = xx$metadata, p_adjust = TRUE)
# sum(p.adjust(ttest_degs(xx$expression, xx$metadata), method ="bonferroni") < 0.05)/1000
# sum(vv$p_value < 0.05)/1000

