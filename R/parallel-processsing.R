#' Starting parallel computing
#'
#' The function `start_parallel()` is used to setup parallel processing by appointing multiple
#' workers for a sequence of similar task.
#'
#' @param multicore Parallel computation (default is `TRUE`)
#' @param ... Additional parameter
#'
#' @return
#' Class of parallel object
#' @export
#'
#' @examples
#' \dontrun{
#' start_parallel(multicore = TRUE)
#' }
start_parallel <- function(multicore = TRUE, error_call = caller_env(),  ...) {
 if (!all(requireNamespace("parallel", quietly = TRUE),
          requireNamespace("doParallel", quietly = TRUE))) {
   cli::cli_abort(c("Packages `parallel` and `doParallel` required for parallelization!",
               "i" = "Please install `parallel` and `doParallel`."),
             call = error_call)
 }
  if (any(class(multicore) == "cluster"))
  {
    cl <- multicore
    multicore <- TRUE
    attr(multicore, "type") <- foreach::getDoParName()
    attr(multicore, "cores") <- foreach::getDoParWorkers()
    attr(multicore, "cluster") <- cl
    return(multicore)
  }
  parallel_type <- ifelse(identical(.Platform$OS.type, "windows"), "snow", "multicores")
  # get the current number of cores available (left one for system)
  numCores <- parallel::detectCores() - 1
  # set parameters for parallelization
  if (is.logical(multicore)) {
    NULL
  } else if (is.numeric(multicore)) {
    numCores <- as.integer(multicore)
    multicore <- TRUE
  } else if (is.character(multicore)) {
    parallel_type <- multicore
    multicore <- TRUE
  } else {
    multicore <- FALSE
  }
  attr(multicore, "type") <- parallel_type
  attr(multicore, "cores") <- numCores
  # start "parallel backend" if needed
  if (multicore) {
    if (parallel_type == "snow") {
      # snow functionality on Unix-like systems & Windows
      cl <- parallel::makeCluster(numCores, type = "PSOCK")
      attr(multicore, "cluster") <- cl
      # export parent environment
      varlist <- ls(envir = parent.frame(), all.names = TRUE)
      varlist <- varlist[varlist != "..."]
      parallel::clusterExport(
        cl,
        varlist = varlist,
        envir = parent.frame()
        )
      # export global environment (workspace)
      parallel::clusterExport(
        cl,
        varlist = ls(envir = globalenv(), all.names = TRUE),
        envir = globalenv()
      )
      # load current packages in workers
      pkgs <- .packages()
      distb <- parallel::clusterCall(cl, function(x) lapply(x, require, character.only = TRUE), pkgs)
      doParallel::registerDoParallel(cl, cores = numCores)
    } else if (parallel_type == "multicores") {
      # multicore functionality on Unix-like systems
      cl <- parallel::makeCluster(numCores, type = "FORK")
      doParallel::registerDoParallel(cl, cores = numCores)
      attr(multicore, "cluster") <- cl
    } else {
      cli::cli_abort(c(
        "Parallel type must be either `snow` or `multicores`.",
        "x" = "Parallel type of  {style_bold(col_red(backtick(parallelType)))} is invalid.",
        "i" = "Please use either `snow`, or `multicores` instread."
        ),
        call = error_call
      )
    }
  }
  return(multicore)
}

#' Stop parallel setup
#'
#' The function `stop_parallel()` is used to stop parallel processing.
#'
#' @param cluster Class of parallel object
#' @param ... Additional parameter
#'
#' @return
#' Exit cluster
#' @export
#'
#' @examples
#' \dontrun{
#' cll <- start_parallel()
#' stop_parallel(attr(cll, "cluster"))
#' }
stop_parallel <- function(cluster, ...)
{
  parallel::stopCluster(cluster)
  foreach::registerDoSEQ()
  invisible()
}
