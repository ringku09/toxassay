toxassay_startup_message <- function()
{
  msg <- c(paste0(
"   _______
  |__   __|           /\\                        Toxicity
     | | ___ _  __   /  \\  ___ ___  __ _ __  __ Assesment
     | |/ _ \\ \\/ /  / /\\ \\/ __/ __|/ _` |\\ \\/ /
     | | (_) >  <  / ____ \\__ \\__ \\ (_| | \\  /
     |_|\\___/_/\\_\\/_/    \\____/___/\\__,_| /_/   version ", utils::packageVersion("toxassay")),
"\n\nTo cite this R package in publications, run the command: \`citation(\"toxassay\")\`.")
  return(msg)
}

.onAttach <- function(libname, pkgname)
{
  msg <- toxassay_startup_message()
  if(!interactive())
    msg[1] <- paste("Package 'toxassay' version", utils::packageVersion("toxassay"))
  packageStartupMessage(msg)
  invisible()
}
