#' run_crimap_map: Run prepare in crimap
#'
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as...
#' @param genfile path to the .gen file for running the chrompic function.
#' @export


run_crimap_map <- function(crimap.path, genfile){

  #~~ parse crimap.file if in another directory

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

  N <- read.table(crimap.file[length(crimap.file)], skip = 1, nrows = 1)[[1]]

  write.table(1, paste0("chr", crimap.stem, ".ord"), quote = F, row.names = F, col.names = F)
  write.table("", paste0("chr", crimap.stem, ".ord"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste(c(1, N), collapse = " "), paste0("chr", crimap.stem, ".ord"), quote = F, row.names = F, col.names = F, append = T)
  write.table(paste(c(0:(N - 1)), collapse = " "), paste0("chr", crimap.stem, ".ord"), quote = F, row.names = F, col.names = F, append = T)


  if(Sys.info()["sysname"] == "Windows") {

    system("cmd", input = paste0(crimap.path, " ", crimap.stem, " build > chr", crimap.stem, ".map"), show.output.on.console = F)

  } else {

    system(paste0(crimap.path, " ", crimap.stem, " build > chr", crimap.stem, ".map"), show.output.on.console = F)

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(paste(rep("..", times = length(crimap.file)), collapse = "/"))


}


