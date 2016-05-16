#' run_crimap_chrompic: Run chrompic in crimap
#'
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as...
#' @param genfile path to the .gen file for running the chrompic function.
#' @export


run_crimap_chrompic <- function(crimap.path, genfile){

  #~~ parse crimap.file if in another directory

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)


  if(Sys.info()["sysname"] == "Windows") {

    system("cmd", input = paste0(crimap.path, " ", crimap.stem, " chrompic > chr", crimap.stem, ".cmp"), show.output.on.console = F)
    system("cmd", input = "del *.cg")


  } else {

    system(paste0(crimap.path, " ", crimap.stem, " chrompic > chr", crimap.stem, ".cmp"), show.output.on.console = F)
    system("rm *.cg")

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(paste(rep("..", times = length(crimap.file)), collapse = "/"))


}


