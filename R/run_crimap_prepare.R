#' run_crimap_prepare: Run prepare in crimap
#'
#' NB. This function will delete all files associated with the
#' analysis suffix, except the .gen file.
#'
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as...
#' @param genfile path for crimap .gen file.
#' @export


run_crimap_prepare <- function(crimap.path, genfile){

  #~~ parse crimap.file if in another directory

  # if(length(grep("..", crimap.file, fixed = T)) > 0) stop("This function will only work if the crimap file is in a sub-directory of the working directory")

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  #if(!crimap.path %in% dir()) stop("This function will only work if the crimap.path is in the same directory as the crimap file.")

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

  if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

  del.vec <- grep(paste0("chr", crimap.stem, "."), dir(), value = T)
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".gen"))]
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mnd"))]
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mndverbose"))]


  if(Sys.info()["sysname"] == "Windows") {

    if(length(del.vec) > 0){

      for(i in del.vec){

        system("cmd", input = paste0("del ", i), show.output.on.console = F)
      }

    }

    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"), show.output.on.console = F)

  } else {

    for(i in del.vec){

      system(paste0("rm ", i), show.output.on.console = F)
    }
    system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"), show.output.on.console = F)

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(paste(rep("..", times = length(crimap.file)), collapse = "/"))


}


