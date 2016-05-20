#' run_crimap_prepare: Run prepare in crimap
#'
#' NB. This function will delete all files associated with the
#' analysis suffix, except the .gen file.
#'
#' @param genfile path for crimap .gen file.
#' @export


run_crimap_prepare <- function(genfile, build = FALSE){


  pwd <- getwd()
  #~~ parse crimap.file if in another directory

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

  if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)
  if(build == TRUE & !file.exists("crimapinput2")) write.table(data.frame(c("n", "n", "n", "n", 1, "y", "y")), "crimapinput2", row.names = F, col.names = F, quote = F)


  del.vec <- grep(paste0("chr", crimap.stem, "."), dir(), value = T)
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".gen"))]
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mnd"))]
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mndverbose"))]


  if(Sys.info()["sysname"] == "Windows") {
    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")

    if(length(del.vec) > 0){

      for(i in del.vec){

        system("cmd", input = paste0("del ", i), show.output.on.console = F)
      }

    }

    if(build == FALSE) system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"), show.output.on.console = F)
    if(build == TRUE)  system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput2 > chr", crimap.stem, ".pre"), show.output.on.console = F)

  } else {
    crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")

    for(i in del.vec){

      system(paste0("rm ", i))
    }
    if(build == FALSE) system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"))
    if(build == TRUE)  system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput2 > chr", crimap.stem, ".pre"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


