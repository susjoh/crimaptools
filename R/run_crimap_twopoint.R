#' run_crimap_twopoint: Run twopoint in crimap
#'
#' @param genfile path to the .gen file for running the twopoint function.
#' @export


run_crimap_twopoint <- function(genfile, crimap.path = NULL){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)


  if(Sys.info()["sysname"] == "Windows") {
    if(is.null(crimap.path)){
      crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")
    }
    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " twopoint > chr", crimap.stem, ".tpt"), show.output.on.console = F)


  } else {
    if(Sys.info()["sysname"] == "Linux"){
      crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/linux/crimap")
      if(!file.exists(crimap.path)){
        crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")
      }
    } else {
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/macos/crimap")
    }
    system(paste0(crimap.path, " ", crimap.stem, " twopoint > chr", crimap.stem, ".tpt"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


