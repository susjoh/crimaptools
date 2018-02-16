#' run_crimap_twopoint: Run twopoint in crimap
#'
#' @param genfile path to the .gen file for running the twopoint function.
#' @export


run_crimap_twopoint <- function(genfile){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)


  if(Sys.info()["sysname"] == "Windows") {
    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")

    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " chrompic > chr", crimap.stem, ".cmp"), show.output.on.console = F)
    system("cmd", input = paste0("del chr", crimap.stem, "_*.cg"), show.output.on.console = F)


  } else {
    if(Sys.info()["sysname"] == "Linux"){
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")
    } else {
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/macos/crimap")
    }
    system(paste0(crimap.path, " ", crimap.stem, " twopoint > chr", crimap.stem, ".tpt"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


