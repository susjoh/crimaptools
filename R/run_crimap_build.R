#' run_crimap_build: Run chrompic in crimap
#'
#' @param genfile path to the .gen file for running the chrompic function.
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as the genfile. Non-windows only at present.
#' @export


run_crimap_build <- function(genfile, crimap.path = NULL){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

  if(!paste0("chr", crimap.stem, ".ord") %in% dir()) stop("Requires ord file - use run_crimap_prepare with build = TRUE")

  if(Sys.info()["sysname"] == "Windows") {
    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")

    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " build > chr", crimap.stem, ".bld"), show.output.on.console = F)


  } else {

    crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")

    system(paste0(crimap.path, " ", crimap.stem, " build > chr", crimap.stem, ".bld"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


