#' run_crimap_flips: Run flips
#'
#' @param genfile path to the .gen file for running the chrompic function.
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as the genfile. Non-windows only at present.
#' @param flips number of loci to flip.
#' @export


run_crimap_flips <- function(genfile, crimap.path = NULL, flips = 2){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)


  if(Sys.info()["sysname"] == "Windows") {

    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")

    eval(
      parse(
        text = paste0("system(\"cmd\", input = \"", "\"", crimap.path, "\" ", crimap.stem, " flips", flips, " > chr", crimap.stem, ".fl", flips, "\", show.output.on.console = F)")
      )
    )


  } else {

    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/linux/crimap")

    eval(
      parse(
        text = paste0(crimap.path, " ", crimap.stem, " flips", flips, " > chr", crimap.stem, ".fl", flips, "\")")
      )
    )
  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


