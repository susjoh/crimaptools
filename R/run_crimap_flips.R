#' run_crimap_flips: Run flips
#'
#' @param crimap.path path to run crimap. This should be given relative to the
#'   same directory as...
#' @param genfile path to the .gen file for running the chrompic function.
#' @param flips number of loci to flip.
#' @export


run_crimap_flips <- function(crimap.path, genfile, flips = 2){

  #~~ parse crimap.file if in another directory

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)


  if(Sys.info()["sysname"] == "Windows") {

    eval(
      parse(
        text = paste0("system(\"cmd\", input = \"", "\"", crimap.path, "\" ", crimap.stem, " flips", flips, " > chr", crimap.stem, ".fl", flips, "\", show.output.on.console = F)")
      )
    )


  } else {

    eval(
      parse(
        text = paste0(crimap.path, " ", crimap.stem, " flips", flips, " > chr", crimap.stem, ".fl", flips, "\", show.output.on.console = F)")
      )
    )
  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(paste(rep("..", times = length(crimap.file)), collapse = "/"))


}


