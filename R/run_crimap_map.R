#' run_crimap_map: Run prepare in crimap
#'
#' @param genfile path to the .gen file for running the chrompic function.
#' @export


run_crimap_map <- function(genfile, crimap.path = NULL){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

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
    if(is.null(crimap.path)){
      crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")
    }    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " build > chr", crimap.stem, ".map"), show.output.on.console = F)

  } else {
    if(Sys.info()["sysname"] == "Linux"){
      crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/linux/crimap")
      if(!file.exists(crimap.path)){
        crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")
      }
    } else {
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/macos/crimap")
    }
    system(paste0(crimap.path, " ", crimap.stem, " build > chr", crimap.stem, ".map"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


