#' run_crimap_build: Run chrompic in crimap
#'
#' @param genfile path to the .gen file for running the chrompic function.
#' @param PLT numeric. Optional. Specify PUK_LIKE_TOL and PK_LIKE_TOL value.
#' @export


run_crimap_build <- function(genfile, PLT = NULL){

  #~~ parse crimap.file if in another directory
  pwd <- getwd()

  crimap.file <- gsub("\\\\", "/", genfile)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

#   if(!paste0("chr", crimap.stem, ".ord") %in% dir()){
#
#     setwd(pwd)
#     stop("Requires ord file - use run_crimap_prepare with build = TRUE")
#
#   }

  bld.run <- 1

  if(length(grep(paste0("chr", crimap.stem, ".bld"), dir())) > 0){
    prev.run <- grep(paste0("chr", crimap.stem, ".bld"), dir(), value = T)
    prev.run <- as.numeric(substr(prev.run, nchar(prev.run), nchar(prev.run)))
    bld.run <- max(prev.run) + 1

    if(!is.null(PLT)){

      par.file <- readLines(paste0("chr", crimap.stem, ".par"))
      par.file[grep("_LIKE_TOL", par.file)] <- c(paste0("PUK_LIKE_TOL  ", PLT, " *"),
                                                 paste0("PK_LIKE_TOL  ", PLT, " *"))
      writeLines(par.file, paste0("chr", crimap.stem, ".par"))
    }

  }


  if(Sys.info()["sysname"] == "Windows") {

    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")
    system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " build > chr", crimap.stem, ".bld", bld.run), show.output.on.console = F)


  } else {

    if(Sys.info()["sysname"] == "Linux"){
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")
    } else {
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/macos/crimap")
    }
    system(paste0(crimap.path, " ", crimap.stem, " build > chr", crimap.stem, ".bld", bld.run))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


