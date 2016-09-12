#' run_crimap_prepare: Run prepare in crimap
#'
#' NB. This function will delete all files associated with the
#' analysis suffix, except the .gen file.
#'
#' @param genfile path for crimap .gen file.
#' @param build logical. Tells crimap whether to prepare for build analysis.
#' @param snplist character vector. Ordered loci for build analysis.
#' @param snpinset character vector. Loci to be added to the build.
#' @import plyr
#' @export



run_crimap_prepare <- function(genfile, build = FALSE, snplist = NULL, snpinsert = NULL){

  if(!is.null(snplist) & build == FALSE){
    message("snplist is specified - changed to build = TRUE")
    build <- TRUE
  }

  if(!is.null(snpinsert) & is.null(snplist)){
    stop("snplist must be specified if snpinsert is specified.")
  }

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

  if(paste0("chr", crimap.stem, ".mnd") %in% del.vec){

    del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mnd"))]

  }

  if(paste0("chr", crimap.stem, ".mndverbose") %in% del.vec){

    del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".mndverbose"))]

  }

  if(!is.null(snplist)){

    #~~ read the SNP loci from the gen file to get indices



    x <- readLines(paste0("chr", crimap.stem, ".gen"))
    x.order <- x[4:(as.numeric(x[2])+3)]
    rm(x)

    x.order <- data.frame(SNP.Name = x.order,
                          Index = 0:(length(x.order)-1))

    snplist.work <- data.frame(SNP.Name = snplist)
    suppressMessages(snplist.work <- join(snplist.work, x.order))
    snplist.index <- snplist.work$Index
    snplist.index <- paste(c(snplist.work$Index, "*"), collapse = " ")


    if(!is.null(snpinsert)){
      snpinsert.work <- data.frame(SNP.Name = snpinsert)
      suppressMessages(snpinsert.work <- join(snpinsert.work, x.order))
      snpinsert.index <- paste(c(snpinsert.work$Index, "*"), collapse = " ")
    } else {

      snpinsert.index <- "*"
    }

  }

  #~~ Create the crimap input for the analysis

  if(build == TRUE & !is.null(snplist)){

    write.table(data.frame(c("n", "n", "n", "n", 1, "n", snplist.index, snpinsert.index, "y", "y")),
                "crimapinput3", row.names = F, col.names = F, quote = F)


  }

  if(Sys.info()["sysname"] == "Windows") {
    crimap.path <- paste0(.libPaths()[1], "/crimaptools/bin/windows64/crimap2504.exe")

    if(length(del.vec) > 0){

      for(i in del.vec){

        system("cmd", input = paste0("del ", i), show.output.on.console = F)
      }

    }

    if(build == FALSE) system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"), show.output.on.console = F)
    if(build == TRUE & is.null(snplist))  system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput2 > chr", crimap.stem, ".pre"), show.output.on.console = F)
    if(build == TRUE & !is.null(snplist))  system("cmd", input = paste0("\"", crimap.path, "\" ", crimap.stem, " prepare < crimapinput3 > chr", crimap.stem, ".pre"), show.output.on.console = F)

  } else {
    
    if(Sys.info()["sysname"] == "Linux"){
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/linux/crimap")
    } else {
      crimap.path <- paste0(.libPaths()[length(.libPaths())], "/crimaptools/bin/macos/crimap")
    }

    for(i in del.vec){

      system(paste0("rm ", i))
    }
    if(build == FALSE) system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"))
    if(build == TRUE & is.null(snplist))  system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput2 > chr", crimap.stem, ".pre"))
    if(build == TRUE & !is.null(snplist))  system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput3 > chr", crimap.stem, ".pre"))

  }


  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(pwd)


}


