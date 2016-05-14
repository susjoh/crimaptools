crimap.path <- "crimap2504.exe"
crimap.file <- "crimap\\chr1a.gen"

run_crimap_prepare <- function(crimap.path, crimap.file){

  #~~ parse crimap.file if in another directory

  message("This function will only work if the crimap file is in a sub-directory of the working directory, and that the crimap.path is in the same directory as the crimap file.")

  crimap.file <- gsub("\\\\", "/", crimap.file)

  crimap.file <- strsplit(crimap.file, split = "/")[[1]]

  if(length(crimap.file) > 1) setwd(paste(crimap.file[1:(length(crimap.file)-1)], collapse = "/"))

  crimap.stem <- gsub("^chr", "", crimap.file[length(crimap.file)])
  crimap.stem <- gsub(".gen$", "", crimap.stem)

  if(!file.exists("crimapinput1")) write.table(data.frame(c("n", "n", "n", "n", 7, "y", "y")), "crimapinput1", row.names = F, col.names = F, quote = F)

  del.vec <- grep(paste0("chr", crimap.stem, "."), dir(), value = T)
  del.vec <- del.vec[-which(del.vec == paste0("chr", crimap.stem, ".gen"))]

  if(Sys.info()["sysname"] == "Windows") {

    if(length(del.vec) > 0){

      for(i in del.vec){

        system("cmd", input = paste0("del ", i))
      }

    }

    system("cmd", input = paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"))

  } else {

    for(i in del.vec){

      system(paste0("rm ", i))
    }

  }

  system(paste0(crimap.path, " ", crimap.stem, " prepare < crimapinput1 > chr", crimap.stem, ".pre"))

  crimap.file  <- crimap.file[-length(crimap.file)]

  setwd(paste(rep("..", times = length(crimap.file)), collapse = "/"))


}


