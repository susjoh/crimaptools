#' parse_crossovers: Parse crossover information from CriMAP chrompic output
#' @param chrompicfile File with output from chrompic
#' @param familyPedigree data.frame containing columns ANIMAL, FATHER, MOTHER
#'   and FAMILY. FAMILY defines family groups for crimap. IDs in ANIMAL can be
#'   repeated within and between families if necessary. Should be that used in
#'   create_crimap_input.
#' @param remove.zero.inf.loci Logical, default = TRUE. Remove IDs with no
#'   informative loci.
#' @import plyr
#' @export


parse_crossovers <- function(chrompicfile, familyPedigree, remove.zero.inf.loci = TRUE){
  #~~ read lines from the chrompic file

  x <- readLines(chrompicfile)

  nmarkers <- nrow(parse_map_chrompic(chrompicfile))

  x <- x[1:grep("Sex_averaged", x)]

  rowstokeep <- sort(c(grep("^ ", x),
                       grep("phase likelihood", x)))

  x <- x[rowstokeep]
  if(nmarkers < 101) x <- x[-c(1:nmarkers)]
  if(nmarkers > 100) x <- x[-c(1:100)]


  x <- x[-grep("         ", x)]

  #~~ create a table with family information to add to the full table later

  famtab <- data.frame(FamilyOrder = 1:length(grep("Family", x)),
                       FamilyID = grep("Family", x, value = T),
                       LineNo = grep("Family", x),
                       stringsAsFactors = F)

  famtab$startline <- famtab$LineNo - (famtab$FamilyOrder - 1)
  famtab$FamilyShortID <- NA

  famtab$FamilyShortID <- unlist(lapply(famtab$FamilyID, function (x) strsplit(x, split = " ")[[1]][2]))
  x <- x[-grep("Family", x)]

  famtab$Count <- diff(c(famtab$startline, (length(x) + 1)))

  famvec <- rep(famtab$FamilyShortID, famtab$Count)

  #~~ create a data frame to put all information in

  recombframe <- data.frame(data = x,
                            ANIMAL = NA,
                            RecombCount = NA,
                            parent = NA,
                            Family = famvec,
                            stringsAsFactors=F)

  


  #~~ Get the IDs

  recombframe$ANIMAL <- unlist(lapply(recombframe$data, function(foo){
    z <- strsplit(foo, split = "\\s+")[[1]]
    if(length(z) == ceiling((nmarkers/10)+2)){
      return(z[1])
    } else {
      return(z[2])
    }
  }))
  
  recombframe$parent[which(recombframe$ANIMAL != "")] <- "MOTHER"
  recombframe$ANIMAL[which(recombframe$ANIMAL == "")] <- recombframe$ANIMAL[which(recombframe$ANIMAL != "")]
  recombframe$parent[which(is.na(recombframe$parent))] <- "FATHER"



  #~~ get the number of recombinations for each individual

  
  recExtract <- function(CMPstring){
    x <- strsplit(CMPstring, split = " ")[[1]][length(strsplit(CMPstring, split = " ")[[1]])]
  }

  recombframe$RecombCount <- sapply(recombframe$data, recExtract, simplify = TRUE)

  #~~ format the data column so that we can get out the number of informative loci
  removeID <- function(CMPstring){

    z <- strsplit(CMPstring, split = "\\s+")[[1]]
    if(length(z) == ceiling((nmarkers/10)+2)){
      return(paste0(z[2:(length(z)-1)], collapse = ""))
    } else {
      return(paste0(z[3:(length(z)-1)], collapse = ""))
    }

  }

  system.time(recombframe$data <- unlist(lapply(recombframe$data, removeID)))

  recombframe$data2 <- gsub("-" , "", recombframe$data)
  recombframe$data2 <- gsub("c" , "", recombframe$data2)
  recombframe$data2 <- gsub(":" , "", recombframe$data2)

  recombframe$No.Inf.Loci <- nchar(recombframe$data2)

  recombframe <- subset(recombframe, select = -data2)

  recombframe$RecombCount <- as.numeric(recombframe$RecombCount)

  if(remove.zero.inf.loci == TRUE) recombframe <- subset(recombframe, No.Inf.Loci != 0)

  #~~ Get first and list informative positions

  recombframe$First.Inf.Order <- unlist(lapply(recombframe$data, InfLengthFunc))
  recombframe$Last.Inf.Order <-  unlist(lapply(recombframe$data, function(x) InfLengthFunc(x, position = "Last")))

  #~~ add the RRID information - the individual in which the recombination took place

  suppressMessages(recombframe <- join(recombframe, familyPedigree))
  recombframe$RRID <- NA
  recombframe$RRID[which(recombframe$parent == "MOTHER")] <- recombframe$MOTHER[which(recombframe$parent == "MOTHER")]
  recombframe$RRID[which(recombframe$parent == "FATHER")] <- recombframe$FATHER[which(recombframe$parent == "FATHER")]

  #~~ add analysisID

  analysisID.val <- gsub("\\\\", "/", chrompicfile)

  analysisID.val <- strsplit(analysisID.val, split = "/")[[1]]
  analysisID.val <- analysisID.val[length(analysisID.val)]
  analysisID.val <- gsub("chr", "", analysisID.val)
  analysisID.val <- gsub(".cmp", "", analysisID.val, fixed = T)

  recombframe$analysisID <- analysisID.val

  recombframe$UniqueID <- paste(recombframe$analysisID, recombframe$Family, recombframe$RRID, recombframe$parent, sep = "_")


  recombframe

}
