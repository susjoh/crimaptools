#' parse_crossovers: Parse crossover information from CriMAP chrompic output
#' @param chrompicfile File with output from chrompic


parse_crossovers <- function(chrompicfile){
  #~~ read lines from the chrompic file

  x <- readLines(chrompicfile)

  nmarkers <- nrow(MapFromChrompic(chrompicfile))

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
  for(i in 1:nrow(famtab)) famtab$FamilyShortID[i] <- strsplit(famtab$FamilyID[i], split = " ")[[1]][2]

  x <- x[-grep("Family", x)]

  famtab$Count <- diff(c(famtab$startline, (length(x) + 1)))

  famvec <- rep(famtab$FamilyShortID, famtab$Count)

  #~~ create a data frame to put all information in

  recombframe <- data.frame(data = x,
                            id = NA,
                            RecombCount = NA,
                            parent = NA,
                            Family = famvec,
                            stringsAsFactors=F)


  recombframe$data <- gsub("^  ", " ",  recombframe$data)
  recombframe$data <- gsub("^   ", " ",  recombframe$data)

  #~~ Get the IDs


  recombframe$id <- unlist(lapply(recombframe$data, function(foo) strsplit(foo, split = " ")[[1]][2]))
  recombframe$parent[which(recombframe$id != "")] <- "MOTHER"
  recombframe$id[which(recombframe$id == "")] <- recombframe$id[which(recombframe$id != "")]
  recombframe$parent[which(is.na(recombframe$parent))] <- "FATHER"



  #~~ get the number of recombinations for each individual

  recExtract <- function(CMPstring){
    x <- strsplit(CMPstring, split = " ")[[1]][length(strsplit(CMPstring, split = " ")[[1]])]
  }

  recombframe$RecombCount <- sapply(recombframe$data, recExtract, simplify = TRUE)

  #~~ format the data column so that we can get out the number of informative loci
  removeID <- function(CMPstring){

    x <- strsplit(CMPstring, split = " ")[[1]]
    paste0(x[3:(length(x)-1)], collapse = "")

  }

  system.time(recombframe$data <- unlist(lapply(recombframe$data, removeID)))

  recombframe$data2 <- gsub("-" , "", recombframe$data)
  recombframe$data2 <- gsub("c" , "", recombframe$data2)
  recombframe$data2 <- gsub(":" , "", recombframe$data2)

  recombframe$No.Inf.Loci <- nchar(recombframe$data2)

  recombframe <- subset(recombframe, select = -data2)

  head(recombframe)

  recombframe$RecombCount <- as.numeric(recombframe$RecombCount)

  recombframe

}
