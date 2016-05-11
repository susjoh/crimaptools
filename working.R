# Prerequisites (need to be introduced for data checking)
# pedigree must be in format ANIMAL FATHER MOTHER with those headers
# unix tools required on system: fgrep, sed



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ genabelSubset: Make a chromosome object with pedigree
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



# FUNCTION: Make a subset of a gwaa.data object by chromosome and ID.

# INPUT:
# gwaa.data: GenABEL gwaa.data object
# chr      : Chromosome Number(s) as vector
# pedigree : data.frame with pedigree with column names c("ANIMAL", "FATHER", "MOTHER")
# sampleno : (default = NULL) number of SNP loci to sample from the dataset
# ordsample: (default = TRUE) whether to order sampled SNP loci by their position on the chromosome

# OUTPUT:
# chrobj: GenABEL gwaa.data object


genabelSubset <- function(gwaa.data, chr = NULL, pedigree = NULL, sampleno = NULL, snplist = NULL, ordsample = TRUE){

  require(GenABEL)

  chrobj <- gwaa.data

  if(!is.null(chr)){

    chrobj <- gwaa.data[,which(chromosome(gwaa.data) %in% chr)]    #~~ Subset the data by Chromosome

  }

  if(!is.null(pedigree)){

    chrobj <- chrobj[which(idnames(chrobj) %in% pedigree[,1]),]  #~~ subset the data further by IDs in the pedigree object

  }

  #~~ extract snps if list is specified

  if(!is.null(snplist)){

    chrobj <- chrobj[,snplist]

  }

  #~~ conduct SNP sampling if specified

  if(!is.null(sampleno)){

    subsnps <- sample(1:nsnps(chrobj), size = sampleno, replace=F)

    if(ordsample == TRUE) subsnps <- sort(subsnps)

    chrobj <- chrobj[,subsnps]

  }



  #~~ return the new gwaa.data object

  chrobj

}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# genabel2PreCrimap: Create a data frame in pre-Crimap format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


genabel2PreCrimap <- function(gwaa.data, pedigree){

  pedigree <- pedigree[,1:3]
  # FUNCTION: Create a data frame in pre-Crimap format

  # INPUT:
  # gwaa.data: GenABEL gwaa.data object
  # pedigree : data.frame with pedigree with column names c("ANIMAL", "FATHER", "MOTHER")


  # OUTPUT:
  # test: data.frame in pre-Crimap format

  require(GenABEL)
  require(data.table)

  #~~ select only individuals which are within the pedigree provided

  ids <- as.character(unique(c(pedigree[,1], pedigree[,2], pedigree[,3])))
  ids <- ids[which(ids %in% idnames(gwaa.data))]

  gwaa.data <- gwaa.data[ids,]

  #~~ Recode ACGT as 1234 with space delimited between alleles

  test <- as.character.snp.data(gtdata(gwaa.data))
  test <- as.matrix(test)

  print("Recoding alleles to numeric values...")

  test[which(is.na(test))] <- "0 0"
  test[which(test == "A/A")] <- "1 1"
  test[which(test == "A/C")] <- "1 2"
  test[which(test == "A/G")] <- "1 3"
  test[which(test == "A/T")] <- "1 4"
  test[which(test == "C/A")] <- "2 1"
  test[which(test == "C/C")] <- "2 2"
  test[which(test == "C/G")] <- "2 3"
  test[which(test == "C/T")] <- "2 4"
  test[which(test == "G/A")] <- "3 1"
  test[which(test == "G/C")] <- "3 2"
  test[which(test == "G/G")] <- "3 3"
  test[which(test == "G/T")] <- "3 4"
  test[which(test == "T/A")] <- "4 1"
  test[which(test == "T/C")] <- "4 2"
  test[which(test == "T/G")] <- "4 3"
  test[which(test == "T/T")] <- "4 4"
  test[which(test == "1/1")] <- "1 1"
  test[which(test == "1/2")] <- "1 2"
  test[which(test == "2/1")] <- "2 1"
  test[which(test == "2/2")] <- "2 2"



  print("...done.")


  #~~ Create a data frame with ID, FATHER, MOTHER, SEX, and Genotypes.

  test <- data.frame(test, stringsAsFactors=F)

  test$ANIMAL <- idnames(gwaa.data)
  test$SEX <- phdata(gwaa.data)$sex


  pedigree$ANIMAL <- as.character(pedigree$ANIMAL)
  test2 <- data.table(test, key = "ANIMAL")
  ped2 <- data.table(pedigree, key = "ANIMAL")


  print("Merging pedigree and genotype information...")
  system.time(test2 <- merge(test2, ped2, all.y = T))
  print("...done.")

  test2 <- data.frame(test2)
  test2 <- unique(test2)


  system.time(test3 <- cbind(test2[,c("ANIMAL", "MOTHER", "FATHER", "SEX")], test2[,which(!names(test2) %in% c("ANIMAL", "MOTHER", "FATHER", "SEX"))]))

  #~~ put in missing sexes

  #   test3$SEX[which(test3$ANIMAL %in% pedigree$MOTHER)] <- 0
  #   test3$SEX[which(test3$ANIMAL %in% pedigree$FATHER)] <- 1
  #   test3$SEX[which(is.na(test3$SEX))] <- 3
  #
  #~~ Deal with NA's

  corner(test3)
  test3$MOTHER[which(is.na(test3$MOTHER))] <- 0
  test3$FATHER[which(is.na(test3$FATHER))] <- 0

  for(i in 1:ncol(test3)) test3[which(is.na(test3[,i])),i] <- "0 0"

  #~~ return object

  test3
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create Crimap Files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# gwaa.data <- genabelSubset(abeldata, chr = 24)
# chromosomeid <- paste(i, "c", sep = "")
# familyPedigree <- famped
# chr = 24
# snp.list <- snp.names(gwaa.data)[c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9)]

createCrimapFiles <- function(gwaa.data, chromosomeid, familyPedigree, pseudoautoSNPs, menderrtab = NULL, chr = NULL, outdir = NULL) {

  nfamilies <- length(unique(familyPedigree$Family))
  nloci <- nsnps(gwaa.data)
  locus.names <- snp.names(gwaa.data)

  if(!is.null(outdir))  outfile <- paste0(outdir, "/chr", chromosomeid, ".gen") else outfile <- paste0("chr", chromosomeid, ".gen")

  write.table(nfamilies, outfile, row.names = F, quote = F, col.names = F)
  write.table(nloci, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(locus.names, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)


  #~~ create a master file of genotypes
  system.time(genotab <- genabel2PreCrimap(gwaa.data, pedigree=familyPedigree[,1:3]))

  genotab <- unique(genotab)

  genotab <- merge(familyPedigree, genotab, all.x = T)

  genotab <- genotab[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", names(genotab)[6:ncol(genotab)])]


  if(length(grep("^X", chromosomeid)) > 0){
    xmarkers <- which(!names(genotab) %in% pseudoautoSNPs)
    xmarkers <- xmarkers[6:length(xmarkers)]

    heterozygotes <- c("1 2", "1 3", "1 4", "2 1", "2 3", "2 4", "3 1", "3 2", "3 4", "4 1", "4 2", "4 3")

    for(j in xmarkers){

      #~~ remove heterozygotes

      genotab[which(genotab$SEX == 1 & genotab[,j] %in% heterozygotes),j] <- "0 0"

      genotab[which(genotab$SEX == 1 & genotab[,j] == "1 1"),j] <- "1 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "2 2"),j] <- "2 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "3 3"),j] <- "3 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "4 4"),j] <- "4 0"

    }

    rm(xmarkers, heterozygotes)

  }


  #~~ deal with mendelian errors

  if(!is.null(menderrtab)){
    menderrtab <- subset(menderrtab, Chr == chr)
    if(nrow(menderrtab) > 0){
      for(i in 1:nrow(menderrtab)){
        genotab[which(genotab$ANIMAL == menderrtab$ANIMAL[i]), 6 + as.numeric(menderrtab$Locus[i])] <- "0 0"
      }
    }
  }

  genotab <- as.matrix(genotab)


  counter <- 0

  system.time({

    #~~ determine lengths of families

    test <- data.frame(FamSize = tapply(familyPedigree$Family, familyPedigree$Family, length))
    test$Family <- row.names(test)

    familyPedigree <- arrange(familyPedigree, Family)

    all(unique(test$Family) == unique(familyPedigree$Family))

    final.frame <- matrix("", ncol = (ncol(genotab)-1), nrow = (sum(test$FamSize) + (nrow(test)*4)))



    counter <- 1

    for(i in unique(familyPedigree$Family)){

      famsize <- test$FamSize[which(test$Family == i)]

      final.frame[counter, 1] <- i
      final.frame[counter+1, 1] <- famsize
      final.frame[(counter+3):(counter+2+famsize),] <- genotab[which(genotab[,"Family"] == i),-1]

      counter <- counter + 4 + famsize

    }
  })


  final.frame2 <- apply(final.frame, 1, function (x) paste(x, collapse = " "))
  final.frame2 <- gsub("\\s+$", "", final.frame2)

  write.table(final.frame2, outfile, row.names = F, quote = F, col.names = F, append=T)


}



createCrimapFilesList <- function(gwaa.data, chromosomeid, familyPedigree, pseudoautoSNPs, snp.list = NULL, chr = NULL, outdir = NULL) {

  print("This function does not deal with Mendelian errors.")

  if(!is.null(snp.list)) gwaa.data <- gwaa.data[,snp.list]
  nfamilies <- length(unique(familyPedigree$Family))
  nloci <- nsnps(gwaa.data)
  locus.names <- snp.names(gwaa.data)

  if(!is.null(outdir))  outfile <- paste0(outdir, "/chr", chromosomeid, ".gen") else outfile <- paste0("chr", chromosomeid, ".gen")

  write.table(nfamilies, outfile, row.names = F, quote = F, col.names = F)
  write.table(nloci, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(locus.names, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)


  #~~ create a master file of genotypes
  system.time(genotab <- genabel2PreCrimap(gwaa.data, pedigree=familyPedigree[,1:3]))

  genotab <- unique(genotab)

  genotab <- merge(familyPedigree, genotab, all.x = T)

  if(is.null(snp.list)){
    genotab <- genotab[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", names(genotab)[6:ncol(genotab)])]
  }

  if(!is.null(snp.list)){
    genotab <- genotab[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", snp.list)]
  }


  if(length(grep("^X", chromosomeid)) > 0){
    xmarkers <- which(!names(genotab) %in% pseudoautoSNPs)
    xmarkers <- xmarkers[6:length(xmarkers)]

    heterozygotes <- c("1 2", "1 3", "1 4", "2 1", "2 3", "2 4", "3 1", "3 2", "3 4", "4 1", "4 2", "4 3")

    for(j in xmarkers){

      #~~ remove heterozygotes

      genotab[which(genotab$SEX == 1 & genotab[,j] %in% heterozygotes),j] <- "0 0"

      genotab[which(genotab$SEX == 1 & genotab[,j] == "1 1"),j] <- "1 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "2 2"),j] <- "2 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "3 3"),j] <- "3 0"
      genotab[which(genotab$SEX == 1 & genotab[,j] == "4 4"),j] <- "4 0"

    }

    rm(xmarkers, heterozygotes)

  }



  genotab <- as.matrix(genotab)

  print("Writing crimap file...")

  system.time({

    #~~ determine lengths of families

    test <- data.frame(FamSize = tapply(familyPedigree$Family, familyPedigree$Family, length))
    test$Family <- row.names(test)

    familyPedigree <- arrange(familyPedigree, Family)

    all(unique(test$Family) == unique(familyPedigree$Family))

    final.frame <- matrix("", ncol = (ncol(genotab)-1), nrow = (sum(test$FamSize) + (nrow(test)*4)))



    counter <- 1

    for(i in unique(familyPedigree$Family)){

      famsize <- test$FamSize[which(test$Family == i)]

      final.frame[counter, 1] <- i
      final.frame[counter+1, 1] <- famsize
      final.frame[(counter+3):(counter+2+famsize),] <- genotab[which(genotab[,"Family"] == i),-1]

      counter <- counter + 4 + famsize

    }
  })


  final.frame2 <- apply(final.frame, 1, function (x) paste(x, collapse = " "))
  final.frame2 <- gsub("\\s+$", "", final.frame2)

  write.table(final.frame2, outfile, row.names = F, quote = F, col.names = F, append=T)

  print("...done.")

}


#~~ Parse flips files

flipsToTable <- function(flipsfile){

  test2 <- readLines(flipsfile)
  if(length(grep("LOG_LIKE DECREASED IN MLE", test2) > 0)) test2 <- test2[-grep("LOG_LIKE DECREASED IN MLE", test2)]
  test2 <- test2[-grep("^$", test2)]
  test2 <- test2[(grep("(= log10_like[orig] - log10_like[curr])", test2, fixed = T)+1):length(test2)]


  z <- lapply(test2, function(x){
    y <- strsplit(x, split = "\\s+")[[1]][-1]
    y <- matrix(y, ncol = length(y))
    data.frame(y)
  }
  )

  data.table::rbindlist(z)
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Pull Map from Chrompic                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


MapFromChrompic <- function(chrompicfile){

  x <- readLines(chrompicfile)
  map <- x[grep("Sex_averaged", x):length(x)]
  map <- map[3:length(map)]
  map <- map[1:(length(map)-6)]

  map <- data.frame(map, stringsAsFactors=F)

  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]

  for(i in 1:10) map$map  <- gsub("  ", " ", map$map)
  for(i in 1:10) map$map2 <- gsub("  ", " ", map$map2)

  rm(x)

  map$Order <- NA
  map$SNP.Name <- NA
  map$Position <- NA
  map$r <- NA
  map$cMdiff <- NA

  for(i in 1:nrow(map)){
    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Order[i] <- x[1]
    map$SNP.Name[i] <- x[2]
    map$Position[i] <- x[3]
    rm(x)

    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$r[i] <- x[1]
    map$cMdiff[i] <- x[2]

    rm(x)
  }

  map <- subset(map, select = -c(map, map2))

  map


}



MapFromMap <- function(mapfile){

  x <- readLines(mapfile)
  map <- x[grep("Sex-specific map", x):length(x)]


  map <- map[3:length(map)]
  map <- map[1:(length(map)-7)]
  map <- map[-which(map == "")]

  map <- data.frame(map, stringsAsFactors=F)

  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]

  for(i in 1:1e4) if(length(grep("  ", map$map)  > 0)) map$map   <- gsub("  ", " ", map$map)  else break
  for(i in 1:1e4) if(length(grep("  ", map$map2) > 0)) map$map2  <- gsub("  ", " ", map$map2) else break

  rm(x)

  map$Order <- NA
  map$SNP.Name <- NA
  map$cMPosition.Female <- NA
  map$cMPosition.Male <- NA
  map$Female.r <- NA
  map$cMdiff.Female <- NA
  map$Male.r <- NA
  map$cMdiff.Male <- NA

  for(i in 1:nrow(map)){

    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Order[i]             <- x[1]
    map$SNP.Name[i]          <- x[2]
    map$cMPosition.Female[i] <- x[3]
    map$cMPosition.Male[i]   <- x[4]

    rm(x)

    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Female.r[i]      <- x[1]
    map$cMdiff.Female[i] <- x[2]
    map$Male.r[i]        <- x[3]
    map$cMdiff.Male[i]   <- x[4]


    rm(x)
  }

  map <- subset(map, select = -c(map, map2))

  map
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Determine the recombination rate form chrompic files       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


RecombRateFromChrompic <- function(chrompicfile){
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











#' Format Pedigree
#' This function formats the pedigree for downstream analysis.
#' @param ped Pedigree object in "simple" format (Three columns for ANIMAL,
#' MOTHER and FATHER) or in "plink" format (Five to Six columns for FAMILY,
#' ANIMAL, FATHER, MOTHER, SEX and Phenotype, where the phenotype column is
#' optional). The simple argument can recognise the order of parents if
#' they are named sensibly. Run simple.ped.name.rules() for an example.
#' @param pedigree.type Defaults to "simple", can also accept "plink" which
#' is equivalent for for first 5 to 6 columns of a PLINK .ped file.


pedigree.format <- function(ped, pedigree.type = "simple"){   # "plink"

  simple.ped.name.rules <- function(){
    writeLines("Pedigree columns must be named as follows:
               ID should be named ID or ANIMAL
               Mother should be MOTHER, MUM, MOM or DAM
               Father should be FATHER, DAD, POP or SIRE")
  }

  if(pedigree.type == "simple"){

    ped <- ped[,1:3]
    names(ped) <- toupper(names(ped))
    if(!all(names(ped) %in% c("ID", "ANIMAL", "MUM", "MOM",
                              "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE"))) stop(simple.ped.name.rules())

    names(ped)[which(names(ped) %in% c("ID", "ANIMAL"))] <- "ANIMAL"
    names(ped)[which(names(ped) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"
    names(ped)[which(names(ped) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"

    ped <- ped[,c("ANIMAL", "MOTHER", "FATHER")]

    ped[,1:3] <- lapply(ped[,1:3], as.character)


    for(i in 1:3) ped[which(is.na(ped[,i])),i] <- 0

    if(any(!ped$MOTHER %in% ped$ANIMAL)){
      ped <- rbind(data.frame(ANIMAL = ped$MOTHER[which(!ped$MOTHER %in% ped$ANIMAL)],
                              MOTHER = 0, FATHER = 0,
                              stringsAsFactors = F),
                   ped)
    }

    if(any(!ped$FATHER %in% ped$ANIMAL)){
      ped <- rbind(data.frame(ANIMAL = ped$FATHER[which(!ped$FATHER %in% ped$ANIMAL)],
                              MOTHER = 0, FATHER = 0,
                              stringsAsFactors = F),
                   ped)
    }

    ped <- subset(ped, ANIMAL != 0)
    ped <- droplevels(unique(ped))

  }

  if(pedigree.type == "plink"){

    if(ncol(ped) == 5){
      names(ped) <- c("FAMILY", "ANIMAL", "FATHER", "MOTHER", "SEX") #(1=male; 2=female; other=unknown)
      print(paste("Assuming columns ordered as:", paste(names(ped), collapse = " ")))
    }

    if(ncol(ped) == 6){
      names(ped) <- c("FAMILY", "ANIMAL", "FATHER", "MOTHER", "PHENOTYPE") #(1=male; 2=female; other=unknown)
      print(paste("Assuming columns ordered as:", paste(names(ped), collapse = " ")))
    }

    if(!ncol(ped) %in% c(5, 6)){
      stop("Number of columns does not match those expected of PLINK format.")
    }

    # re-code missing values

    for(i in 2:4) ped[which(is.na(ped[,i])),i] <- 0
    for(i in 2:4) ped[which(ped[,i] == -9),i]  <- 0

  }

  ped

  }


