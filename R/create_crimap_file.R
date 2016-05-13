# Prerequisites (need to be introduced for data checking)
# pedigree must be in format ANIMAL FATHER MOTHER with those headers
# unix tools required on system: fgrep, sed

#' create_crimap_input: Subset a genabel dataset based on specified criteria.
#' @param gwaa.data GenABEL gwaa.data to subset
#' @param chr vector, optional: retain SNPs on particular chromosome(s)
#' @param ped data.frame with pedigree.
#' @param snplist vector of ordered SNP loci.

# gwaa.data <- deer.abel
# familyPedigree <- deer.famped
# mend.errors <- NULL
# analysisID <- "1a"
# snplist <- snpnames(gwaa.data[,chromosome(gwaa.data) ==3])
# is.X = FALSE
# is.Z = FALSE
# pseudoautoSNPs = NULL
# mend.errors = NULL
# outdir = "crimap"
# verbose = TRUE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create Crimap Files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


create_crimap_input <- function(gwaa.data,
                              familyPedigree,
                              analysisID,
                              snplist = NULL,
                              chr = NULL,
                              is.X = FALSE,
                              is.Z = FALSE,
                              pseudoautoSNPs = NULL,
                              mend.errors = NULL,
                              outdir = NULL,
                              verbose = TRUE) {

  #~~ Check that analysisID starts with a number


  if(!substr(analysisID, 1, 1) %in% 0:9) stop("analysisID must start with a numeric character. Please rename it and try again.")
  if(is.X == TRUE & is.Z == TRUE) stop("Cannot specify as both X and Z chromosome")
  if(is.X == TRUE) message("Running as X chromosome: specify pseudoautosomalSNPs to retain some male PAR heterozygotes")
  if(is.Z == TRUE) message("Running as Z chromosome: specify pseudoautosomalSNPs to retain some female PAR heterozygotes")
  if(any(which(!c(familyPedigree[,1], familyPedigree[,2], familyPedigree[,3]) %in% c(idnames(gwaa.data), 0, NA)))){
    stop("familyPedigree must not contain any IDs that are not in the gwaa.data object.")
  }

  #~~ Check family pedigree format

  names(familyPedigree) <- toupper(names(familyPedigree))
  if(!all(names(familyPedigree) %in% c("ID", "ANIMAL", "MUM", "MOM",
                            "MOTHER", "DAM", "DAD", "POP", "FATHER", "SIRE", "FAM", "FAMILY", "IID", "FID"))) stop(simple.ped.name.rules())

  names(familyPedigree)[which(names(familyPedigree) %in% c("ID", "ANIMAL", "IID"))] <- "ANIMAL"
  names(familyPedigree)[which(names(familyPedigree) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"
  names(familyPedigree)[which(names(familyPedigree) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"
  names(familyPedigree)[which(names(familyPedigree) %in% c("FAM", "FID", "FAMILY"))] <- "Family"

  familyPedigree <- familyPedigree[,c("ANIMAL", "FATHER", "MOTHER", "Family")]

  #~~ subset based on snplist

  if(!is.null(snplist)) gwaa.data <- gwaa.data[,snplist]

  nfamilies <- length(unique(familyPedigree$Family))
  nloci <- nsnps(gwaa.data)
  locus.names <- snp.names(gwaa.data)

  #~~ define outfile and write initial information to file

  if(!is.null(outdir))  outfile <- paste0(outdir, "/chr", analysisID, ".gen") else outfile <- paste0("chr", analysisID, ".gen")

  write.table(nfamilies, outfile, row.names = F, quote = F, col.names = F)
  write.table(nloci, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table(locus.names, outfile, row.names = F, quote = F, col.names = F, append=T)
  write.table("", outfile, row.names = F, quote = F, col.names = F, append=T)

  #~~ Subset gwaa.data to include only individuals which are within the pedigree provided

  gwaa.data <- gwaa.data[which(idnames(gwaa.data) %in% c(familyPedigree[,1], familyPedigree[,2], familyPedigree[,3])),]

  #~~ Recode ACGT as 1234 with space delimited between alleles

  temp.geno <- as.character.snp.data(gtdata(gwaa.data))
  temp.geno <- as.matrix(temp.geno)

  if(verbose == TRUE) message("Recoding alleles to numeric values...")

  temp.geno[which(is.na(temp.geno))] <- "0 0"
  temp.geno[which(temp.geno == "A/A")] <- "1 1"
  temp.geno[which(temp.geno == "A/C")] <- "1 2"
  temp.geno[which(temp.geno == "A/G")] <- "1 3"
  temp.geno[which(temp.geno == "A/T")] <- "1 4"
  temp.geno[which(temp.geno == "C/A")] <- "2 1"
  temp.geno[which(temp.geno == "C/C")] <- "2 2"
  temp.geno[which(temp.geno == "C/G")] <- "2 3"
  temp.geno[which(temp.geno == "C/T")] <- "2 4"
  temp.geno[which(temp.geno == "G/A")] <- "3 1"
  temp.geno[which(temp.geno == "G/C")] <- "3 2"
  temp.geno[which(temp.geno == "G/G")] <- "3 3"
  temp.geno[which(temp.geno == "G/T")] <- "3 4"
  temp.geno[which(temp.geno == "T/A")] <- "4 1"
  temp.geno[which(temp.geno == "T/C")] <- "4 2"
  temp.geno[which(temp.geno == "T/G")] <- "4 3"
  temp.geno[which(temp.geno == "T/T")] <- "4 4"
  temp.geno[which(temp.geno == "1/1")] <- "1 1"
  temp.geno[which(temp.geno == "1/2")] <- "1 2"
  temp.geno[which(temp.geno == "2/1")] <- "2 1"
  temp.geno[which(temp.geno == "2/2")] <- "2 2"


  if(verbose == TRUE) message("...done")


  #~~ Create a data frame with ID, FATHER, MOTHER, SEX, and Genotypes.

  temp.geno <- data.frame(temp.geno, stringsAsFactors=F)

  temp.geno$ANIMAL <- idnames(gwaa.data)
  temp.geno$SEX <- phdata(gwaa.data)$sex


  familyPedigree$ANIMAL <- as.character(familyPedigree$ANIMAL)
  temp.geno <- data.table(temp.geno, key = "ANIMAL")
  ped2 <- data.table(unique(familyPedigree[,1:3]), key = "ANIMAL")


  if(verbose == TRUE) message("Merging pedigree and genotype information...")
  temp.geno <- join(ped2, temp.geno, ped2)
  if(verbose == TRUE) message("...done.")

  temp.geno <- data.frame(temp.geno)

  temp.geno <- cbind(temp.geno[,c("ANIMAL", "MOTHER", "FATHER", "SEX")], temp.geno[,which(!names(temp.geno) %in% c("ANIMAL", "MOTHER", "FATHER", "SEX"))])

  #~~ Deal with NA's

  temp.geno$MOTHER[which(is.na(temp.geno$MOTHER))] <- 0
  temp.geno$FATHER[which(is.na(temp.geno$FATHER))] <- 0

  for(i in 1:ncol(temp.geno)) temp.geno[which(is.na(temp.geno[,i])),i] <- "0 0"

  ####################################


  temp.geno <- unique(temp.geno)

  temp.geno <- join(familyPedigree, temp.geno)

  temp.geno <- temp.geno[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", names(temp.geno)[6:ncol(temp.geno)])]


  if(is.X == TRUE){
    xmarkers <- which(!names(temp.geno) %in% pseudoautoSNPs)
    xmarkers <- xmarkers[6:length(xmarkers)]

    heterozygotes <- c("1 2", "1 3", "1 4", "2 1", "2 3", "2 4", "3 1", "3 2", "3 4", "4 1", "4 2", "4 3")

    for(j in xmarkers){

      #~~ remove heterozygotes

      temp.geno[which(temp.geno$SEX == 1 & temp.geno[,j] %in% heterozygotes),j] <- "0 0"

      temp.geno[which(temp.geno$SEX == 1 & temp.geno[,j] == "1 1"),j] <- "1 0"
      temp.geno[which(temp.geno$SEX == 1 & temp.geno[,j] == "2 2"),j] <- "2 0"
      temp.geno[which(temp.geno$SEX == 1 & temp.geno[,j] == "3 3"),j] <- "3 0"
      temp.geno[which(temp.geno$SEX == 1 & temp.geno[,j] == "4 4"),j] <- "4 0"

    }

    rm(xmarkers, heterozygotes)

  }

  if(is.Z == TRUE){
    xmarkers <- which(!names(temp.geno) %in% pseudoautoSNPs)
    xmarkers <- xmarkers[6:length(xmarkers)]

    heterozygotes <- c("1 2", "1 3", "1 4", "2 1", "2 3", "2 4", "3 1", "3 2", "3 4", "4 1", "4 2", "4 3")

    for(j in xmarkers){

      #~~ remove heterozygotes

      temp.geno[which(temp.geno$SEX == 0 & temp.geno[,j] %in% heterozygotes),j] <- "0 0"

      temp.geno[which(temp.geno$SEX == 0 & temp.geno[,j] == "1 1"),j] <- "1 0"
      temp.geno[which(temp.geno$SEX == 0 & temp.geno[,j] == "2 2"),j] <- "2 0"
      temp.geno[which(temp.geno$SEX == 0 & temp.geno[,j] == "3 3"),j] <- "3 0"
      temp.geno[which(temp.geno$SEX == 0 & temp.geno[,j] == "4 4"),j] <- "4 0"

    }

    rm(xmarkers, heterozygotes)

  }


  #~~ deal with mendelian errors

  #   if(!is.null(menderrtab)){
  #     menderrtab <- subset(menderrtab, Chr == chr)
  #     if(nrow(menderrtab) > 0){
  #       for(i in 1:nrow(menderrtab)){
  #         temp.geno[which(temp.geno$ANIMAL == menderrtab$ANIMAL[i]), 6 + as.numeric(menderrtab$Locus[i])] <- "0 0"
  #       }
  #     }
  #   }

  temp.geno <- as.matrix(temp.geno)

  #~~ Create the family and genotype part of the file.

  test <- data.frame(FamSize = tapply(familyPedigree$Family, familyPedigree$Family, length))
  test$Family <- row.names(test)

  familyPedigree <- arrange(familyPedigree, Family)

  all(unique(test$Family) == unique(familyPedigree$Family))

  final.frame <- matrix("", ncol = (ncol(temp.geno)-1), nrow = (sum(test$FamSize) + (nrow(test)*4)))

  counter <- 1

  for(i in unique(familyPedigree$Family)){

    famsize <- test$FamSize[which(test$Family == i)]

    final.frame[counter, 1] <- i
    final.frame[counter+1, 1] <- famsize
    final.frame[(counter+3):(counter+2+famsize),] <- temp.geno[which(temp.geno[,"Family"] == i),-1]

    counter <- counter + 4 + famsize

  }



  final.frame2 <- apply(final.frame, 1, function (x) paste(x, collapse = " "))
  final.frame2 <- gsub("\\s+$", "", final.frame2)

  write.table(final.frame2, outfile, row.names = F, quote = F, col.names = F, append=T)

}

























