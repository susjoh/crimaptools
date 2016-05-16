#' create_crimap_input: Create a CriMAP input file.
#' @param gwaa.data GenABEL gwaa.data object containing genomic data
#' @param familyPedigree data.frame containing columns ANIMAL, FATHER, MOTHER
#'   and FAMILY. FAMILY defines family groups for crimap. IDs in ANIMAL can be
#'   repeated within and between families if necessary.
#' @param analysisID string. Used as an identifier for the chromosome. Must
#'   begin with a numeric character as per crimap convention; for example,
#'   analysisID = "1a" will output the file "chr1a.gen".
#' @param snplist vector, optional. A list of ordered SNPs. UNless build is run,
#'   SNPs will be assumed to be in this order.
#' @param is.X logical. If true, then all heterozygotes are scored as missing in
#'   males, except for pseudoautosomal SNPs when defined (see below)
#' @param is.Z logical. If true, then all heterozygotes are scored as missing in
#'   females, except for pseudoautosomal SNPs when defined (see below)
#' @param pseudoautoSNPs Character vector of pseudoautosomal SNP IDs. Only used
#'   if is.X or is.Z == TRUE.
#' @param use.mnd logical, default = FALSE. Masks parent offspring genotype mismatches based on output from parse_mend_err
#' @param outdir String, optional. Specify the path of the directory in which
#'   the output file should be written.
#' @param verbose Logical. FALSE will suppress messages.
#' @param clear.existing.analysisID Logical, default is TRUE. Overwrites .gen
#'   file and deletes crimap files that have the analysisID in the target
#'   directory.
#' @import GenABEL
#' @import data.table
#' @import plyr
#' @export
#
#
# gwaa.data <- deer.abel
# familyPedigree <- deer.famped
# mend.errors <- NULL
# analysisID <- "1a"
# snplist <- snpnames(gwaa.data[,chromosome(gwaa.data) ==1])
# chr = NULL
# is.X = FALSE
# is.Z = FALSE
# pseudoautoSNPs = NULL
# use.mnd <- TRUE
# outdir = "crimap"
# verbose = TRUE
# clear.existing.analysisID = TRUE

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
                                use.mnd = FALSE,
                                outdir = NULL,
                                verbose = TRUE,
                                clear.existing.analysisID = TRUE) {

  #~~ Check that analysisID starts with a number


  if(!substr(analysisID, 1, 1) %in% 0:9) stop("analysisID must start with a numeric character. Please rename it and try again.")

  if(is.X == TRUE & is.Z == TRUE) stop("Cannot specify as both X and Z chromosome")

  if(!is.null(chr) & !is.null(snplist)) stop("Please only specify only SNP list or chromosome")

  if(is.X == TRUE) message("Running as X chromosome: specify pseudoautosomalSNPs to retain some male PAR heterozygotes")

  if(is.Z == TRUE) message("Running as Z chromosome: specify pseudoautosomalSNPs to retain some female PAR heterozygotes")

  if(any(which(!c(familyPedigree[,1], familyPedigree[,2], familyPedigree[,3]) %in% c(idnames(gwaa.data), 0, NA)))){

    stop("familyPedigree must not contain any IDs that are not in the gwaa.data object.")

  }

  if(use.mnd == TRUE & !paste0("chr", analysisID, ".mnd") %in% dir(outdir)){
    stop(paste0("use.mnd == TRUE, but there is no file named chr",
                analysisID, ".mnd. Change to FALSE and/or run parse_mend_err()."))
  }

  #~~ Delete existing files if specified

  if(!is.null(outdir)) out.path.stem <- paste0(outdir, "/chr", analysisID)
  if( is.null(outdir)) out.path.stem <- paste0("chr", analysisID)


  if(clear.existing.analysisID == TRUE){

    del.vec <- grep(paste0("chr", analysisID, "."), dir(outdir), value = T)

    if(use.mnd == TRUE & paste0("chr", analysisID, ".mnd") %in% del.vec){

        del.vec <- del.vec[-which(del.vec == paste0("chr", analysisID, ".mnd"))]

        }


    if(Sys.info()["sysname"] == "Windows") {

      if(length(del.vec) > 0){

        if(!is.null(outdir)) del.vec <- paste0(outdir, "\\", del.vec)

        for(i in del.vec){

          system("cmd", input = paste0("del ", i), show.output.on.console = F)
        }

      }

    } else {

      if(!is.null(outdir)) del.vec <- paste0(outdir, "/", del.vec)


      for(i in del.vec){

        system(paste0("rm ", i), show.output.on.console = F)
      }

    }

  }


  #~~ Check family pedigree format

  familyPedigree <- format_family_pedigree(familyPedigree)

  #~~ subset based on snplist

  if(!is.null(snplist)) gwaa.data <- gwaa.data[,snplist]
  if(!is.null(chr))     gwaa.data <- gwaa.data[,which(chromosome(gwaa.data) == chr)]


  nfamilies <- length(unique(familyPedigree$Family))

  nloci <- nsnps(gwaa.data)

  locus.names <- snp.names(gwaa.data)

  #~~ define outfile and write initial information to file


  outfile <- paste0(out.path.stem, ".gen")

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

  suppressMessages(temp.geno <- join(ped2, temp.geno))

  if(verbose == TRUE) message("...done.")

  temp.geno <- data.frame(temp.geno)

  temp.geno <- cbind(temp.geno[,c("ANIMAL", "MOTHER", "FATHER", "SEX")], temp.geno[,which(!names(temp.geno) %in% c("ANIMAL", "MOTHER", "FATHER", "SEX"))])

  #~~ Deal with NA's

  temp.geno$MOTHER[which(is.na(temp.geno$MOTHER))] <- 0

  temp.geno$FATHER[which(is.na(temp.geno$FATHER))] <- 0

  for(i in 1:ncol(temp.geno)) temp.geno[which(is.na(temp.geno[,i])),i] <- "0 0"

  #~~ Format temp.geno

  temp.geno <- unique(temp.geno)

  suppressMessages(temp.geno <- join(familyPedigree, temp.geno))

  temp.geno <- temp.geno[,c("Family", "ANIMAL", "MOTHER", "FATHER", "SEX", names(temp.geno)[6:ncol(temp.geno)])]

  #~~ Deal with sex chromosomes

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

  if(verbose == TRUE) message(paste0("Parsing and writing to ", outfile, "..."))


  #~~ deal with mendelian errors if specified

  if(use.mnd == TRUE){

    menderrtab <- read.table(paste0(out.path.stem, ".mnd"), header = T, stringsAsFactors = F)

    if(nrow(menderrtab) > 0){


      for(i in 1:nrow(menderrtab)){

        temp.geno[which(temp.geno$ANIMAL == menderrtab$ANIMAL[i]),
                  which(names(temp.geno) == menderrtab$SNP.Name[i])] <- "0 0"

      }

    }

  }

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

  #~~ remove white space

  final.frame2 <- apply(final.frame, 1, function (x) paste(x, collapse = " "))
  final.frame2 <- gsub("\\s+$", "", final.frame2)

  #~~ Write to file

  write.table(final.frame2, outfile, row.names = F, quote = F, col.names = F, append=T)

  if(verbose == TRUE) message("...done")

}

























