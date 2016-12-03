#' parse_mend_err: Parse output from CriMAP .loc files
#' @param prefile .pre file (output from prepare function)
#' @param genfile .gen file
#' @param save.mendfile logical. Default = TRUE, saves .mnd file
#' @param familyPedigree data.frame containing columns ANIMAL, FATHER, MOTHER
#'   and FAMILY. FAMILY defines family groups for crimap. IDs in ANIMAL can be
#'   repeated within and between families if necessary.
#' @param is.X logical. If true, then Mendelian errors are dealt with
#'   differently for the homogametic sex (see below).
#' @param is.Z logical. If true, then Mendelian errors are dealt with
#'   differently for the homogametic sex (see below).
#' @param pseudoautoSNPs Character vector of pseudoautosomal SNP IDs. Only used
#'   if is.X or is.Z == TRUE.
#' @param genabel.phdata GenABEL phdata object, or data frame with fields
#'   information on "id" and sex"
#' @import plyr
#' @export

# prefile <- "crimap/crimap_b/chr30b.pre"
# genfile <- "crimap/crimap_b/chr30b.gen"
# save.mendfile <- TRUE
# familyPedigree <- famped
# is.X = TRUE
# pseudoautoSNPs = pseudoautoSNPs
# genabel.phdata = phdata(abeldata)

parse_mend_err <- function(prefile,
                           genfile,
                           save.mendfile = TRUE,
                           familyPedigree,
                           is.X = NULL,
                           is.Z = NULL,
                           pseudoautoSNPs = NULL,
                           genabel.phdata = NULL){

  if(is.X == TRUE) message("Running as X chromosome: analysis may take time. Specify pseudoautosomalSNPs to retain some male PAR heterozygotes")

  if(is.Z == TRUE) message("Running as Z chromosome: analysis may take time. Specify pseudoautosomalSNPs to retain some female PAR heterozygotes")

  if(is.X == TRUE || is.Z == TRUE){
    if(is.null(pseudoautoSNPs)){
      message("No pseudoautosomal SNPs specified")
    }
    if(is.null(phdata)) stop("Running Mendelian check requires information on offspring sex. Please include
                             genabel.phdata object")
  }

  #~~ read in Mendelian errors

  x <- readLines(prefile)

  x.0 <- x[grep("NONINHERITANCE", x)]
  x.1 <- x[grep("NONINHERITANCE", x) + 1]
  x.2 <- x[grep("NONINHERITANCE", x) + 2]
  x.3 <- x[grep("NONINHERITANCE", x) + 3]

  if(length(x.0) == 0){

    if(file.exists(gsub(".gen", ".mnd", genfile, fixed = T))) {
      message("No Mendelian errors detected. No changes made to .mnd file.")
    }

    if(!file.exists(gsub(".gen", ".mnd", genfile, fixed = T))) {
      message("No Mendelian errors detected. Blank .mnd file created.")
      writeLines("ANIMAL\tSNP.Name", gsub(".gen", ".mnd", genfile, fixed = T))
    }

  } else {

    menderr <- data.frame(Family     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][3]),
                          ANIMAL     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][5]),
                          Locus      = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][7]),
                          Mat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][7]),
                          Pat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][9]),
                          Mat.UniqueAlleles = x.2,
                          Pat.UniqueAlleles = x.3)

    row.names(menderr) <- 1:nrow(menderr)
    head(menderr)

    menderr$Family     <- gsub(",", "", menderr$Family)
    menderr$ANIMAL     <- gsub(",", "", menderr$ANIMAL)
    menderr$Mat.Allele <- gsub(",", "", menderr$Mat.Allele)
    menderr$Locus      <- gsub(":", "", menderr$Locus)
    menderr$Mat.UniqueAlleles <- gsub("Maternal alleles: ", "", menderr$Mat.UniqueAlleles)
    menderr$Pat.UniqueAlleles <- gsub("Paternal alleles: ", "", menderr$Pat.UniqueAlleles)

    head(menderr)

    rm(x.0, x.1, x.2, x.3)

    #~~ Read in locus information

    nloci <- read.table(genfile, nrows = 1, skip = 1)[1,1]
    snplist <- read.table(genfile, skip = 3, nrows = nloci)
    names(snplist) <- "SNP.Name"
    snplist$Locus <- 0:(nloci - 1)

    suppressMessages(menderr <- join(menderr, snplist))

    #~~ Merge with pedigree information

    suppressMessages(menderr <- join(menderr, familyPedigree))

    #~~ Find the parental mismatch

    menderr$Mat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Mat.Allele, menderr$Mat.UniqueAlleles)
    menderr$Pat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Pat.Allele, menderr$Pat.UniqueAlleles)

    #~~ If sex linked, deal with incorrect mismatches. Sex 1 = male, 0 = female


    if(is.X == TRUE || is.Z == TRUE){

      names(genabel.phdata)[which(names(genabel.phdata) == "id")] <- "ANIMAL"
      menderr <- join(menderr, genabel.phdata[,c("ANIMAL", "sex")])

    }

    # for X, if a father's allele is not in a daughter, it is an error.
    #        if a father's allele is not in a son, it is not an error.

    if(is.X == TRUE){

      table(menderr$sex, menderr$Pat.Mismatch)

      menderr <- menderr[-which(menderr$Pat.Mismatch == "yes" & menderr$sex == 1),]


    }

    if(is.Z == TRUE){

      table(menderr$sex, menderr$Pat.Mismatch)

      menderr <- menderr[-which(menderr$Mat.Mismatch == "yes" & menderr$sex == 0),]


    }



    newerrtab <- menderr[,c("ANIMAL", "SNP.Name")]

    mattab <- menderr[which(menderr$Mat.Mismatch == "yes"),c("MOTHER", "SNP.Name")]
    names(mattab) <- c("ANIMAL", "SNP.Name")
    pattab <- menderr[which(menderr$Pat.Mismatch == "yes"),c("FATHER", "SNP.Name")]
    names(pattab) <- c("ANIMAL", "SNP.Name")

    newerrtab <- rbind(newerrtab, mattab, pattab)
    newerrtab <- unique(newerrtab)

    newerrtab

    if(save.mendfile == TRUE){
      message(paste0("Writing new file ", gsub(".gen", ".mnd", genfile, fixed = T)))
      write.table(newerrtab,
                  gsub(".gen", ".mnd", genfile, fixed = T),
                  row.names = F, sep = "\t", quote = F)

      write.table(menderr,
                  gsub(".gen", ".mndverbose", genfile, fixed = T),
                  row.names = F, sep = "\t", quote = F)
    }

  }
}
