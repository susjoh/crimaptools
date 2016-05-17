#' InfLengthFunc: Determine positions of first and last informative loci from chrompic output
#' @param CMPstring String of crossover events and informative loci from chrompic output
#' @param position "First" or "Last", return either the first or last informative locus position.
#' @export

InfLengthFunc <- function(CMPstring, position = "First"){
  test <- gsub(" ", "", CMPstring)
  temp <- which(unlist(strsplit(test, split = "")) %in% c("i", "1","0","o", ":", "c"))
  if(position == "First") return(temp[1])
  if(position == "Last"){
    if(length(temp) == 0) return(NA)
    if(length(temp) > 0) return(temp[length(temp)])
  }
}

#' RecCountFunc: Determine recombination count from chrompic output strings
#' @param CMPstring.vec vector of strings of crossover events and informative
#'   loci from chrompic output
#' @export


RecCountFunc <- function(CMPstring.vec){

  unlist(lapply(CMPstring.vec, function(test){

    y <- unlist(strsplit(test, split = ""))
    y <- y[which(y != "-")]
    y[which(y %in% c("o", "c"))] <- 0
    y[which(y %in% c("i", ":"))] <- 1

    y1 <- length(which(c(-9, y) != c(y, -9))) - 2
    if(y1 == -2) y1 <- 0

    return(y1)
  }
  )
  )
}


#' SingletonStatusFunc: Flag if there is a double crossover spanning a single SNP
#' @param CMPstring.vec vector of strings of crossover events and informative
#'   loci from chrompic output
#' @export


SingletonStatusFunc <- function(CMPstring.vec){

  unlist(lapply(CMPstring.vec, function(string){

    y <- unlist(strsplit(string, split = ""))
    y <- y[which(y != "-")]
    y[which(y %in% c("o", "c"))] <- 0
    y[which(y %in% c("i", ":"))] <- 1

    x <- which(c(-9, y) != c(y, -9))
    ifelse(1 %in% diff(x), "Singleton.Present", "Singleton.Absent")
  }
  )
  )
}

#' format_family_pedigree: Format the family pedigree for constistency
#' @param familyPedigree data.frame containing columns ANIMAL, FATHER, MOTHER
#'   and FAMILY. FAMILY defines family groups for crimap. IDs in ANIMAL can be
#'   repeated within and between families if necessary.
#' @export


format_family_pedigree <- function(familyPedigree){

  #~~ Check family pedigree format

  names(familyPedigree) <- toupper(names(familyPedigree))
  if(!all(names(familyPedigree) %in% c("ID", "ANIMAL", "IID",
                                       "MUM", "MOM", "MOTHER", "DAM",
                                       "DAD", "POP","FATHER", "SIRE",
                                       "FAM", "FAMILY", "FID"))) stop(simple.ped.name.rules())

  names(familyPedigree)[which(names(familyPedigree) %in% c("ID", "ANIMAL", "IID"))] <- "ANIMAL"

  names(familyPedigree)[which(names(familyPedigree) %in% c("MUM", "MOM", "MOTHER", "DAM"))] <- "MOTHER"

  names(familyPedigree)[which(names(familyPedigree) %in% c("DAD", "POP", "FATHER", "SIRE"))] <- "FATHER"

  names(familyPedigree)[which(names(familyPedigree) %in% c("FAM", "FID", "FAMILY"))] <- "Family"

  familyPedigree <- familyPedigree[,c("ANIMAL", "FATHER", "MOTHER", "Family")]

  familyPedigree

}

