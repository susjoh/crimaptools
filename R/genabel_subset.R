#' genabel_subset: Subset a genabel dataset based on specified criteria.
#' @param gwaa.data GenABEL gwaa.data to subset
#' @param chr vector, optional: retain SNPs on particular chromosome(s)
#' @param ped data.frame with pedigree.
#' @param snplist vector of ordered SNP loci.


genabel_subset <- function(gwaa.data, chr = NULL, ped = NULL, snplist = NULL){

  if(!is.null(chr)){

    gwaa.data <- gwaa.data[,which(chromosome(gwaa.data) %in% chr)]    #~~ Subset the data by Chromosome

  }

  if(!is.null(ped)){

    id.vec <- as.character(unique(c(pedigree[,1], pedigree[,2], pedigree[,3])))
    gwaa.data <- gwaa.data[which(idnames(gwaa.data) %in% id.vec),]  #~~ subset the data further by IDs in the pedigree object

  }

  #~~ extract snps if list is specified

  if(!is.null(snplist)){

    gwaa.data <- gwaa.data[,snplist]

  }

  #~~ return the new gwaa.data object

  gwaa.data

}
