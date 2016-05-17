#' revise_double_crossovers: remove erroneously phased loci from parsed
#' crossover data.
#' @param parsed.xovers output from parse_crossovers
#' @param removesections data.frame that contains the following fields: StartPos, StopPos and
#'   UniqueID from check_double_crossovers function. Frame should only include
#'   values that should be removed
#' @export

revise_double_crossovers <- function(parsed.xovers, removesections){

  for(i in 1:nrow(removesections)){

    if(i %in% seq(1, nrow(removesections), 100)) print(paste("Fixing Problem", i, "of", nrow(removesections)))

    x <- parsed.xovers$data[which(parsed.xovers$UniqueID == removesections$UniqueID[i])]
    x <- unlist(strsplit(x, split = ""))
    x[removesections$StartPos[i]:removesections$StopPos[i]] <- "-"
    x <- paste(x, collapse = "")

    parsed.xovers$data[parsed.xovers$UniqueID == removesections$UniqueID[i]] <- x

  }

  parsed.xovers$RecombCount <- RecCountFunc(parsed.xovers$data)

  parsed.xovers
}
