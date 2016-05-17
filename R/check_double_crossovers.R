#' check_double_crossovers: Split chromosomes into segments of runs of SNPs
#' inherited from particular grandparents. Used as quality control to determine
#' double crossovers over short genomic distances.
#' @param parsed.xovers Output from parse_crossovers
#' @param physical.map Optional data.frame. Must have columns as follows, or
#'   will throw error: c("SNP.Name", "Position", "Order", "analysisID").
#' @import plyr
#' @import data.table
#' @export


#
# parsed.xovers <- xovertab
# physical.map <- NULL

check_double_crossovers <- function(parsed.xovers, physical.map = NULL){

  if(!is.null(physical.map) & all(names(physical.map) != c("SNP.Name", "Position", "Order", "analysisID")))
  {
    stop("physical.map should be in the format c(\"SNP.Name\", \"Position\", \"Order\", \"analysisID\").")
  }


  head(parsed.xovers[,-which(names(parsed.xovers) %in% c("data"))])

  #~~ create a smaller table to determine the recombination break points

  distab <- subset(parsed.xovers, select = c(Family, RRID, parent, UniqueID, data, analysisID))

  #~~ Apply function using lapply and create a table with a unique line for each segment

  message("Splitting chromosome into segments of shared grandparental origin")
  test <- mapply(switch_position_cmpstring, distab$data, distab$UniqueID, SIMPLIFY = F)

  switchtab <- data.frame(data.table::rbindlist(test))

  #~~ Merge with distab and add map information

  suppressMessages(switchtab <- join(switchtab, distab))
  switchtab <- subset(switchtab, select = -data)

  #~~ Merge with map information

  if(!is.null(physical.map)){

    suppressMessages({
      names(physical.map)[2:3] <- c("StartPos.GenomePos", "StartPos")
      switchtab <- join(switchtab, physical.map[,2:4])

      names(physical.map)[2:3] <- c("StopPos.GenomePos", "StopPos")
      switchtab <- join(switchtab, physical.map[,2:4])

      names(physical.map)[2:3] <- c("StartSpan.GenomePos", "StartSpan")
      switchtab <- join(switchtab, physical.map[,2:4])


      names(physical.map)[2:3] <- c("StopSpan.GenomePos", "StopSpan")
      switchtab <- join(switchtab, physical.map[,2:4])
    })

    switchtab$PosLength  <- switchtab$StopPos.GenomePos - switchtab$StartPos.GenomePos
    switchtab$SpanLength <- switchtab$StopSpan.GenomePos - switchtab$StartSpan.GenomePos
  } else {

#     switchtab$StartPos.GenomePos <- NA
#     switchtab$StopPos.GenomePos <- NA
#     switchtab$StartSpan.GenomePos <- NA
#     switchtab$StopSpan.GenomePos <- NA
#     switchtab$PosLength <- NA
#     switchtab$SpanLength <- NA

  }

  switchtab$Singleton <- ifelse(switchtab$InfCount == 1, "yes", "no")

  switchtab

}


