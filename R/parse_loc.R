#' parse_log: Parse output from CriMAP .loc files
#' @param locfile
#' @export


parse_loc <- function(locfile){

  loc.tab <- read.table(locfile, skip = 5)
  names(loc.tab) <- c("CrimapOrder", "SNP.Name", "inf.mei", "inf.mei.PK", "tot_f", "tot_m", "pk_f", "pk_m")
  loc.tab

}
