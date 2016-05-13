#' parseFlips: Parse output from CriMAP flips
#' @param flipsfile File with output from flips.


parseFlips <- function(flipsfile){

  test2 <- readLines(flipsfile)
  if(length(grep("LOG_LIKE DECREASED IN", test2) > 0)) test2 <- test2[-grep("LOG_LIKE DECREASED IN", test2)]
  test2 <- test2[-grep("^$", test2)]
  test2 <- test2[(grep("(= log10_like[orig] - log10_like[curr])", test2, fixed = T)+1):length(test2)]


  z <- lapply(test2, function(x){
    y <- strsplit(x, split = "\\s+")[[1]][-1]
    y <- matrix(y, ncol = length(y))
    data.frame(y)
  }
  )

  data.frame(data.table::rbindlist(z))
}
