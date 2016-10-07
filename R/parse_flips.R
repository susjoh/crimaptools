#' parse_flips: Parse output from CriMAP flips
#' @param flipsfile File with output from flips.
#' @export


parse_flips <- function(flipsfile){

  test2 <- readLines(flipsfile)
  if(length(grep("(= log10_like[orig] - log10_like[curr])", test2, fixed = T)) > 0){

    if(length(grep("LOG_LIKE DECREASED IN", test2) > 0)) test2 <- test2[-grep("LOG_LIKE DECREASED IN", test2)]
    test2 <- test2[-grep("^$", test2)]
    test2 <- test2[(grep("(= log10_like[orig] - log10_like[curr])", test2, fixed = T)+1):length(test2)]


    z <- lapply(test2, function(x){
      y <- strsplit(x, split = "\\s+")[[1]][-1]
      y <- matrix(y, ncol = length(y))
      data.frame(y)
    }
    )

    z <- data.frame(data.table::rbindlist(z))
    names(z) <- paste0("X", 0:(ncol(z)-1))
    names(z)[ncol(z)] <- "LogLi"
    z

  } else {
    message(paste0("No flips present in ", flipsfile, ": has the flips run finished?"))
    return(NULL)
  }
}
