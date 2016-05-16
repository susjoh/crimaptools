#' switch_position_cmpstring: function to pull out chromosome segments of
#' different grandparental origin. Used within the function
#' check_double_crossovers
#' @param CMPstring a chrompic strings
#' @param UniqueID unique identifier
#' @export

switch_position_cmpstring <- function(CMPstring, UniqueID = NULL){

  #~~ recode phased alleles as 0 and 1 and remove uninformative loci

  y <- unlist(strsplit(CMPstring, split = ""))
  y[which(y %in% c("o", "c"))] <- 0
  y[which(y %in% c("i", ":"))] <- 1

  y <- data.frame(Phase = y)
  y$Order <- 1:nrow(y)

  y <- y[-which(y$Phase == "-"),]
  y$Phase <- as.numeric(as.character(y$Phase))

  #~~ Determine the first and last positions of each segment between crossover points

  y <- rbind(y, c(2, -999))

  y$Temp  <- c(2, y$Phase[1:(nrow(y)-1)])

  head(y)

  sortvec <- sort(c(which(y$Phase != y$Temp)))

  sortvec <- sort(c(sortvec, sortvec[2:length(sortvec)]-1))

  y2 <- y[sortvec[-length(sortvec)],-ncol(y)]
  y2$Inf.Order <- sortvec[-length(sortvec)]

  #~~ Create a table of the segments

  if(nrow(y2) > 0){
    x <- data.frame(Phase = y2$Phase[seq(1, nrow(y2), 2)],
                    StartPos = y2$Order[seq(1, nrow(y2), 2)],
                    StopPos = y2$Order[seq(2, nrow(y2), 2)],
                    StartInf = y2$Inf.Order[seq(1, nrow(y2), 2)],
                    StopInf = y2$Inf.Order[seq(2, nrow(y2), 2)])
    x$StartSpan <- c(1, x$StopPos[-nrow(x)])
    x$StopSpan  <- c(x$StartPos[-1], x$StopPos[nrow(x)])


    x$InfCount <- x$StopInf - x$StartInf + 1
    x <- subset(x, select = -c(StartInf, StopInf))


    x$Segment <- 1:nrow(x)
    x$Segment.Count <- nrow(x)
    x$Type <- "Mid"
    x$Type[1] <- "First"
    x$Type[nrow(x)] <- "Last"
    if(nrow(x) == 1) x$Type <- "Only"
    x$UniqueID <- UniqueID

  }

  if(nrow(y2) == 0){
    x <- data.frame(Phase = NA,
                    StartPos = NA,
                    StopPos = NA,
                    StartSpan = NA,
                    StopSpan = NA,
                    InfCount = NA,
                    Segment = NA,
                    Segment.Count = NA,
                    Type = NA,
                    UniqueID = UniqueID)
  }

  x

}
