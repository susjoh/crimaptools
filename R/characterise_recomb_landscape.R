#' characterise_recomb_landscape: Characterise recombination landscape using ordered linkage map and physical map data.
#' @param LG.vec vector of chromosome/linkage map identifiers.
#' @param LG.pos vector of base pair positions.
#' @param LG.cM  vector of linkage map positions.
#' @param window.size integer Window size in which recombination rate is estimated (base pairs).
#' @param verbose logical. Default is TRUE.
#' @export
#

characterise_recomb_landscape <- function(LG.vec, LG.pos, LG.cM, window.size, verbose = TRUE){

  #~~ Make mapdata frame

  mapdata <- data.frame(LG = LG.vec,
                        cM = LG.cM,
                        Position = LG.pos)

  mapdata$Order <- NA
  mapdata$Order[1] <- 1
  for(i in 2:nrow(mapdata)) mapdata$Order[i] <- ifelse(mapdata$LG[i] == mapdata$LG[i-1], mapdata$Order[i-1] + 1, 1)

  mapdata$cMdiff <- c(diff(mapdata$cM), NA)
  mapdata$cMdiff[which(mapdata$cMdiff < 0)] <- NA

  #~~ Create a table of maximum values per chromosome

  max.vals <- data.frame(Est.Length = tapply(LG.pos, LG.vec, max))
  max.vals$LG <- row.names(max.vals)

  #~~ Create a table to save information

  bin.tab <- NULL

  for(i in 1:nrow(max.vals)){
    x <- data.frame(LG = max.vals$LG[i],
                    Window = 1:ceiling(max.vals$Est.Length[i]/window.size))
    x$Start = ((x$Window-1)*window.size) + 1
    x$Stop = x$Start + window.size - 1
    bin.tab <- rbind(bin.tab, x)
    rm(x)
  }

  bin.tab$Locus.Count <- NA
  bin.tab$cM <- NA

  #~~ Extract map information within each bin, and calculate the probability of
  #   crossover between bin boundaries

  for(i in 1:nrow(bin.tab)){

    if(i %in% seq(1, nrow(bin.tab), 100) & verbose == TRUE) print(paste("Running line", i, "of", nrow(bin.tab)))

    # Extract the lines in the window

    x.lines <- which(mapdata$LG == bin.tab$LG[i] &
                       mapdata$Position >= bin.tab$Start[i] &
                       mapdata$Position <= bin.tab$Stop[i])

    bin.tab$Locus.Count[i] <- length(x.lines)

    if(length(x.lines) > 0){

      x <- mapdata[x.lines,]


      # Extract last line before and first line after x.lines (Make NULL of not on same chromosome)

      if(x$Order[1] > 1){

        x.before <- mapdata[min(x.lines) - 1,]

      } else {

        x.before <- NULL
      }

      if(max(x.lines) < nrow(mapdata)) {
        x.after <- mapdata[max(x.lines) + 1,]
        if(x.after$Order[1] == 1) x.after <- NULL

      } else {
        x.after <- NULL
      }


      # Get Recombination length in the window

      recomb.rate <- max(x$cM) - min(x$cM)


      # Get Pr of recombination at the start of the bin

      if(!is.null(x.before)){
        recomb.rate <- recomb.rate + (x.before$cMdiff * (x$Position[1] - bin.tab$Start[i])/(x$Position[1] - x.before$Position))
      }

      if(!is.null(x.after)){
        recomb.rate <- recomb.rate +  (x$cMdiff[nrow(x)] * (bin.tab$Stop[i] - x$Position[nrow(x)])/(x.after$Position - x$Position[nrow(x)]))
      }

      bin.tab$cM[i] <- recomb.rate

      rm(recomb.rate, x, x.before, x.after, x.lines)

    }

  }

  bin.tab

}
