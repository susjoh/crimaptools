#' parse_map_chrompic: Parse map from CriMAP chrompic output
#' @param chrompicfile File with output from chrompic
#' @export


parse_map_chrompic <- function(chrompicfile){

  x <- readLines(chrompicfile)
  map <- x[grep("Sex_averaged", x):length(x)]
  map <- map[3:length(map)]
  map <- map[1:(length(map)-6)]

  map <- data.frame(map, stringsAsFactors=F)

  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]

  for(i in 1:10) map$map  <- gsub("  ", " ", map$map)
  for(i in 1:10) map$map2 <- gsub("  ", " ", map$map2)

  rm(x)

  map$Order <- NA
  map$SNP.Name <- NA
  map$cMPosition <- NA
  map$r <- NA
  map$cMdiff <- NA

  for(i in 1:nrow(map)){
    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Order[i] <- x[1]
    map$SNP.Name[i] <- x[2]
    map$cMPosition[i] <- x[3]
    rm(x)

    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$r[i] <- x[1]
    map$cMdiff[i] <- x[2]

    rm(x)
  }

  map <- subset(map, select = -c(map, map2))

  analysisID.val <- gsub("\\\\", "/", chrompicfile)

  analysisID.val <- strsplit(analysisID.val, split = "/")[[1]]
  analysisID.val <- analysisID.val[length(analysisID.val)]
  analysisID.val <- gsub("chr", "", analysisID.val)
  analysisID.val <- gsub(".cmp", "", analysisID.val, fixed = T)

  map$analysisID <- analysisID.val

  convert.cols <- which(names(map) %in%
                          c("Order", "cMPosition", "r", "cMdiff"))

  map[,convert.cols] <- sapply(map[,convert.cols], as.numeric)


  map


}
