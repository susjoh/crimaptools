#' parse_map: Parse map from CriMAP map output
#' @param mapfile File with output from map
#' @export


parse_map <- function(mapfile){

  x <- readLines(mapfile)

  map1 <- x[grep("Sex_averaged map", x):(grep("Sex-specific map", x)-1)]

  map1 <- map1[3:length(map1)]
  map1 <- map1[-which(map1 == "")]
  map1 <- map1[1:(length(map1)-2)]

  map1 <- data.frame(map1, stringsAsFactors=F)

  map1$map2 <- NA
  for(i in seq(1, nrow(map1), 2)) map1$map2[i] <- map1$map1[i+1]
  map1 <- map1[seq(1, nrow(map1), 2),]

  for(i in 1:1e4) if(length(grep("  ", map1$map1) > 0)) map1$map1  <- gsub("  ", " ", map1$map1) else break
  for(i in 1:1e4) if(length(grep("  ", map1$map2) > 0)) map1$map2  <- gsub("  ", " ", map1$map2) else break


  map <- x[grep("Sex-specific map", x):length(x)]

  map <- map[3:length(map)]
  map <- map[1:(length(map)-7)]
  map <- map[-which(map == "")]

  map <- data.frame(map, stringsAsFactors=F)

  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]

  for(i in 1:1e4) if(length(grep("  ", map$map)  > 0)) map$map   <- gsub("  ", " ", map$map)  else break
  for(i in 1:1e4) if(length(grep("  ", map$map2) > 0)) map$map2  <- gsub("  ", " ", map$map2) else break

  map$Order <- NA
  map$SNP.Name <- NA
  map$cMPosition.Female <- NA
  map$cMPosition.Male <- NA
  map$Female.r <- NA
  map$cMdiff.Female <- NA
  map$Male.r <- NA
  map$cMdiff.Male <- NA
  map$cMPosition <- NA
  map$r <- NA
  map$cMdiff <- NA

  rm(x)

  for(i in 1:nrow(map)){

    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    x1 <- strsplit(map1$map1[i], split = " ")[[1]]
    if(length(-which(x1 == "")) > 0) x1 <- x1[-which(x1 == "")]

    map$Order[i]             <- x[1]
    map$SNP.Name[i]          <- x[2]
    map$cMPosition.Female[i] <- x[3]
    map$cMPosition.Male[i]   <- x[4]
    map$cMPosition[i]        <- x1[3]

    rm(x, x1)

    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]
    x1 <- strsplit(map1$map2[i], split = " ")[[1]]
    if(length(-which(x1 == "")) > 0) x1 <- x1[-which(x1 == "")]

    map$Female.r[i]      <- x[1]
    map$cMdiff.Female[i] <- x[2]
    map$Male.r[i]        <- x[3]
    map$cMdiff.Male[i]   <- x[4]
    map$r[i]             <- x1[1]
    map$cMdiff[i]        <- x1[2]


    rm(x, x1)
  }

  map <- subset(map, select = -c(map, map2))

  analysisID.val <- gsub("\\\\", "/", mapfile)

  analysisID.val <- strsplit(analysisID.val, split = "/")[[1]]
  analysisID.val <- analysisID.val[length(analysisID.val)]
  analysisID.val <- gsub("chr", "", analysisID.val)
  analysisID.val <- gsub(".map", "", analysisID.val, fixed = T)

  map$analysisID <- analysisID.val

  convert.cols <- which(names(map) %in%
                          c("Order", "cMPosition.Female", "cMPosition.Male",
                            "Female.r", "cMdiff.Female", "Male.r", "cMdiff.Male",
                            "cMPosition", "r", "cMdiff"))

  map[,convert.cols] <- sapply(map[,convert.cols], as.numeric)

  map
}
