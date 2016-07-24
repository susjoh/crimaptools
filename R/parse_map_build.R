#' parse_map_build: Parse map from CriMAP chrompic output
#' @param buildfile File with output from build
#' @export


parse_map_build <- function(buildfile){

  x <- readLines(buildfile)
  map <- x[grep("Sex-specific map", x):length(x)]


  map <- map[3:length(map)]
  map <- map[1:(grep("denotes recomb. frac. held fixed in this analysis", map) - 1)]
  map <- map[-which(map == "")]

  map <- data.frame(map, stringsAsFactors=F)

  map$map2 <- NA
  for(i in seq(1, nrow(map), 2)) map$map2[i] <- map$map[i+1]
  map <- map[seq(1, nrow(map), 2),]

  for(i in 1:1e4) if(length(grep("  ", map$map)  > 0)) map$map   <- gsub("  ", " ", map$map)  else break
  for(i in 1:1e4) if(length(grep("  ", map$map2) > 0)) map$map2  <- gsub("  ", " ", map$map2) else break

  rm(x)

  map$Order <- NA
  map$SNP.Name <- NA
  map$cMPosition.Female <- NA
  map$cMPosition.Male <- NA
  map$Female.r <- NA
  map$cMdiff.Female <- NA
  map$Male.r <- NA
  map$cMdiff.Male <- NA

  for(i in 1:nrow(map)){

    x <- strsplit(map$map[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Order[i]             <- x[1]
    map$SNP.Name[i]          <- x[2]
    map$cMPosition.Female[i] <- x[3]
    map$cMPosition.Male[i]   <- x[4]

    rm(x)

    x <- strsplit(map$map2[i], split = " ")[[1]]
    if(length(-which(x == "")) > 0) x <- x[-which(x == "")]

    map$Female.r[i]      <- x[1]
    map$cMdiff.Female[i] <- x[2]
    map$Male.r[i]        <- x[3]
    map$cMdiff.Male[i]   <- x[4]


    rm(x)
  }

  map <- subset(map, select = -c(map, map2))

  analysisID.val <- gsub("\\\\", "/", buildfile)

  analysisID.val <- strsplit(analysisID.val, split = "/")[[1]]
  analysisID.val <- analysisID.val[length(analysisID.val)]
  analysisID.val <- gsub("chr", "", analysisID.val)
  analysisID.val <- gsub(".map", "", analysisID.val, fixed = T)

  map$analysisID <- analysisID.val

  convert.cols <- which(names(map) %in%
                          c("Order", "cMPosition.Female", "cMPosition.Male",
                            "Female.r", "cMdiff.Female", "Male.r", "cMdiff.Male"))

  map[,convert.cols] <- sapply(map[,convert.cols], as.numeric)

  map


}
