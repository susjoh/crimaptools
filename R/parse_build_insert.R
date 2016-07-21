#' parse_build_insert: Parse crossover information from CriMAP chrompic output
#' @param buildfile File with output from build with inserted loci
#' @param snpinsert Vector of SNPs that were inserted in prepare and build
#' @import plyr
#' @import stringr
#' @export
#'
#'

# buildfile <- "chr2test.bld1"
# snpinsert <- c("cela1_red_29_4391579", "cela1_red_29_4608019")


parse_built_insert <- function(buildfile, snpinsert){

  buildinfo <- readLines(buildfile)

  #~~ parse the SNP locus index information

  snp.info <- buildinfo[(grep("  0   ", buildinfo)[1]):(grep("ordered.loci", buildinfo)-1)]
  snp.info <- data.frame(Temp = snp.info[which(snp.info != "")])

  snp.info$Index <- sapply(snp.info$Temp, function (x){
    y <- strsplit(as.character(x), split = " ")[[1]]
    y <- y[which(y != "")]
    y[1]
  })

  snp.info$SNP.Name <- sapply(snp.info$Temp, function (x){
    y <- strsplit(as.character(x), split = " ")[[1]]
    y <- y[which(y != "")]
    y[2]
  })

  snp.info <- subset(snp.info, select = -Temp)

  #~~ Find the lines where the SNP is inserted into the map

  insert.lines <- match(snpinsert, buildinfo)

  insert.info <- data.frame(SNP.Name = snpinsert,
                            Start = insert.lines + 1,
                            Stop = c(insert.lines[2:length(insert.lines)], (length(buildinfo)+1))-1)

  #~~ Now determine the best position for each locus

  best.position <- NULL

  for(i in 1:length(snpinsert)){

    insert.ext <- buildinfo[(insert.info$Start[i]):(insert.info$Stop[i])]
    insert.ext <- insert.ext[which(insert.ext != "")]

    snp.ordered <- data.frame(str_locate(insert.ext[1], paste0(" ", snp.info$Index[which(!snp.info$SNP.Name %in% snpinsert)], " ")))
    snp.ordered$Index <- snp.info$Index[which(!snp.info$SNP.Name %in% snpinsert)]

    snp.inserted <- data.frame(Position = str_locate_all(insert.ext[2], "X")[[1]][,1],
                               Li = insert.ext[3:length(insert.ext)])
    snp.inserted$BeforeIndex <- sapply(snp.inserted$Position, function(x){
      snp.ordered$Index[which(snp.ordered$end < x+1)[max(which(snp.ordered$end < x+1))]]
    })

    snp.inserted$InsertedSNP <- snpinsert[i]

    best.position <- rbind(best.position, snp.inserted)
  }

  best.position$Index <- best.position$BeforeIndex

  suppressMessages(best.position <- join(best.position, snp.info))
  best.position <- subset(best.position, select = -Index)

}
