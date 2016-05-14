

prefile <- "crimap/chr1a.pre"
genfile <- "crimap/chr1a.gen"
familyPedigree <- deer.famped
exclusion.threshold <- 0.01

#~~ read in Mendelian errors

x <- readLines(prefile)

x.0 <- x[grep("NONINHERITANCE", x)]
x.1 <- x[grep("NONINHERITANCE", x) + 1]
x.2 <- x[grep("NONINHERITANCE", x) + 2]
x.3 <- x[grep("NONINHERITANCE", x) + 3]

if(length(x.0) > 0){
  menderr <- data.frame(Family     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][3]),
                        ANIMAL     = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][5]),
                        Locus      = sapply(x.0, function (y) strsplit(y, split = " ")[[1]][7]),
                        Mat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][7]),
                        Pat.Allele = sapply(x.1, function (y) strsplit(y, split = " ")[[1]][9]),
                        Mat.UniqueAlleles = x.2,
                        Pat.UniqueAlleles = x.3)

  row.names(menderr) <- 1:nrow(menderr)
  head(menderr)

  menderr$Family     <- gsub(",", "", menderr$Family)
  menderr$ANIMAL     <- gsub(",", "", menderr$ANIMAL)
  menderr$Mat.Allele <- gsub(",", "", menderr$Mat.Allele)
  menderr$Locus      <- gsub(":", "", menderr$Locus)
  menderr$Mat.UniqueAlleles <- gsub("Maternal alleles: ", "", menderr$Mat.UniqueAlleles)
  menderr$Pat.UniqueAlleles <- gsub("Paternal alleles: ", "", menderr$Pat.UniqueAlleles)

  head(menderr)

  rm(x.0, x.1, x.2, x.3)
}

#~~ Read in locus information

nloci <- read.table(genfile, nrows = 1, skip = 1)[1,1]
snplist <- read.table(genfile, skip = 3, nrows = nloci)
names(snplist) <- "SNP.Name"
snplist$Locus <- 0:(nloci - 1)

menderr <- join(menderr, snplist)



#~~ Merge with pedigree information

menderr <- join(menderr, familyPedigree)

#~~ Identify individuals with many mismatches

menderr$Mat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Mat.Allele, menderr$Mat.UniqueAlleles)
menderr$Pat.Mismatch <- mapply(function(x, y) ifelse(length(grep(x, y)) > 0, "no", "yes"), menderr$Pat.Allele, menderr$Pat.UniqueAlleles)

parent.mismatch <- data.frame(table(menderr$ANIMAL, menderr$Mat.Mismatch, menderr$Pat.Mismatch))

parent.mismatch <- unique(parent.mismatch[-which(parent.mismatch$Freq == 0),])

head(parent.mismatch)
names(parent.mismatch) <- c("ANIMAL", "Mat.Mismatches", "Pat.Mismatches", "Count")

test <- rbind(cbind(parent.mismatch[which(parent.mismatch$Mat.Mismatches == "yes"), c("ANIMAL", "Count")],
                    Parent = "MOTHER"),
              cbind(parent.mismatch[which(parent.mismatch$Pat.Mismatches == "yes"), c("ANIMAL", "Count")],
                    Parent = "FATHER"))

head(test)

parent.mismatch <- test

#~~ As there are repeated IDs in there, extract the highest incidence of the count

parent.mismatch$Parent.ID.Code <- paste(parent.mismatch$ANIMAL, parent.mismatch$Parent, sep = "_")
parent.dup.tab <- data.frame(table(parent.mismatch$Parent.ID.Code))
parent.dup.tab <- subset(parent.dup.tab, Freq > 1)

for(i in 1:nrow(parent.dup.tab)){
  count.temp <- parent.mismatch$Count[which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i])]

  if(count.temp[1] != count.temp[2]){
    parent.mismatch <- parent.mismatch[-which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i] &
                                                parent.mismatch$Count == min(count.temp)),]
  }

  if(count.temp[1] == count.temp[2]) {

    parent.mismatch <- parent.mismatch[-which(parent.mismatch$Parent.ID.Code == parent.dup.tab$Var1[i] &
                                                parent.mismatch$Count == min(count.temp))[1],]
  }
}


exclusion.threshold <- nsnps(abeldata) * 0.001

ggplot(parent.mismatch[which(parent.mismatch$Count < exclusion.threshold),], aes(Count)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~Parent)

remove.ids <- parent.mismatch[which(parent.mismatch$Count > exclusion.threshold),]

remove.ids$ANIMAL <- as.numeric(as.character(remove.ids$ANIMAL))
remove.ids$Parent <- as.character(remove.ids$Parent)


#~~ revise the pedigree and error table

pedigree.new <- pedigree

if(nrow(remove.ids) > 0){
  for(i in 1:nrow(remove.ids)){

    pedigree.new[which(pedigree.new$ANIMAL == remove.ids$ANIMAL[i]), remove.ids$Parent[i]] <- 0

    if(remove.ids$Parent[i] == "MOTHER") menderr <- menderr[-which(menderr$ANIMAL == remove.ids$ANIMAL[i] & menderr$Mat.Mismatch == "yes"),]
    if(remove.ids$Parent[i] == "FATHER") menderr <- menderr[-which(menderr$ANIMAL == remove.ids$ANIMAL[i] & menderr$Pat.Mismatch == "yes"),]

  }
}

#~~ Deal with problematic loci within parents by removing genotype in parents and offspring.

head(menderr)

newerrtab <- menderr[,c("ANIMAL", "Locus", "Chr")]

mattab <- menderr[which(menderr$Mat.Mismatch == "yes"),c("MOTHER", "Locus", "Chr")]
names(mattab) <- c("ANIMAL", "Locus", "Chr")
pattab <- menderr[which(menderr$Pat.Mismatch == "yes"),c("FATHER", "Locus", "Chr")]
names(pattab) <- c("ANIMAL", "Locus", "Chr")

newerrtab <- rbind(newerrtab, mattab, pattab)

tapply(newerrtab$Locus, newerrtab$Chr, max)

head(newerrtab)

