
devtools::document()
data(deer)

create_crimap_input (deer.abel, deer.famped, analysisID = "1a", snplist = snpnames(deer.abel)[1:3], outdir = "crimap")


library(crimaptools)

nids(deer.abel)
nsnps(deer.abel)
table(chromosome(deer.abel))

sub1 <- genabelSubset(gwaa.data = deer.abel, ped = deer.ped, chr = 1)
nids(sub1)
nsnps(sub1)
table(chromosome(sub1))

chr.vec <- snp.names(sub1)[1:3]
chr.vec
sub1 <- genabelSubset(gwaa.data = deer.abel, ped = deer.ped, snplist = chr.vec)
