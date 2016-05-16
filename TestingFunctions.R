
devtools::document()
data(deer)



create_crimap_input (deer.abel, deer.famped, analysisID = "3a",
                     chr = 3, outdir = "crimap", clear.existing.analysisID = TRUE)

run_crimap_prepare("crimap2504.exe", "crimap/chr3a.gen")

parse_mend_err("crimap/chr3a.pre", "crimap/chr3a.gen", familyPedigree = deer.famped)

create_crimap_input (deer.abel, deer.famped, analysisID = "3a",
                     chr = 3, outdir = "crimap", clear.existing.analysisID = TRUE,
                     use.mnd = TRUE)

run_crimap_prepare("crimap2504.exe", "crimap/chr3a.gen")

parse_mend_err("crimap/chr3a.pre", "crimap/chr3a.gen", familyPedigree = deer.famped)

run_crimap_chrompic("crimap2504.exe", "crimap/chr3a.gen")


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
