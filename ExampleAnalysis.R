library(GenABEL)

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

run_crimap_map("crimap2504.exe", "crimap/chr3a.gen")

maptab <- parse_map_chrompic("crimap/chr3a.cmp")
str(maptab)

maptab2 <- parse_map("crimap/chr3a.map")
str(maptab2)


xovertab <- parse_crossovers("crimap/chr3a.cmp", familyPedigree = deer.famped)
str(xovertab)

physmap <- data.frame(SNP.Name = snpnames(deer.abel)[chromosome(deer.abel) == 3],
                      Position = map(deer.abel)[chromosome(deer.abel) == 3],
                      Order = 1:length(which(chromosome(deer.abel) == 3)),
                      analysisID = "3a")



doub.xover <- check_double_crossovers(xovertab, physmap)
head(doub.xover)

