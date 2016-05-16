
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
head(maptab)

maptab2 <- parse_map("crimap/chr3a.map")
head(maptab2)

xovertab <- parse_crossovers("crimap/chr3a.cmp")
head(xovertab)
