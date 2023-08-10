chap1_ccg = read.delim("ccg_in_GEPs_genotypes.txt", header = F)[,1]

chap2_ccg = read.delim("Chapter_3/ccg_in_GEPs_genotypes.txt", header = F)[,1]

length(intersect(chap1_ccg, chap2_ccg))

common_ccg_chap1_chap2 = intersect(chap1_ccg, chap2_ccg)
writeLines(common_ccg_chap1_chap2, "common_ccg_chap1_chap2.txt")
