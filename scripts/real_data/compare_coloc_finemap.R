require(readr)
require(dplyr)

coloc = read_delim("/users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/coloc_results.txt", delim="\t")
finemap = read_delim("/users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_results.txt", delim="\t")

combo = inner_join(finemap, coloc, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))
combo = combo[combo$n_snps.x >= 50 & combo$n_snps.y >= 50,]

cor(combo$clpp, combo$clpp_h4, use="complete.obs")
cor(combo$clpp, combo$clpp_h4, use="complete.obs", method="spearman")
plot(combo$clpp, combo$clpp_h4)
combo$clpp_rank = rank(combo$clpp)
combo$h4_rank = rank(combo$clpp_h4)
plot(combo$clpp_rank, combo$h4_rank, xlab="Order of FINEMAP results", ylab = "Order of COLOC results", pch=18)

# Find points that are very different between the two methods
disagreements = combo[na.omit(abs(combo$clpp_rank - combo$h4_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$clpp_rank - disagreements$h4_rank))),]
points(disagreements$clpp_rank, disagreements$h4_rank, col="red", pch=18)


