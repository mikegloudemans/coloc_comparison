require(readr)
require(dplyr)
require(UpSetR)

# Load all colocalization results for different methods
coloc = read_delim("/users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/coloc_results.txt", delim="\t")
finemap = read_delim("/users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_results.txt", delim="\t")
finemap2 = read_delim("/users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_c2_results.txt", delim="\t")
baseline = read_delim("/users/mgloud/projects/coloc_comparisons/output/baseline/baseline_results.txt", delim="\t")
smr = read_delim("/users/mgloud/projects/coloc_comparisons/output/smr/smr_results.txt", delim="\t")
gsmr = read_delim("/users/mgloud/projects/coloc_comparisons/output/gsmr/gsmr_results.txt", delim="\t")

names(baseline)[6] = "base_gwas_file"

idx = which(!(names(baseline) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(baseline)[idx] = paste0(names(baseline[idx]), "_baseline")
idx = which(!(names(finemap) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(finemap)[idx] = paste0(names(finemap[idx]), "_finemap_c1")
idx = which(!(names(finemap2) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(finemap2)[idx] = paste0(names(finemap2[idx]), "_finemap_c2")
idx = which(!(names(coloc) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(coloc)[idx] = paste0(names(coloc[idx]), "_coloc")
idx = which(!(names(smr) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(smr)[idx] = paste0(names(smr[idx]), "_smr")
idx = which(!(names(gsmr) %in% c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file")))
names(gsmr)[idx] = paste0(names(gsmr[idx]), "_gsmr")

# Join all results for the same locus into a single row with scores for each method
combo = full_join(baseline, coloc, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))
combo = full_join(combo, finemap2, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))
combo = full_join(combo, finemap, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))
combo = full_join(combo, smr, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))
combo = full_join(combo, gsmr, by=c("ref_snp", "eqtl_file", "gwas_trait", "feature", "base_gwas_file"))


#######################################################
# Compare COLOC and FINEMAP with one causal variant max
#######################################################
subcombo = combo[combo$n_snps_finemap_c1 >= 50 & combo$n_snps_coloc >= 50,]

cor(subcombo$clpp_finemap_c1, subcombo$clpp_h4_coloc, use="complete.obs")
cor(subcombo$clpp_finemap_c1, subcombo$clpp_h4_coloc, use="complete.obs", method="spearman")
plot(subcombo$clpp_finemap_c1, subcombo$clpp_h4_coloc)
subcombo$clpp_finemap_c1_rank = rank(subcombo$clpp_finemap_c1)
subcombo$h4_rank = rank(subcombo$clpp_h4_coloc)
plot(subcombo$clpp_finemap_c1_rank, subcombo$h4_rank, xlab="Order of FINEMAP results", ylab = "Order of COLOC results", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$clpp_finemap_c1_rank - subcombo$h4_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$clpp_finemap_c1_rank - disagreements$h4_rank))),]
points(disagreements$clpp_finemap_c1_rank, disagreements$h4_rank, col="red", pch=18)

#######################################################
# Compare FINEMAP with 1 vs. 2 causal variants max
#######################################################
subcombo = combo[combo$n_snps_finemap_c1 >= 50 & combo$n_snps_finemap_c2 >= 50,]
cor(subcombo$clpp_finemap_c1, subcombo$clpp_finemap_c2, use="complete.obs")
cor(subcombo$clpp_finemap_c1, subcombo$clpp_finemap_c2, use="complete.obs", method="spearman")
plot(subcombo$clpp_finemap_c1, subcombo$clpp_finemap_c2)
subcombo$clpp_finemap_c1_rank = rank(subcombo$clpp_finemap_c1)
subcombo$clpp_finemap_c2_rank = rank(subcombo$clpp_finemap_c2)
plot(subcombo$clpp_finemap_c1_rank, subcombo$clpp_finemap_c2_rank, xlab="FINEMAP (max 1 causal variant)", ylab = "FINEMAP (max 2 causal variants)", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$clpp_finemap_c1_rank - subcombo$clpp_finemap_c2_rank) > 1000),]
# Sort them from most to least disagreeing
signed_disagreements = disagreements[rev(order(disagreements$clpp_finemap_c1_rank - disagreements$clpp_finemap_c2_rank)),]
disagreements = disagreements[rev(order(abs(disagreements$clpp_finemap_c1_rank - disagreements$clpp_finemap_c2_rank))),]
points(disagreements$clpp_finemap_c1_rank, disagreements$clpp_finemap_c2_rank, col="red", pch=18)


############################################
# Compare FINEMAP w/ baseline
############################################

# Baseline 1

subcombo = combo[combo$n_snps_finemap_c1 >= 50 & combo$n_snps_baseline >= 50,]

cor(subcombo$clpp_finemap_c1, subcombo$baseline_pval_baseline, use="complete.obs")
cor(subcombo$clpp_finemap_c1, subcombo$baseline_pval_baseline, use="complete.obs", method="spearman")
plot(subcombo$clpp_finemap_c1, subcombo$baseline_pval_baseline)
subcombo$clpp_finemap_c1_rank = rank(subcombo$clpp_finemap_c1)
subcombo$baseline_rank = rank(subcombo$baseline_pval_baseline)
plot(subcombo$clpp_finemap_c1_rank, subcombo$baseline_rank, xlab="Order of FINEMAP results", ylab = "Order of baseline results", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$clpp_finemap_c1_rank - subcombo$baseline_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$clpp_finemap_c1_rank - disagreements$baseline_rank))),]
points(disagreements$clpp_finemap_c1_rank, disagreements$baseline_rank, col="red", pch=18)


# Baseline 2 (right now it turns out to just be the same thing...probably b/c of selection criteria
subcombo = combo[combo$n_snps_finemap_c1 >= 50 & combo$n_snps_baseline >= 50 & combo$baseline_pval_baseline > 1,]

cor(subcombo$clpp_finemap_c1, subcombo$baseline_pval2_baseline, use="complete.obs")
cor(subcombo$clpp_finemap_c1, subcombo$baseline_pval2_baseline, use="complete.obs", method="spearman")
plot(subcombo$clpp_finemap_c1, subcombo$baseline_pval2_baseline)
subcombo$clpp_finemap_c1_rank = rank(subcombo$clpp_finemap_c1)
subcombo$baseline2_rank = rank(subcombo$baseline_pval2_baseline)
plot(subcombo$clpp_finemap_c1_rank, subcombo$baseline2_rank, xlab="Order of FINEMAP results", ylab = "Order of baseline results", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$clpp_finemap_c1_rank - subcombo$baseline2_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$clpp_finemap_c1_rank - disagreements$baseline2_rank))),]
points(disagreements$clpp_finemap_c1_rank, disagreements$baseline2_rank, col="red", pch=18)

# Baseline 2 and FINEMAP-c2

subcombo = combo[combo$n_snps_finemap_c2 >= 50 & combo$n_snps_baseline >= 50,]

cor(subcombo$clpp_finemap_c2, subcombo$baseline_pval_baseline, use="complete.obs")
cor(subcombo$clpp_finemap_c2, subcombo$baseline_pval_baseline, use="complete.obs", method="spearman")
plot(subcombo$clpp_finemap_c2, subcombo$baseline_pval_baseline)
subcombo$clpp_finemap_c2_rank = rank(subcombo$clpp_finemap_c2)
subcombo$baseline_rank = rank(subcombo$baseline_pval_baseline)
plot(subcombo$clpp_finemap_c2_rank, subcombo$baseline_rank, xlab="Order of FINEMAP-c2 results", ylab = "Order of baseline results", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$clpp_finemap_c2_rank - subcombo$baseline_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$clpp_finemap_c2_rank - disagreements$baseline_rank))),]
points(disagreements$clpp_finemap_c2_rank, disagreements$baseline_rank, col="red", pch=18)



############################################
# Compare COLOC w/ baseline
############################################

subcombo = combo[combo$n_snps_coloc >= 50 & combo$n_snps_baseline >= 50,]

cor(subcombo$clpp_h4_coloc, subcombo$baseline_pval_baseline, use="complete.obs")
cor(subcombo$clpp_h4_coloc, subcombo$baseline_pval_baseline, use="complete.obs", method="spearman")
plot(subcombo$clpp_h4_coloc, subcombo$baseline_pval_baseline)
subcombo$baseline_rank = rank(subcombo$baseline_pval_baseline)
subcombo$h4_rank = rank(subcombo$clpp_h4_coloc)
plot(subcombo$baseline_rank, subcombo$h4_rank, xlab="Order of baseline results", ylab = "Order of COLOC results", pch=18)

# Find points that are very different between the two methods
disagreements = subcombo[na.omit(abs(subcombo$baseline_rank - subcombo$h4_rank) > 1000),]
# Sort them from most to least disagreeing
disagreements = disagreements[rev(order(abs(disagreements$baseline_rank - disagreements$h4_rank))),]
points(disagreements$baseline_rank, disagreements$h4_rank, col="red", pch=18)

###########################################
# A bit of other processing
###########################################

# Apply HEIDI test to throw away SMR tests that don't pass
combo$heidi_pval_smr[is.na(combo$heidi_pval_smr)] = 1
combo$smr_heidi_adjusted = combo$smr_neg_log_pval_smr
combo$smr_heidi_adjusted[combo$heidi_pval_smr < 0.05] = 0

###########################################
# UpSet plot for all
###########################################


# Make an array to show the different sets that we can have...
# Note: Right now this plot's pretty sloppy because we're not sure
# that the combo was actually tested in all methods...for some, it may
# have been dropped due to bugs in our pipeline rather than deficiencies
# of the method itself

combo$test_names = paste(combo$ref_snp, combo$eqtl_file, combo$gwas_trait, combo$base_gwas_file, combo$feature, sep="_")
stopifnot(length(combo$test_names) == length(unique(combo$test_names)))	# Each test should have a unique ID

upset_matrix = array(-1, dim=c(dim(combo)[1],8))
dimnames(upset_matrix)[[1]] = combo$test_names
dimnames(upset_matrix)[[2]] = c("FINEMAP_c1", "FINEMAP_c2", "COLOC", "Baseline_simple", "Baseline_smart", "SMR_no_HEIDI", "SMR_with_HEIDI", "GSMR")

combo$ranked_finemap_c1 = rank(-combo$clpp_finemap_c1)
combo$ranked_finemap_c2 = rank(-combo$clpp_finemap_c2)
combo$ranked_coloc = rank(-combo$clpp_h4_coloc)
combo$ranked_baseline_simple = rank(-combo$baseline_pval_baseline)
combo$ranked_baseline_smart = rank(-combo$baseline_pval5_baseline)
combo$ranked_smr_no_heidi = rank(-combo$smr_neg_log_pval_smr)
combo$ranked_smr_with_heidi = rank(-combo$smr_heidi_adjusted)
combo$ranked_gsmr_no_heidi = rank(-combo$smr_neg_log_pval_gsmr)

rank_threshold = 1000
for (i in 1:dim(upset_matrix)[1])
{
	upset_matrix[i, 1] = (combo[i,]$ranked_finemap_c1 <= rank_threshold)
	upset_matrix[i, 2] = (combo[i,]$ranked_finemap_c2 <= rank_threshold)
	upset_matrix[i, 3] = (combo[i,]$ranked_coloc <= rank_threshold)
	upset_matrix[i, 4] = (combo[i,]$ranked_baseline_simple <= rank_threshold)
	upset_matrix[i, 5] = (combo[i,]$ranked_baseline_smart <= rank_threshold)
	upset_matrix[i, 6] = (combo[i,]$ranked_smr_no_heidi <= rank_threshold)
	upset_matrix[i, 7] = (combo[i,]$ranked_smr_with_heidi <= rank_threshold)
	upset_matrix[i, 8] = (combo[i,]$ranked_gsmr_no_heidi <= rank_threshold)
}

upset_matrix[is.na(upset_matrix)] = 0
upset_matrix = upset_matrix[rowSums(upset_matrix) > 0,]
upset_matrix = data.frame(upset_matrix)
#upset_matrix$Names = rownames(upset_matrix)
upset_matrix_names = rownames(upset_matrix)

upset(upset_matrix, 
      sets = dimnames(upset_matrix)[[2]], 
      order.by="freq", matrix.color="blue", point.size=5,
      sets.bar.color=c("maroon","blue","orange", "green", "red", "violet", "cyan", "forestgreen"))

# The above plot may be misleading if we count COLOC and SMR results on methods for which these tests weren't even run
# due to insufficient data. Try subsetting down to the common denominator here.

# If at least one of the methods wasn't able to run on the study AT ALL, then get rid of that study
fair_studies = unique(combo[!is.na(combo$smr_neg_log_pval_smr),]$base_gwas_file)
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$clpp_h4_coloc),]$base_gwas_file)]
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$clpp_finemap_c1),]$base_gwas_file)]
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$clpp_finemap_c2),]$base_gwas_file)]
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$baseline_pval_baseline),]$base_gwas_file)]
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$baseline_pval5_baseline),]$base_gwas_file)]
fair_studies = fair_studies[fair_studies %in% unique(combo[!is.na(combo$smr_neg_log_pval_gsmr),]$base_gwas_file)]

fair_combo = combo[combo$base_gwas_file %in% fair_studies,]

upset_matrix = array(-1, dim=c(dim(fair_combo)[1],8))
dimnames(upset_matrix)[[1]] = fair_combo$test_names
dimnames(upset_matrix)[[2]] = c("FINEMAP_c1", "FINEMAP_c2", "COLOC", "Baseline_simple", "Baseline_smart", "SMR_no_HEIDI", "SMR_with_HEIDI", "GSMR")

fair_combo$ranked_finemap_c1 = rank(-fair_combo$clpp_finemap_c1)
fair_combo$ranked_finemap_c2 = rank(-fair_combo$clpp_finemap_c2)
fair_combo$ranked_coloc = rank(-fair_combo$clpp_h4_coloc)
fair_combo$ranked_baseline_simple = rank(-fair_combo$baseline_pval_baseline)
fair_combo$ranked_baseline_smart = rank(-fair_combo$baseline_pval5_baseline)
fair_combo$ranked_smr_no_heidi = rank(-fair_combo$smr_neg_log_pval_smr)
fair_combo$ranked_smr_with_heidi = rank(-fair_combo$smr_heidi_adjusted)
fair_combo$ranked_gsmr_no_heidi = rank(-fair_combo$smr_neg_log_pval_gsmr)

rank_threshold = 1000
for (i in 1:dim(upset_matrix)[1])
{
	upset_matrix[i, 1] = (fair_combo[i,]$ranked_finemap_c1 <= rank_threshold)
	upset_matrix[i, 2] = (fair_combo[i,]$ranked_finemap_c2 <= rank_threshold)
	upset_matrix[i, 3] = (fair_combo[i,]$ranked_coloc <= rank_threshold)
	upset_matrix[i, 4] = (fair_combo[i,]$ranked_baseline_simple <= rank_threshold)
	upset_matrix[i, 5] = (fair_combo[i,]$ranked_baseline_smart <= rank_threshold)
	upset_matrix[i, 6] = (fair_combo[i,]$ranked_smr_no_heidi <= rank_threshold)
	upset_matrix[i, 7] = (fair_combo[i,]$ranked_smr_with_heidi <= rank_threshold)
	upset_matrix[i, 8] = (fair_combo[i,]$ranked_gsmr_no_heidi <= rank_threshold)
}

upset_matrix[is.na(upset_matrix)] = 0
upset_matrix = upset_matrix[rowSums(upset_matrix) > 0,]
upset_matrix = data.frame(upset_matrix)
#upset_matrix$Names = rownames(upset_matrix)
upset_matrix_names = rownames(upset_matrix)

upset(upset_matrix, 
      sets = dimnames(upset_matrix)[[2]], 
      order.by="freq", matrix.color="blue", point.size=5,
      sets.bar.color=c("maroon","blue","orange", "green", "red", "violet", "cyan", "forestgreen"))






##########################################
# SVD to determine eigenmethods,
# using ranks of individual methods 
# for each site as the features
##########################################

# NOTE: I'm not 100% sure yet that SVD is the best
# way to do this since there are U and V matrices...
# think about this

rank_pairs = combo[c("ranked_finemap_c1", "ranked_finemap_c2", "ranked_coloc", "ranked_baseline_simple", "ranked_baseline_smart", "ranked_smr_no_heidi", "ranked_smr_with_heidi", "ranked_gsmr_no_heidi")]
eigenmethods = svd(rank_pairs)
em = eigenmethods$v
rownames(em) = colnames(rank_pairs)

fair_rank_pairs = fair_combo[c("ranked_finemap_c1", "ranked_finemap_c2", "ranked_coloc", "ranked_baseline_simple", "ranked_baseline_smart", "ranked_smr_no_heidi", "ranked_smr_with_heidi")]
fair_eigenmethods = svd(fair_rank_pairs)
fair_em = fair_eigenmethods$v
rownames(fair_em) = colnames(fair_rank_pairs)


##########################################
# Pairs plot for all methods
##########################################

pairs(rank_pairs)
pairs(fair_rank_pairs)

