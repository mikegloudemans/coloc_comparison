require(readr)
require(dplyr)
require(UpSetR)

# Load all colocalization results for different methods
combo = read.csv("/users/j29tien/colocalization_ML/real_matrix.tsv", sep='\t')
metadata = na.omit(read.csv("/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix_real_data.tsv", sep='\t'))[c(1:5)]
ensemble = read.csv("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/realdata_scores.tsv", sep='\t')
###########################################
# UpSet plot for all
###########################################


# Make an array to show the different sets that we can have...
# Note: Right now this plot's pretty sloppy because we're not sure
# that the combo was actually tested in all methods...for some, it may
# have been dropped due to bugs in our pipeline rather than deficiencies
# of the method itself

combo = cbind(ensemble, combo)
combo$test_names = paste(metadata$ref_snp, metadata$eqtl_file, metadata$gwas_trait, metadata$base_gwas_file, metadata$feature, sep="_")
stopifnot(length(combo$test_names) == length(unique(combo$test_names))) # Each test should have a unique ID

upset_matrix = array(-1, dim=c(dim(combo)[1],10))
dimnames(upset_matrix)[[1]] = combo$test_names
dimnames(upset_matrix)[[2]] = c('ENSEMBLE', 'COLOC', 'RTC', 'FINEMAP_STANDARD', 'FINEMAP_MODIFIED', 'BASELINE', 'SMARTBASELINE', 'SMR', 'GSMR', 'TWAS')

combo$ranked_ENSEMBLE = rank(-combo$ENSEMBLE)
combo$ranked_COLOC = rank(-combo$COLOC)
combo$ranked_RTC = rank(-combo$RTC)
combo$ranked_FINEMAP_STANDARD = rank(-combo$FINEMAP_STANDARD)
combo$ranked_FINEMAP_MODIFIED = rank(-combo$FINEMAP_MODIFIED)
combo$ranked_BASELINE = rank(-combo$BASELINE)
combo$ranked_SMARTBASELINE = rank(-combo$SMARTBASELINE)
combo$ranked_SMR = rank(-combo$SMR)
combo$ranked_GSMR = rank(-combo$GSMR)
combo$ranked_TWAS = rank(-combo$TWAS)

rank_threshold = 1000
for (i in 1:dim(upset_matrix)[1])
{
        upset_matrix[i, 1] = (combo[i,]$ranked_ENSEMBLE<= rank_threshold)
        upset_matrix[i, 2] = (combo[i,]$ranked_COLOC<= rank_threshold)
        upset_matrix[i, 3] = (combo[i,]$ranked_RTC<= rank_threshold)
        upset_matrix[i, 4] = (combo[i,]$ranked_FINEMAP_STANDARD<= rank_threshold)
        upset_matrix[i, 5] = (combo[i,]$ranked_FINEMAP_MODIFIED<= rank_threshold)
        upset_matrix[i, 6] = (combo[i,]$ranked_BASELINE<= rank_threshold)
        upset_matrix[i, 7] = (combo[i,]$ranked_SMARTBASELINE<= rank_threshold)
        upset_matrix[i, 8] = (combo[i,]$ranked_SMR<= rank_threshold)
        upset_matrix[i, 9] = (combo[i,]$ranked_GSMR<= rank_threshold)
        upset_matrix[i, 10] = (combo[i,]$ranked_TWAS<= rank_threshold)
}

upset_matrix[is.na(upset_matrix)] = 0
upset_matrix = upset_matrix[rowSums(upset_matrix) > 0,]
upset_matrix = data.frame(upset_matrix)
#upset_matrix$Names = rownames(upset_matrix)
upset_matrix_names = rownames(upset_matrix)

pdf(file="/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/allmethods_upset.pdf", width=15, height=10)

upset(upset_matrix, 
      sets = dimnames(upset_matrix)[[2]], 
      order.by="freq", matrix.color="blue", point.size=5,
      sets.bar.color=c("gold", "maroon","blue","orange", "green", "red", "violet", "cyan", "forestgreen", "black"))

dev.off()
