# Author: Jeremy Tien

require(readr)
require(dplyr)
require(UpSetR)

# Load all colocalization results for different methods
combo = read.csv("/users/j29tien/colocalization_ML/real_matrix.tsv", sep='\t')
metadata = na.omit(read.csv("/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix_real_data.tsv", sep='\t'))[c(1:5)]
ensemble = read.csv("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/realdata_scores.tsv", sep='\t')
###########################################
# Pairs plot for all
###########################################

combo = cbind(ensemble, combo)
combo$test_names = paste(metadata$ref_snp, metadata$eqtl_file, metadata$gwas_trait, metadata$base_gwas_file, metadata$feature, sep="_")
stopifnot(length(combo$test_names) == length(unique(combo$test_names))) # Each test should have a unique ID

upset_matrix = array(-1, dim=c(dim(combo)[1],10))
dimnames(upset_matrix)[[1]] = combo$test_names
dimnames(upset_matrix)[[2]] = c('ENSEMBLE', 'COLOC', 'RTC', 'FINEMAP_STANDARD', 'FINEMAP_MODIFIED', 'BASELINE', 'SMARTBASELINE', 'SMR', 'GSMR', 'TWAS')

# rank scores from each method, greatest to least (hence the '-' sign)
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

upset_matrix[is.na(upset_matrix)] = 0
upset_matrix = upset_matrix[rowSums(upset_matrix) > 0,]
upset_matrix = data.frame(upset_matrix)
upset_matrix_names = rownames(upset_matrix)

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=3)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

# Plot and save
pdf(file="/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/allmethods_pairs.pdf", width=20, height=20)

pairs(combo[c("ranked_ENSEMBLE", "ranked_COLOC", "ranked_RTC", "ranked_FINEMAP_STANDARD", "ranked_FINEMAP_MODIFIED", "ranked_BASELINE", "ranked_SMARTBASELINE", "ranked_SMR", "ranked_GSMR", "ranked_TWAS")], lower.panel=panel.cor, upper.panel = upper.panel)

dev.off()
