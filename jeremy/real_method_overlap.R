require(ggplot2)
require(reshape2)

data = read.table("colocalization_matrix_real_data.tsv", header=TRUE)
data = data[rowSums(is.na(data)) == 0,]
data = data[,6:14]
#data$twas_log_pval_twas[data$twas_log_pval_twas == 0] = NA
#data$smr_neg_log_pval_gsmr[data$smr_neg_log_pval_gsmr == 0] = NA

#data[data == apply(data, 2, min)] = NA
#data[data == apply(data, 2, max, na.rm=TRUE)] = NA

rank_matrix = apply(data, 2, function(x){rank(x, na.last="keep", ties.method="random")})
rank_matrix = rank_matrix[,which(colnames(rank_matrix) %in% c("clpp_h4_coloc", "clpp_finemap_c1", "smr_neg_log_pval_smr", "smr_neg_log_pval_gsmr", "twas_log_pval_twas", "rtc_score_rtc"))]
colnames(rank_matrix) = c("COLOC", "eCAVIAR", "SMR", "GSMR", "TWAS", "RTC")


top_sets = sapply(1:dim(rank_matrix)[2], function(x)
			{
				return(rank_matrix[,x] > quantile(rank_matrix[,x], 0.7, na.rm=TRUE))
			})
colnames(top_sets) = colnames(rank_matrix)

heat = list(method1 = rep("",dim(top_sets)[2]^2), method2 = rep("", dim(top_sets)[2]^2), overlap=rep(0, dim(top_sets)[2]^2))
index = 0
for (i in 1:(dim(top_sets)[2]-1))
{
	for (j in (i+1):dim(top_sets)[2])
	{
		index = index+1
		heat$method1[index] = colnames(top_sets)[i]
		heat$method2[index] = colnames(top_sets)[j]
		heat$overlap[index] = sum(which(top_sets[,i]) %in% which(top_sets[,j])) / sum(top_sets[,j], na.rm=TRUE)
	}
}
heat = as.data.frame(heat)
heat = heat[heat$method1 != "",]
heat$method1 = factor(heat$method1, levels = colnames(top_sets))
heat$method2 = factor(heat$method2, levels = colnames(top_sets))

png("real_rank_comparisons.png", units="in", height=6.5, width=6.5, res=360)
pairs(rank_matrix, pch=16, upper.panel = points, cex=0.3)
dev.off()

