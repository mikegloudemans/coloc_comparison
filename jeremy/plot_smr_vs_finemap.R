data = read.table("colocalization_matrix_real_data.tsv", header=TRUE)

png("ecaviar-smr.png", res=600, units="in", width=5, height=5.2)
plot(rank(data$clpp_finemap_c1 + runif(dim(data)[1]) * 0.0000001), rank(data$smr_neg_log_pval_smr + runif(dim(data)[1]) * 0.0000001), pch=16, cex=0.8, xlab="", ylab="")
cor(rank(data$clpp_finemap_c1), rank(data$smr_neg_log_pval_smr))
dev.off()

