
# For each effect size...
for (es in c(0.3, 0.005, 0.05, 0.1, 0.2, 0.3))
{
	# For each locus...
	for (locus in 0:4)
	{
		png(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/es_sims/effectsize", es, "_locus", locus, ".png"), height = 4, width = 15, units="in", res=360)
		par(mfcol=c(2,5))

		# For each hypothesis...
		for (hypothesis in 0:4)
		{
			timestamp=dir(paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/vary_effect_size/h", hypothesis,"/es_", es, "/"))

			answer_key = read.table(paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/vary_effect_size/h", hypothesis, "/es_", es, "/", timestamp, "/answer_key.txt"), header=TRUE, stringsAsFactors=FALSE, sep="\t")

			eqtl = read.table(paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/vary_effect_size/h", hypothesis, "/es_", es, "/", timestamp, "/hg19/eqtl/eqtl_sumstats", locus, ".txt.gz"),header=TRUE)
			gwas = read.table(paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/vary_effect_size/h", hypothesis, "/es_", es, "/", timestamp, "/hg19/gwas/gwas_sumstats", locus, ".txt.gz"),header=TRUE)

			plot(eqtl$snp_pos, -log10(eqtl$pvalue), pch=18, xlab="Position", ylab="eQTL -log10(p-value)", main=paste0("H", hypothesis))

			if (hypothesis == 1)
			{
				gwas_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				gwas_pos = gwas[gwas$rsid == gwas_causal,]$snp_pos
				abline(v=gwas_pos, col="blue", lwd=2)
			} else if (hypothesis == 2)
			{
				eqtl_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "eqtl:")[[1]][2], ":")[[1]][1]
				eqtl_pos = gwas[eqtl$rsid == eqtl_causal,]$snp_pos
				abline(v=eqtl_pos, col="red", lwd=2)
			} else if (hypothesis == 3)
			{
				gwas_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				gwas_pos = gwas[gwas$rsid == gwas_causal,]$snp_pos
				abline(v=gwas_pos, col="blue", lwd=2)
				eqtl_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "eqtl:")[[1]][2], ":")[[1]][1]
				eqtl_pos = gwas[eqtl$rsid == eqtl_causal,]$snp_pos
				abline(v=eqtl_pos, col="red", lwd=2)
			} else if (hypothesis == 4)
			{
				both_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				both_pos = gwas[gwas$rsid == both_causal,]$snp_pos
				abline(v=both_pos, col="black", lwd=2)
			}

			plot(gwas$snp_pos, -log10(gwas$pvalue), pch=18, xlab="Position", ylab="GWAS -log10(p-value)")

			if (hypothesis == 1)
			{
				gwas_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				gwas_pos = gwas[gwas$rsid == gwas_causal,]$snp_pos
				abline(v=gwas_pos, col="blue", lwd=2)
			} else if (hypothesis == 2)
			{
				eqtl_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "eqtl:")[[1]][2], ":")[[1]][1]
				eqtl_pos = gwas[eqtl$rsid == eqtl_causal,]$snp_pos
				abline(v=eqtl_pos, col="red", lwd=2)
			} else if (hypothesis == 3)
			{
				gwas_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				gwas_pos = gwas[gwas$rsid == gwas_causal,]$snp_pos
				abline(v=gwas_pos, col="blue", lwd=2)
				eqtl_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "eqtl:")[[1]][2], ":")[[1]][1]
				eqtl_pos = gwas[eqtl$rsid == eqtl_causal,]$snp_pos
				abline(v=eqtl_pos, col="red", lwd=2)
			} else if (hypothesis == 4)
			{
				both_causal = strsplit(strsplit(answer_key$causal_variants[locus+1], "gwas:")[[1]][2], ":")[[1]][1]
				both_pos = gwas[gwas$rsid == both_causal,]$snp_pos
				abline(v=both_pos, col="black", lwd=2)
			}

		}
		#readline("")
		dev.off()
	}
}
