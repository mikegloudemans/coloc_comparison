require(tidyverse)

# Evaluate robustness of this before presenting on them
#
# For example, will want to try including ranks rather than scores probably,
# because scores aren't normalized
#
# Can also try normalized scores though
#
# And throw away NA's and see to what extent that affects things...
#
# And what about heatmaps, which heatmaps are present the most in which PCs
#

max_pcs = 5

#####################################################
# Real
#####################################################

for (na_handling in c("na_fill_zero", "na_exclude"))
{
	real_data = read_delim("/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix_real_data.tsv", delim="\t")
	feats = real_data[,6:dim(real_data)[2]]

	if (na_handling == "na_fill_zero")
	{
		feats[is.na(feats)] = 0		# If it couldn't run, we've just got to fill it with 0 for now
	}
	if (na_handling == "na_exclude")
	{
		feats = feats[apply(feats, 1, function(x) {sum(is.na(x))}) == 0,]
	}

	stopifnot(sum(names(feats) != c("baseline_pval_baseline", "baseline_pval5_baseline", "clpp_h4_coloc", "clpp_finemap_c1",
					"clpp_mod_finemap_c1", "smr_neg_log_pval_smr", "smr_neg_log_pval_gsmr", "twas_log_pval_twas",
					"rtc_score_rtc"))== 0)
	names(feats) = c("baseline", "baseline++", "coloc", "finemap", "finemap-ld", "smr", "gsmr", "twas", "rtc")

	normed_feats = apply(feats, 2, scale)
	ranked_feats = apply(feats, 2, rank)
	normed_ranked_feats = apply(ranked_feats, 2, scale)
	# But be sure results are robust to this filling, even it we look at subsets of the results or whatever...

	pca = prcomp(feats)
	normed_pca = prcomp(normed_feats)
	ranked_pca = prcomp(ranked_feats)
	normed_ranked_pca = prcomp(normed_ranked_feats)

	pca_imp = round(summary(pca)$importance[2,]*100, 3)
	normed_pca_imp = round(summary(normed_pca)$importance[2,]*100, 3)
	ranked_pca_imp = round(summary(ranked_pca)$importance[2,]*100, 3)
	normed_ranked_pca_imp = round(summary(normed_ranked_pca)$importance[2,]*100, 3)

	print(pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/real_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(pca$rotation[,i], pca$rotation[,j], xlim = c(min(pca$rotation[,i]) - ((max(pca$rotation[,i]) - min(pca$rotation[,i]))*0.1),
										  max(pca$rotation[,i]) + ((max(pca$rotation[,i]) - min(pca$rotation[,i]))*0.1)),
									 ylim = c(min(pca$rotation[,j]) - ((max(pca$rotation[,j]) - min(pca$rotation[,j]))*0.1),
										  max(pca$rotation[,j]) + ((max(pca$rotation[,j]) - min(pca$rotation[,j]))*0.1)), main = paste("Real GWAS PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(pca_imp[j])))
				text(pca$rotation[,i], pca$rotation[,j] + 0.025, rownames(pca$rotation))
			dev.off()
		}
	}

	print(normed_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/real_normed_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(normed_pca$rotation[,i], normed_pca$rotation[,j], xlim = c(min(normed_pca$rotation[,i]) - ((max(normed_pca$rotation[,i]) - min(normed_pca$rotation[,i]))*0.1),
										  max(normed_pca$rotation[,i]) + ((max(normed_pca$rotation[,i]) - min(normed_pca$rotation[,i]))*0.1)),
									 ylim = c(min(normed_pca$rotation[,j]) - ((max(normed_pca$rotation[,j]) - min(normed_pca$rotation[,j]))*0.1),
										  max(normed_pca$rotation[,j]) + ((max(normed_pca$rotation[,j]) - min(normed_pca$rotation[,j]))*0.1)), main = paste("Real GWAS Normalized PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(normed_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(normed_pca_imp[j])))
				text(normed_pca$rotation[,i], normed_pca$rotation[,j] + 0.025, rownames(normed_pca$rotation))
			dev.off()
		}
	}

	print(ranked_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/real_ranked_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(ranked_pca$rotation[,i], ranked_pca$rotation[,j], xlim = c(min(ranked_pca$rotation[,i]) - ((max(ranked_pca$rotation[,i]) - min(ranked_pca$rotation[,i]))*0.1),
										  max(ranked_pca$rotation[,i]) + ((max(ranked_pca$rotation[,i]) - min(ranked_pca$rotation[,i]))*0.1)),
									 ylim = c(min(ranked_pca$rotation[,j]) - ((max(ranked_pca$rotation[,j]) - min(ranked_pca$rotation[,j]))*0.1),
										  max(ranked_pca$rotation[,j]) + ((max(ranked_pca$rotation[,j]) - min(ranked_pca$rotation[,j]))*0.1)), main = paste("Real GWAS Ranked PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(ranked_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(ranked_pca_imp[j])))
				text(ranked_pca$rotation[,i], ranked_pca$rotation[,j] + 0.025, rownames(ranked_pca$rotation))
			dev.off()
		}
	}

	print(normed_ranked_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/real_normed_ranked_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(normed_ranked_pca$rotation[,i], normed_ranked_pca$rotation[,j], xlim = c(min(normed_ranked_pca$rotation[,i]) - ((max(normed_ranked_pca$rotation[,i]) - min(normed_ranked_pca$rotation[,i]))*0.1),
										  max(normed_ranked_pca$rotation[,i]) + ((max(normed_ranked_pca$rotation[,i]) - min(normed_ranked_pca$rotation[,i]))*0.1)),
									 ylim = c(min(normed_ranked_pca$rotation[,j]) - ((max(normed_ranked_pca$rotation[,j]) - min(normed_ranked_pca$rotation[,j]))*0.1),
										  max(normed_ranked_pca$rotation[,j]) + ((max(normed_ranked_pca$rotation[,j]) - min(normed_ranked_pca$rotation[,j]))*0.1)), main = paste("Real GWAS Normalized Ranked PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(normed_ranked_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(normed_ranked_pca_imp[j])))
				text(normed_ranked_pca$rotation[,i], normed_ranked_pca$rotation[,j] + 0.025, rownames(normed_ranked_pca$rotation))
			dev.off()
		}
	}





	#####################################################
	# Simulations
	#####################################################

	sim_data = read_delim("/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix.tsv", delim="\t")
	sim_data = sim_data[,!grepl("rank", colnames(sim_data))]
	sim_data = sim_data[,!grepl("status", colnames(sim_data))]
	# If this isn't true, it means the file format has changed:
	stopifnot(sum(names(sim_data) != c("coloc_h4", "rtc_score", "finemap_clpp", "finemap_clppmod", "caviarbf_clpp",
				       "baseline_neglogpvalue", "smartbaseline_neglogpvalue", "smr_neglogpvalue", "smrheidiadjusted_neglogpvalue", "gsmr_neglogpvalue",
				       "twas_neglogpvalue")) == 0)
	names(sim_data) = c("coloc", "rtc", "finemap", "finemap_ld", "caviarbf", "baseline", "baseline++", "smr", "smr_heidi", "gsmr", "twas")

	sim_feats = sim_data
	if (na_handling == "na_fill_zero")
	{
		sim_feats[is.na(sim_feats)] = 0		# If it couldn't run, we've just got to fill it with 0 for now
	}

	if (na_handling == "na_exclude")
	{
		sim_feats = sim_feats[apply(sim_feats, 1, function(x) {sum(is.na(x))}) == 0,]
	}
	normed_sim_feats = apply(sim_feats, 2, scale)
	ranked_sim_feats = apply(sim_feats, 2, rank)
	normed_ranked_sim_feats = apply(ranked_sim_feats, 2, scale)
	# But be sure results are robust to this filling, even it we look at subsets of the results or whatever...

	sim_pca = prcomp(sim_feats)
	normed_sim_pca = prcomp(normed_sim_feats)
	ranked_sim_pca = prcomp(ranked_sim_feats)
	normed_ranked_sim_pca = prcomp(normed_ranked_sim_feats)

	sim_pca_imp = round(summary(sim_pca)$importance[2,]*100, 3)
	normed_sim_pca_imp = round(summary(normed_sim_pca)$importance[2,]*100, 3)
	ranked_sim_pca_imp = round(summary(ranked_sim_pca)$importance[2,]*100, 3)
	normed_ranked_sim_pca_imp = round(summary(normed_ranked_sim_pca)$importance[2,]*100, 3)

	print(sim_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/sim_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(sim_pca$rotation[,i], sim_pca$rotation[,j], xlim = c(min(sim_pca$rotation[,i]) - ((max(sim_pca$rotation[,i]) - min(sim_pca$rotation[,i]))*0.1),
										  max(sim_pca$rotation[,i]) + ((max(sim_pca$rotation[,i]) - min(sim_pca$rotation[,i]))*0.1)),
									 ylim = c(min(sim_pca$rotation[,j]) - ((max(sim_pca$rotation[,j]) - min(sim_pca$rotation[,j]))*0.1),
										  max(sim_pca$rotation[,j]) + ((max(sim_pca$rotation[,j]) - min(sim_pca$rotation[,j]))*0.1)), main = paste("Simulated GWAS PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(sim_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(sim_pca_imp[j])))
				text(sim_pca$rotation[,i], sim_pca$rotation[,j] + 0.025, rownames(sim_pca$rotation))
			dev.off()
		}
	}

	print(normed_sim_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/normed_sim_pca_pc", as.character(i), "_pc", as.character(j),  "_", na_handling, ".pdf"), height=5, width=5)		
				plot(normed_sim_pca$rotation[,i], normed_sim_pca$rotation[,j], xlim = c(min(normed_sim_pca$rotation[,i]) - ((max(normed_sim_pca$rotation[,i]) - min(normed_sim_pca$rotation[,i]))*0.1),
										  max(normed_sim_pca$rotation[,i]) + ((max(normed_sim_pca$rotation[,i]) - min(normed_sim_pca$rotation[,i]))*0.1)),
									 ylim = c(min(normed_sim_pca$rotation[,j]) - ((max(normed_sim_pca$rotation[,j]) - min(normed_sim_pca$rotation[,j]))*0.1),
										  max(normed_sim_pca$rotation[,j]) + ((max(normed_sim_pca$rotation[,j]) - min(normed_sim_pca$rotation[,j]))*0.1)), main = paste("Simulated GWAS Normalized PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(normed_sim_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(normed_sim_pca_imp[j])))
				text(normed_sim_pca$rotation[,i], normed_sim_pca$rotation[,j] + 0.025, rownames(normed_sim_pca$rotation))
			dev.off()
		}
	}

	print(ranked_sim_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/ranked_sim_pca_pc", as.character(i), "_pc", as.character(j),  "_", na_handling, ".pdf"), height=5, width=5)		
				plot(ranked_sim_pca$rotation[,i], ranked_sim_pca$rotation[,j], xlim = c(min(ranked_sim_pca$rotation[,i]) - ((max(ranked_sim_pca$rotation[,i]) - min(ranked_sim_pca$rotation[,i]))*0.1),
										  max(ranked_sim_pca$rotation[,i]) + ((max(ranked_sim_pca$rotation[,i]) - min(ranked_sim_pca$rotation[,i]))*0.1)),
									 ylim = c(min(ranked_sim_pca$rotation[,j]) - ((max(ranked_sim_pca$rotation[,j]) - min(ranked_sim_pca$rotation[,j]))*0.1),
										  max(ranked_sim_pca$rotation[,j]) + ((max(ranked_sim_pca$rotation[,j]) - min(ranked_sim_pca$rotation[,j]))*0.1)), main = paste("Simulated GWAS Ranked PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(ranked_sim_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(ranked_sim_pca_imp[j])))
				text(ranked_sim_pca$rotation[,i], ranked_sim_pca$rotation[,j] + 0.025, rownames(ranked_sim_pca$rotation))
			dev.off()
		}
	}

	print(normed_ranked_sim_pca$sdev)
	for (i in 1:max_pcs)
	{
		for (j in 1:max_pcs)
		{
			pdf(paste0("/users/mgloud/projects/coloc_comparisons/output/plots/pca/normed_ranked_sim_pca_pc", as.character(i), "_pc", as.character(j), "_", na_handling, ".pdf"), height=5, width=5)		
				plot(normed_ranked_sim_pca$rotation[,i], normed_ranked_sim_pca$rotation[,j], xlim = c(min(normed_ranked_sim_pca$rotation[,i]) - ((max(normed_ranked_sim_pca$rotation[,i]) - min(normed_ranked_sim_pca$rotation[,i]))*0.1),
										  max(normed_ranked_sim_pca$rotation[,i]) + ((max(normed_ranked_sim_pca$rotation[,i]) - min(normed_ranked_sim_pca$rotation[,i]))*0.1)),
									 ylim = c(min(normed_ranked_sim_pca$rotation[,j]) - ((max(normed_ranked_sim_pca$rotation[,j]) - min(normed_ranked_sim_pca$rotation[,j]))*0.1),
										  max(normed_ranked_sim_pca$rotation[,j]) + ((max(normed_ranked_sim_pca$rotation[,j]) - min(normed_ranked_sim_pca$rotation[,j]))*0.1)), main = paste("Simulated GWAS Normalized Ranked PCA"),
									 xlab = paste0("PC", as.character(i), ", variance explained: ", as.character(normed_ranked_sim_pca_imp[i]), "%"), ylab=paste0("PC", as.character(j), ", variance explained: ", as.character(normed_ranked_sim_pca_imp[j])))
				text(normed_ranked_sim_pca$rotation[,i], normed_ranked_sim_pca$rotation[,j] + 0.025, rownames(normed_ranked_sim_pca$rotation))
			dev.off()
		}
	}
}

