# Load all results, concatenate them

# TODO: We also need to make sure the set being compared by the different methods is the same
# (should really just use the results being compared with the ensemble method). Because if a 
# site is missing it's not obvious what info that contains right now...and I suspect that
# it's actually making it easier on something like COLOC if this site is just thrown away
# since it's likely to be a negative

require(pROC)
require(PRROC)

main = function()
{
	timestamp = "2018-07-27_15-23-15"
	answer_file = paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/", timestamp, "/answer_key.txt")

	answer_key = get_answer_key(answer_file)

	# Run with full answer key
	compare_methods(answer_key, timestamp)

	# What if we only allow our "test" colocalizations to be those with GWAS
	# effect size over a certain threshold?	
	other_answer_key = filter_answer_key_by_gwas_effect(answer_key, 0.05)
	compare_methods(other_answer_key, timestamp)

	# Those with eQTL effect size over a certain threshold?
	other_answer_key = filter_answer_key_by_eqtl_effect(answer_key, 0.1)
	compare_methods(other_answer_key, timestamp)

	# TODO: What if MAX GWAS threshold is regulated? After all, we don't often see so
	# many great-looking sites, more often we have mostly junk

	# Those with GWAS/eQTL sample size over a certain threshold?
	other_answer_key = filter_answer_key_by_min_gwas_case_n(answer_key, 500)
	other_answer_key = filter_answer_key_by_min_gwas_control_n(other_answer_key, 500)
	other_answer_key = filter_answer_key_by_min_eqtl_n(other_answer_key, 100)
	compare_methods(other_answer_key, timestamp)

	# Those with GWAS/eQTL sample size UNDER a certain threshold
	# (special case for underpowered GWAS or eQTL studies)
	other_answer_key = filter_answer_key_by_max_eqtl_n(answer_key, 100)
	compare_methods(other_answer_key, timestamp)

	# What about if we change the ratio of positive to negative sites?
	other_answer_key = trim_answer_key(answer_key, 0.5)
	compare_methods(other_answer_key, timestamp)

	other_answer_key = trim_answer_key(answer_key, 0.25)
	compare_methods(other_answer_key, timestamp)

	other_answer_key = trim_answer_key(answer_key, 3)
	compare_methods(other_answer_key, timestamp)
	
}

compare_methods = function(answer_key, timestamp)
{
	finemap_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/finemap-comparisons/", timestamp)
	caviarbf_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/caviarbf-comparisons/", timestamp)
	coloc_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/coloc-comparisons/", timestamp)
	rtc_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/rtc-comparisons/", timestamp)
	twas_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/twas-comparisons/", timestamp)
	bl_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/baseline-comparisons/", timestamp)
	smr_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/smr-comparisons/", timestamp)
	gsmr_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/gsmr-comparisons/", timestamp)

	# 
	# Part 1: Compare methods' performance on simulated data
	#
	answer_indices = answer_key$test_case

	# Get the test sites for which all methods were tried.
	# TODO: Figure out why some sites aren't being tested -- technically this
	# should be a penalty to THAT method; doesn't necessarily mean that every
	# method should throw away that site
	finemap_results = get_finemap_results(finemap_base_dir) 
	finemap_problems = gsub("gwas_sumstats", "", finemap_results$base_gwas_file)
	finemap_problems = as.numeric(gsub("_txt_gz", "", finemap_problems))
	coloc_results = get_coloc_results(coloc_base_dir)
	coloc_problems = gsub("eqtl_sumstats", "", coloc_results$eqtl_file)
	coloc_problems = as.numeric(gsub("_txt_gz", "", coloc_problems))
	rtc_results = get_rtc_results(rtc_base_dir)
	rtc_problems = gsub("eqtl_sumstats", "", rtc_results$eqtl_file)
	rtc_problems = as.numeric(gsub("_txt_gz", "", rtc_problems))
	bf_results = get_caviarbf_results(caviarbf_base_dir)
	bf_problems = gsub("eqtl_sumstats", "", bf_results$eqtl_file)
	bf_problems = as.numeric(gsub("_txt_gz", "", bf_problems))
	bl_results = get_baseline_results(bl_base_dir)
	bl_problems = gsub("eqtl_sumstats", "", bl_results$eqtl_file)
	bl_problems = as.numeric(gsub("_txt_gz", "", bl_problems))
	twas_results = get_twas_results(twas_base_dir)
	twas_problems = gsub("eqtl_sumstats", "", twas_results$eqtl_file)
	twas_problems = as.numeric(gsub("_txt_gz", "", twas_problems))
	smr_results = get_smr_results(smr_base_dir)
	smr_problems = gsub("eqtl_sumstats", "", smr_results$eqtl_file)
	smr_problems = as.numeric(gsub("_txt_gz", "", smr_problems))
	gsmr_results = get_gsmr_results(gsmr_base_dir)
	gsmr_problems = gsub("eqtl_sumstats", "", gsmr_results$eqtl_file)
	gsmr_problems = as.numeric(gsub("_txt_gz", "", gsmr_problems))

	# Apply HEIDI test for SMR
	heidi_fails = smr_results$heidi_pval <= 0.05
	heidi_fails[is.na(heidi_fails)] = FALSE
	smr_results$heidi_adjusted_pval = smr_results$smr_neg_log_pval
	smr_results$heidi_adjusted_pval[heidi_fails] = 0
	

	problems = rtc_problems[rtc_problems %in% finemap_problems]
	problems = problems[problems %in% coloc_problems]
	problems = problems[problems %in% bf_problems]
	problems = problems[problems %in% bl_problems]
	problems = problems[problems %in% smr_problems]
	problems = problems[problems %in% gsmr_problems]
	problems = problems[problems %in% twas_problems]
	problems = problems[problems %in% answer_key$test_case]

	selection = match(answer_indices, problems)
	matched_answers = answer_key[which(!is.na(selection)),]
	answers = matched_answers$is_coloc
	selection = selection[!is.na(selection)]
	problem_set = problems[selection]

	h4pp = scale(coloc_results$clpp_h4[match(problem_set, coloc_problems)])
	rtc = scale(rtc_results$rtc[match(problem_set, rtc_problems)])
	clpp = scale(finemap_results$clpp[match(problem_set, finemap_problems)])
	clpp_mod = scale(finemap_results$clpp_mod[match(problem_set, finemap_problems)])
	bf_clpp = scale(bf_results$clpp[match(problem_set, bf_problems)])
	bl = scale(bl_results$baseline_pval[match(problem_set, bl_problems)])
	bl5 = scale(bl_results$baseline_pval5[match(problem_set, bl_problems)])
	smr = scale(smr_results$smr_neg_log_pval[match(problem_set, smr_problems)])
	smr_heidi = scale(smr_results$heidi_adjusted_pval[match(problem_set, smr_problems)])
	gsmr = scale(gsmr_results$smr_neg_log_pval[match(problem_set, gsmr_problems)])
	twas_p = scale(twas_results$twas_log_pval[match(problem_set, twas_problems)])

	plot(roc(answers, as.numeric(clpp)), print.auc = TRUE, col = "black", print.auc.x = 0.2, print.auc.y = 0.32, main = "Colocalization detection performance")
	plot(roc(answers, as.numeric(clpp_mod)), print.auc = TRUE, col = "red", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.28)
	plot(roc(answers, as.numeric(h4pp)), print.auc = TRUE, col = "blue", add=TRUE, print.auc.x = 0.2, print.auc.y = 0.24)
	plot(roc(answers, as.numeric(rtc)), print.auc = TRUE, col = "green", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.20)
	plot(roc(answers, as.numeric(bf_clpp)), print.auc = TRUE, col = "gray", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.16)
	plot(roc(answers, as.numeric(bl)), print.auc = TRUE, col = "turquoise3", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.04)
	plot(roc(answers, as.numeric(bl5)), print.auc = TRUE, col = "yellowgreen", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.36)
	plot(roc(answers, as.numeric(smr)), print.auc = TRUE, col = "darkgoldenrod4", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.40)
	plot(roc(answers, as.numeric(smr_heidi)), print.auc = TRUE, col = "deeppink2", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.44)
	plot(roc(answers, as.numeric(gsmr)), print.auc = TRUE, col = "forestgreen", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.48)
	plot(roc(answers, as.numeric(twas_p)), print.auc = TRUE, col = "orange", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.12)

	# Make colocalization matrix for Jeremy to mess around with for ensemble approach
	colocalization_matrix = data.frame(list(coloc_h4=coloc_results$clpp_h4[match(problem_set, coloc_problems)],
						coloc_h4_rank=rank(-coloc_results$clpp_h4[match(problem_set, coloc_problems)]),
						rtc_score=rtc_results$rtc[match(problem_set, rtc_problems)],
						rtc_score_rank=rank(-rtc_results$rtc[match(problem_set, rtc_problems)]),
						finemap_clpp=finemap_results$clpp[match(problem_set, finemap_problems)],
						finemap_clpp_rank=rank(-finemap_results$clpp[match(problem_set, finemap_problems)]),
						finemap_clppmod=finemap_results$clpp_mod[match(problem_set, finemap_problems)],
						finemap_clppmod_rank=rank(-finemap_results$clpp_mod[match(problem_set, finemap_problems)]),
						caviarbf_clpp=bf_results$clpp[match(problem_set, bf_problems)],
						caviarbf_clpp_rank=rank(-bf_results$clpp[match(problem_set, bf_problems)]),
						baseline_neglogpvalue=bl_results$baseline_pval[match(problem_set, bl_problems)],
						baseline_neglogpvalue_rank=rank(-bl_results$baseline_pval[match(problem_set, bl_problems)]),
						smartbaseline_neglogpvalue=bl_results$baseline_pval5[match(problem_set, bl_problems)],
						smartbaseline_neglogpvalue_rank=rank(-bl_results$baseline_pval5[match(problem_set, bl_problems)]),
						smr_neglogpvalue=smr_results$smr_neg_log_pval[match(problem_set, smr_problems)],
						smr_neglogpvalue_rank=rank(-smr_results$smr_neg_log_pval[match(problem_set, smr_problems)]),
						smrheidiadjusted_neglogpvalue=smr_results$heidi_adjusted_pval[match(problem_set, smr_problems)],
						smrheidiadjusted_neglogpvalue_rank=rank(-smr_results$heidi_adjusted_pval[match(problem_set, smr_problems)]),
						gsmr_neglogpvalue=gsmr_results$smr_neg_log_pval[match(problem_set, gsmr_problems)],
						gsmr_neglogpvalue_rank=rank(-gsmr_results$smr_neg_log_pval[match(problem_set, gsmr_problems)]),
						twas_neglogpvalue=twas_results$twas_log_pval[match(problem_set, twas_problems)],
						twas_neglogpvalue_rank=rank(-twas_results$twas_log_pval[match(problem_set, twas_problems)]),
						colocalization_status=answers))
	write.table(colocalization_matrix, "/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")	

	# With smart baseline included
	ensemble = h4pp + rtc + clpp + clpp_mod + bf_clpp + bl5 + smr + gsmr + twas_p
	# Without baseline included
	ensemble = h4pp + rtc + clpp + clpp_mod + bf_clpp + smr + gsmr + twas_p
	
	ensemble = ensemble / max(ensemble)
	plot(roc(answers, ensemble), print.auc = TRUE, col = "purple", add=TRUE, print.auc.x = 0.2, print.auc.y = 0.08 )
	legend(0.8, 0.4, legend=c("FINEMAP-CLPP", "FINEMAP-CLPP_mod", "COLOC", "RTC", "CAVIARBF", "TWAS", "ensemble (no baseline)", "baseline", "smart-baseline", "SMR", "SMR-HEIDI", "GSMR"),
	              col=c("black", "red", "blue", "green", "gray", "orange", "purple", "turquoise3", "yellowgreen", "darkgoldenrod4", "deeppink2", "forestgreen"), bg="white", lty=1, cex=0.8, lwd=4)
	# Naive ensemble performance is comparable with the best methods

	readline("Press enter:")

	plot(pr.curve(clpp[which(answers==1)], clpp[which(answers==0)], curve=TRUE), color="black")
	plot(pr.curve(clpp_mod[which(answers==1)], clpp_mod[which(answers==0)], curve=TRUE), add=TRUE, color="red")
	plot(pr.curve(h4pp[which(answers==1)], h4pp[which(answers==0)], curve=TRUE), add=TRUE, color="blue")
	plot(pr.curve(rtc[which(answers==1)], rtc[which(answers==0)], curve=TRUE), add=TRUE, color="green")
	plot(pr.curve(bf_clpp[which(answers==1)], bf_clpp[which(answers==0)], curve=TRUE), add=TRUE, color="gray")
	plot(pr.curve(bl[which(answers==1)], bl[which(answers==0)], curve=TRUE), add=TRUE, color="turquoise3")
	plot(pr.curve(bl5[which(answers==1)], bl5[which(answers==0)], curve=TRUE), add=TRUE, color="yellowgreen")
	#plot(pr.curve(twas_p[which(answers==1)], twas_p[which(answers==0)], curve=TRUE), add=TRUE, color="orange")
	plot(pr.curve(ensemble[which(answers==1)], ensemble[which(answers==0)], curve=TRUE), add=TRUE, color="purple")
	plot(pr.curve(smr[which(answers==1)], ensemble[which(answers==0)], curve=TRUE), add=TRUE, color="darkgoldenrod4")
	plot(pr.curve(smr_heidi[which(answers==1)], ensemble[which(answers==0)], curve=TRUE), add=TRUE, color="deeppink2")
	plot(pr.curve(gsmr[which(answers==1)], ensemble[which(answers==0)], curve=TRUE), add=TRUE, color="forestgreen")

	readline("Press enter:")


	# Note of course that any site causing one or more methods to fail completely
	# will not currently appear in this table.
	#results_table = data.frame(list(coloc_h4pp=h4pp, rtc=rtc, finemap_clpp=clpp, finemap_clpp_mod=clpp_mod, caviar_bf_clpp=bf_clpp, twas_logp=twas_p))
	results_table = data.frame(list(coloc_h4pp=h4pp, rtc=rtc, finemap_clpp=clpp, finemap_clpp_mod=clpp_mod, caviar_bf_clpp=bf_clpp, baseline_p=bl, baseline3_p=bl3))
	# These two plots are illuminating. Goal will be to pick out some of the off-diagonal elements and to figure out
	# why they're scoring differently in the different methods. Particularly for the very different methods.
	pairs(results_table, lower.panel=NULL)
	readline("Press enter:")
	cols = c("black","red")[as.numeric(answers) + 1]
	pairs(apply(results_table, 2, rank), lower.panel=NULL, col=cols, pch=19)
	readline("Press enter:")
	
	return(FALSE)

	#
	# Part 2: Modify parameters of the methods themselves
	#

	# Now let's try it with a few other adjustments (will look better once
	# we have more sites though...)

	# What if we penalize scores with insignificant p-values?
	# Set FINEMAP to a very small value if log_eqtl_pval or log_gwas_pval
	# is below certain thresholds
	# - TODO: Plot: what happens to CLPP score as we go from no cutoff to imposing a gradually increasing cutoff on the scores?
	# - TODO: Plot: what happens to CLPP_mod score?
	# - TODO: Is this effect sensitive to the simulation parameters for generating test cases?
	# - TODO: What about if we apply the same cutoffs with other methods? Coloc, RTC, SMR...?

	# NOTE: This part needs to be fixed to deal with the adaptively sized answer key

	require(grDevices)

	dup_finemap_results = finemap_results
	clpp = dup_finemap_results$clpp
	clpp_mod = dup_finemap_results$clpp_mod
	problems = gsub("gwas_sumstats", "", dup_finemap_results$base_gwas_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answer_key = get_answer_key(answer_file)
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
	for (cutoff in 0:5)
	{
		print(cutoff)
		clpp_color = adjustcolor("red", alpha.f = 0.7^cutoff)
		clpp_mod_color = adjustcolor("blue", alpha.f = 0.7^cutoff)

		gwas_cutoff = cutoff
		eqtl_cutoff = cutoff

		dup_finemap_results = finemap_results
		dup_finemap_results$adjusted = pmin(0.0000001 * dup_finemap_results$X.log_eqtl_pval, 0.0000001 * dup_finemap_results$X.log_gwas_pval)
		dup_finemap_results[(dup_finemap_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_finemap_results$X.log_gwas_pval < gwas_cutoff),]$clpp = dup_finemap_results[(dup_finemap_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_finemap_results$X.log_gwas_pval < gwas_cutoff),]$adjusted
		dup_finemap_results[(dup_finemap_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_finemap_results$X.log_gwas_pval < gwas_cutoff),]$clpp_mod = dup_finemap_results[(dup_finemap_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_finemap_results$X.log_gwas_pval < gwas_cutoff),]$adjusted
		clpp = dup_finemap_results$clpp
		clpp_mod = dup_finemap_results$clpp_mod
		problems = gsub("gwas_sumstats", "", dup_finemap_results$base_gwas_file)
		problems = as.numeric(gsub("_txt_gz", "", problems))
		answer_key = get_answer_key(answer_file)
		answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
		
		# Make eCAVIAR ROC
		plot(roc(answers, clpp), add=cutoff!=0, col=clpp_color, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.60 - 0.04*cutoff)
		plot(roc(answers, clpp_mod), add=TRUE, col=clpp_mod_color, print.auc=TRUE, print.auc.x = 0.2, print.auc.y = 0.36 - 0.04*cutoff)	
	}
	readline("Press enter:")
	# Possibly add legend to the above? This coloring isn't good though because the alphas
	# change color when they overlap
}

# Right now some could be sampled twice, which is admittedly not ideal
# We set the ease factor by including all regular coloc sites as normal
# but by resampling the desired number of negative examples to get the
# proportion correct
trim_answer_key = function(answer_key, factor)
{
	set.seed(0)
	coloc_key = answer_key[answer_key$is_coloc,]
	nocoloc_key = answer_key[!answer_key$is_coloc,]
	easy_nocoloc_key = nocoloc_key[sample(1:dim(nocoloc_key)[1], floor(dim(nocoloc_key)[1] * factor), replace=TRUE),]
	easy_answer_key = rbind(coloc_key, easy_nocoloc_key)
	return(easy_answer_key)
}

get_answer_key = function(answer_file)
{
	answer_key = read.table(answer_file, header=TRUE, sep="\t")
	is_coloc = function(x)
	{
		vars = strsplit(x, ",")
		eqtls = c()
		gwas = c()
		for (v in vars[[1]])
		{
			info = strsplit(v, ":")
			if (info[[1]][1] == "gwas")
			{
				gwas = c(gwas, info[[1]][2])
			} else if (info[[1]][1] == "eqtl")
			{
				eqtls = c(eqtls, info[[1]][2])
			}
		}
		overlap = match(eqtls, gwas)
		overlap = overlap[!is.na(overlap)]
		return(length(overlap) > 0)
	}

	answer_key$is_coloc = sapply(as.character(answer_key$causal_variants), is_coloc)

	return(answer_key)
}

filter_answer_key_by_min_gwas_case_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$cases_n >= N,]
	return(new_answer_key)
}

filter_answer_key_by_min_gwas_control_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$controls_n >= N,]
	return(new_answer_key)
}

filter_answer_key_by_min_eqtl_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$eqtl_n >= N,]
	return(new_answer_key)
}

filter_answer_key_by_max_gwas_case_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$cases_n <= N,]
	return(new_answer_key)
}

filter_answer_key_by_max_gwas_control_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$controls_n <= N,]
	return(new_answer_key)
}

filter_answer_key_by_max_eqtl_n = function(answer_key, N)
{
	new_answer_key = answer_key[answer_key$eqtl_n <= N,]
	return(new_answer_key)
}

filter_answer_key_by_max_cc_ratio = function(answer_key, r)
{
	new_answer_key = answer_key[answer_key$cases_n / (answer_key$controls_n + answer_key$cases_n) <= r,]
	return(new_answer_key)
}

filter_answer_key_by_min_cc_ratio = function(answer_key, r)
{
	new_answer_key = answer_key[answer_key$cases_n / (answer_key$controls_n + answer_key$cases_n) >= r,]
	return(new_answer_key)
}


# Any so-called colocalization that has GWAS effect size below a 
# certain level will be removed
filter_answer_key_by_gwas_effect = function(answer_key, level)
{
	exceeds_level = function(x)
	{
		vars = strsplit(x, ",")
		eqtls = c()
		gwas = c()
		for (v in vars[[1]])
		{
			info = strsplit(v, ":")
			if (info[[1]][1] == "gwas")
			{
				return(abs(as.numeric(info[[1]][3])) > level)
			}
		}
		return(FALSE)
	}
	new_answer_key = answer_key[sapply(as.character(answer_key$causal_variants), exceeds_level) | !answer_key$is_coloc,]
	return(new_answer_key)
}

# Any so-called colocalization that has GWAS effect size below a 
# certain level will be removed
filter_answer_key_by_eqtl_effect = function(answer_key, level)
{
	exceeds_level = function(x)
	{
		vars = strsplit(x, ",")
		eqtls = c()
		gwas = c()
		for (v in vars[[1]])
		{
			info = strsplit(v, ":")
			if (info[[1]][1] == "eqtl")
			{
				return(abs(as.numeric(info[[1]][3])) > level)
			}
		}
		return(FALSE)
	}
	new_answer_key = answer_key[sapply(as.character(answer_key$causal_variants), exceeds_level) | !answer_key$is_coloc,]
	return(new_answer_key)
}


get_finemap_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("clpp_status", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)

}

get_coloc_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("h4pp_status", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_rtc_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("rtc_score_status", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_caviarbf_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("clpp_status", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_smr_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("_smr", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_gsmr_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("gsmr", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_twas_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("twas_clpp_status", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

get_baseline_results = function(base_dir)
{
	results_dirs = dir(base_dir)
	data = list()
	for (rd in results_dirs)
	{
		subdir = dir(paste(base_dir, rd, sep="/"))
		results_file = subdir[grepl("baseline", subdir)]
		data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
	}

	all_results = do.call(rbind, data)
	return(all_results)
}

# This function is never called and is currently used for nothing
# LOL
nothing = function()
{
	#
	# Part 3: Modify parameters of the simulation
	#


	# What if we change the ratio of cases to controls?
	sum(answer_key$is_coloc) / length(answer_key$is_coloc)

	set.seed(0)
	coloc_key = answer_key[answer_key$is_coloc,]
	nocoloc_key = answer_key[!answer_key$is_coloc,]
	easy_nocoloc_key = nocoloc_key[sample(1:dim(nocoloc_key)[1], floor(dim(nocoloc_key)[1] / 2)),]
	easy_answer_key = rbind(coloc_key, easy_nocoloc_key)
	easiest_nocoloc_key = nocoloc_key[sample(1:dim(nocoloc_key)[1], floor(dim(nocoloc_key)[1] / 4)),]
	easiest_answer_key = rbind(coloc_key, easiest_nocoloc_key)

	# eCAVIAR/FINEMAP (should actually run both of them, because why not)
	# Get CLPP scores and actual answers for each test
	problems = gsub("gwas_sumstats", "", finemap_results$base_gwas_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	easy_lineup = match(problems, easy_answer_key$test_case)
	easy_sub_problems = problems[!is.na(easy_lineup)]
	easy_clpp = finemap_results[which(problems %in% easy_sub_problems),]$clpp
	easy_answers = easy_answer_key$is_coloc[easy_lineup[!is.na(easy_lineup)]]
	easiest_lineup = match(problems, easiest_answer_key$test_case)
	easiest_sub_problems = problems[!is.na(easiest_lineup)]
	easiest_clpp = finemap_results[which(problems %in% easiest_sub_problems),]$clpp
	easiest_answers = easiest_answer_key$is_coloc[easiest_lineup[!is.na(easiest_lineup)]]
	sum(easy_answers) / length(easy_answers)
	sum(easiest_answers) / length(easiest_answers)
	
	# Make eCAVIAR ROC
	roc(easy_answers, easy_clpp, plot=TRUE)
	roc(easiest_answers, easiest_clpp, plot=TRUE)

	# Coloc
	problems = gsub("eqtl_sumstats", "", coloc_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	easy_lineup = match(problems, easy_answer_key$test_case)
	easy_sub_problems = problems[!is.na(easy_lineup)]
	easy_h4pp = coloc_results[which(problems %in% easy_sub_problems),]$clpp_h4
	easy_answers = easy_answer_key$is_coloc[easy_lineup[!is.na(easy_lineup)]]
	easiest_lineup = match(problems, easiest_answer_key$test_case)
	easiest_sub_problems = problems[!is.na(easiest_lineup)]
	easiest_h4pp = coloc_results[which(problems %in% easiest_sub_problems),]$clpp_h4
	easiest_answers = easiest_answer_key$is_coloc[easiest_lineup[!is.na(easiest_lineup)]]

	# Make COLOC ROC
	roc(easy_answers, easy_h4pp, plot=TRUE)
	roc(easiest_answers, easiest_h4pp, plot=TRUE)
	

	# What if we change the min strength of variant that's allowed to be
	# considered "truly" causal in GWAS and/or eQTL study?

	# What if we require at least a certain sample size for the GWAS/eQTL studies?

	# TODO: How can we adjust these parameters just right so that they simulate an actual GWAS trait
	# with a certain genetic architecture?
	# - Could perhaps say the genome contains 3000 (or however many) possible regions to simulate...

	# (Other adjustments might have to be done at the simulation and/or colocalization test stage.)

	# Putting these results together, for all the methods...could be a main Figure for the paper.
	# (would do it for every method, but could select a representative method to be the face in the
	# main text). Each of these could be a multi-paneled figure showing how performance changes in
	# response to various SIMULATION-based choices (not pertaining to parameters of the methods themselves)

	# What if we try some ensemble, or even just a weighted classifier, of the different methods?

}

main()

# TODO now:
# - Implement baseline method: max (-log eQTL + -log GWAS pval) across all SNPs at the site...
# - Replot with more test cases
# - Penalize sites where the eQTL/GWAS p-value aren't high enough (set them to 0 or lower or something)
# - Figure out what to do with sites where one of the methods fails running (right now, we throw them away)
# - Compare each method when at its best (optimizing params)

