# Load all results, concatenate them

# NOTE: TODO:
# - Right now, we skip some test loci simply because the pipeline
#   sees a top SNP somewhere on the edge of the region and anchors on that SNP.
#   For most methods this means the site is skipped.
#   I need to fix it so that we always anchor on the central (original GWAS) SNP.
#   This seed SNP is now being recorded at the time of the simulations

require(pROC)

main = function()
{
	answer_file = paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/", timestamp, "/answer_key.txt")
	answer_key = get_answer_key(answer_file)
	compare_methods(answer_key, timestamp)
}

compare_methods = function(answer_key, timestamp)
{
	timestamp = "2018-07-27_15-23-15"
	finemap_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/finemap-comparisons/", timestamp)
	caviarbf_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/caviarbf-comparisons/", timestamp)
	coloc_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/coloc-comparisons/", timestamp)
	rtc_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/rtc-comparisons/", timestamp)
	twas_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/twas-comparisons/", timestamp)

	#
	# Part 1: Evaluate individual methods on full data set
	#

	#
	# FINEMAP
	#

	# Get CLPP scores and actual answers for each test
	finemap_results = get_finemap_results(finemap_base_dir) 
	clpp = finemap_results$clpp
	clpp_mod = finemap_results$clpp_mod
	problems = gsub("gwas_sumstats", "", finemap_results$base_gwas_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
	
	# Make FINEMAP ROC
	plot(roc(answers, clpp), print.auc = TRUE, col = "black", print.auc.x = 0.2, print.auc.y = 0.32, main = "Colocalization detection performance")
	plot(roc(answers, clpp_mod), print.auc = TRUE, col = "red", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.28)

	#
	# COLOC
	#
	coloc_results = get_coloc_results(coloc_base_dir)
	h4pp = coloc_results$clpp_h4
	problems = gsub("eqtl_sumstats", "", coloc_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

	# Make COLOC ROC
	plot(roc(answers, h4pp), print.auc = TRUE, col = "blue", add=TRUE, print.auc.x = 0.2, print.auc.y = 0.24)

	#
	# RTC
	#
	rtc_results = get_rtc_results(rtc_base_dir)
	rtc = rtc_results$rtc_score
	problems = gsub("eqtl_sumstats", "", rtc_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
	
	# Make RTC ROC
	plot(roc(answers, rtc), print.auc = TRUE, col = "green", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.20)

	#
	# CAVIARBF
	#
	caviarbf_results = get_caviarbf_results(caviarbf_base_dir) 
	clpp = caviarbf_results$clpp
	problems = gsub("eqtl_sumstats", "", caviarbf_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

	# Make CAVIARBF ROC
	plot(roc(answers, clpp), print.auc = TRUE, col = "gray", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.16)

	#
	# TWAS
	#
	twas_results = get_twas_results(twas_base_dir) 
	twas_p = twas_results$twas_log_pval
	problems = gsub("eqtl_sumstats", "", twas_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

	# Make RTC ROC
	plot(roc(answers, twas_p), print.auc = TRUE, col = "orange", add = TRUE, print.auc.x = 0.2, print.auc.y = 0.12)


	# 
	# Part 1.5: Try an ensemble method
	#

	# Get the test sites for which all methods were tried.	
	# NOTE: Technically I should only be comparing on these anyway,
	# unless I can explain why certain sites weren't tested with certain algorithms.
	# Fix that later
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
	twas_results = get_twas_results(twas_base_dir)
	twas_problems = gsub("eqtl_sumstats", "", twas_results$eqtl_file)
	twas_problems = as.numeric(gsub("_txt_gz", "", twas_problems))

	problems = rtc_problems[rtc_problems %in% finemap_problems]
	problems = problems[problems %in% coloc_problems]
	problems = problems[problems %in% bf_problems]
	problems = problems[problems %in% twas_problems]
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

	# For now we're basically just ranking the results of each method and combining 
	# the ranks
	# In the future though, there are ways we can get more nuance than that
	# if one method is very confident
	h4pp = coloc_results$clpp_h4[match(problems, coloc_problems)]
	h4pp = scale(h4pp)
	rtc = rtc_results$rtc_score[match(problems, rtc_problems)]
	rtc = scale(rtc)
	clpp = finemap_results$clpp[match(problems, finemap_problems)]
	clpp = scale(clpp)
	clpp_mod = finemap_results$clpp_mod[match(problems, finemap_problems)]
	clpp_mod = scale(clpp_mod)
	bf_clpp = bf_results$clpp[match(problems, bf_problems)]
	bf_clpp = scale(bf_clpp)
	twas_p = twas_results$twas_log_pval[match(problems, twas_problems)]
	twas_p = scale(twas_p)

	ensemble = h4pp + rtc + clpp + clpp_mod + bf_clpp + twas_p
	ensemble = ensemble / max(ensemble)
	plot(roc(answers, ensemble), print.auc = TRUE, col = "purple", add=TRUE, print.auc.x = 0.2, print.auc.y = 0.08 )
	legend(0.9, 0.3, legend=c("FINEMAP-CLPP", "FINEMAP-CLPP_mod", "COLOC", "RTC", "CAVIARBF", "TWAS", "ensemble"),
	              col=c("black", "red", "blue", "green", "gray", "orange", "purple"), bg="white", lty=1, cex=0.8, lwd=4)
	# Naive ensemble performance is comparable with the best methods

	readline("Press enter:")
	
	# Note of course that any site causing one or more methods to fail completely
	# will not currently appear in this table.
	results_table = data.frame(list(coloc_h4pp=h4pp, rtc=rtc, finemap_clpp=clpp, finemap_clpp_mod=clpp_mod, caviar_bf_clpp=bf_clpp, twas_logp=twas_p))
	# These two plots are illuminating. Goal will be to pick out some of the off-diagonal elements and to figure out
	# why they're scoring differently in the different methods. Particularly for the very different methods.
	pairs(results_table, lower.panel=NULL)
	readline("Press enter:")
	cols = c("black","red")[as.numeric(answers) + 1]
	pairs(apply(results_table, 2, rank), lower.panel=NULL, col=cols, pch=19)
	readline("Press enter:")
	
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





main()

# TODO now:
# - Replot with more test cases
# - Stratify by GWAS/eQTL effect sizes and sample sizes
# - Vary the percentage of all input sites that are controls vs. cases -- in a real application there will be many more controls than cases
# - Penalize sites where the eQTL/GWAS p-value aren't high enough (set them to 0 or lower or something)
# - Try it with other methods (beyond coloc and RTC)
# - Compare each method when at its best

