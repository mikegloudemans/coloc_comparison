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
	timestamp = "2018-07-24_18-07-19"
	ecaviar_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/ecaviar-comparisons/", timestamp)
	coloc_base_dir = paste0("/users/mgloud/projects/brain_gwas/output/coloc-comparisons/2018-07-24_18-07-19")
	answer_file = paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/", timestamp, "/answer_key.txt")


	# eCAVIAR/FINEMAP (should actually run both of them, because why not)
	# Get CLPP scores and actual answers for each test
	ecaviar_results = get_ecaviar_results(ecaviar_base_dir) 
	clpp = ecaviar_results$clpp
	clpp_mod = ecaviar_results$clpp_mod
	problems = gsub("gwas_sumstats", "", ecaviar_results$base_gwas_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answer_key = get_answer_key(answer_file)
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
	
	# Make eCAVIAR ROC
	roc(answers, clpp, plot=TRUE)
	roc(answers, clpp_mod, plot=TRUE)

	# Now do it with coloc
	coloc_results = get_coloc_results(coloc_base_dir)
	h4pp = coloc_results$clpp_h4
	problems = gsub("eqtl_sumstats", "", coloc_results$eqtl_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

	# Make COLOC ROC
	roc(answers, h4pp, plot=TRUE)
	
	# Now let's try it with a few other adjustments (will look better once
	# we have more sites though...)

	# What if we penalize scores with insignificant p-values?
	# Set FINEMAP to a very small value if log_eqtl_pval or log_gwas_pval
	# is below certain thresholds
	# - TODO: Plot: what happens to CLPP score as we go from no cutoff to imposing a gradually increasing cutoff on the scores?
	# - TODO: Plot: what happens to CLPP_mod score?
	# - TODO: Is this effect sensitive to the simulation parameters for generating test cases?
	# - TODO: What about if we apply the same cutoffs with other methods? Coloc, RTC, SMR...?

	for (cutoff in 0:5)
	{
		gwas_cutoff = cutoff
		eqtl_cutoff = cutoff

		dup_ecaviar_results = ecaviar_results
		dup_ecaviar_results$adjusted = pmin(0.0000001 * dup_ecaviar_results$X.log_eqtl_pval, 0.0000001 * dup_ecaviar_results$X.log_gwas_pval)
		dup_ecaviar_results[(dup_ecaviar_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_ecaviar_results$X.log_gwas_pval < gwas_cutoff),]$clpp = dup_ecaviar_results[(dup_ecaviar_results$X.log_eqtl_pval < eqtl_cutoff) | (dup_ecaviar_results$X.log_gwas_pval < gwas_cutoff),]$adjusted
		clpp = dup_ecaviar_results$clpp
		clpp_mod = dup_ecaviar_results$clpp_mod
		problems = gsub("gwas_sumstats", "", dup_ecaviar_results$base_gwas_file)
		problems = as.numeric(gsub("_txt_gz", "", problems))
		answer_key = get_answer_key(answer_file)
		answers = answer_key$is_coloc[match(problems, answer_key$test_case)]
		
		# Make eCAVIAR ROC
		roc(answers, clpp, plot=TRUE, main=paste("CLPP with min GWAS/eQTL p-value", cutoff))
		readline("Press enter:")
		#roc(answers, clpp_mod, plot=TRUE, main=paste("CLPP_mod with min GWAS/eQTL p-value", cutoff))
		#readline("Press enter:")
	}


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
	problems = gsub("gwas_sumstats", "", ecaviar_results$base_gwas_file)
	problems = as.numeric(gsub("_txt_gz", "", problems))
	easy_lineup = match(problems, easy_answer_key$test_case)
	easy_sub_problems = problems[!is.na(easy_lineup)]
	easy_clpp = ecaviar_results[which(problems %in% easy_sub_problems),]$clpp
	easy_answers = easy_answer_key$is_coloc[easy_lineup[!is.na(easy_lineup)]]
	easiest_lineup = match(problems, easiest_answer_key$test_case)
	easiest_sub_problems = problems[!is.na(easiest_lineup)]
	easiest_clpp = ecaviar_results[which(problems %in% easiest_sub_problems),]$clpp
	easiest_answers = easiest_answer_key$is_coloc[easiest_lineup[!is.na(easiest_lineup)]]
	
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

get_ecaviar_results = function(base_dir)
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

main()

# TODO now:
# - Replot with more test cases
# - Stratify by GWAS/eQTL effect sizes and sample sizes
# - Vary the percentage of all input sites that are controls vs. cases -- in a real application there will be many more controls than cases
# - Penalize sites where the eQTL/GWAS p-value aren't high enough (set them to 0 or lower or something)
# - Try it with other methods (coloc and RTC)
# - Compare each method when at its best

