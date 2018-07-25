# Load all results, concatenate them

timestamp = "2018-07-24_18-07-19"

base_dir = paste0("/users/mgloud/projects/brain_gwas/output/ecaviar-comparisons/", timestamp)
answer_file = paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/", timestamp, "/answer_key.txt")

results_dirs = dir(base_dir)
data = list()
for (rd in results_dirs)
{
	subdir = dir(paste(base_dir, rd, sep="/"))
	results_file = subdir[grepl("clpp_status", subdir)] 
	data[[rd]] = read.table(paste(base_dir, rd, results_file, sep="/"), header=TRUE)
}

all_results = do.call(rbind, data)

answer_key = read.table(answer_file, header=TRUE, sep="\t")


# Let's make the first one a simple test...

# Get vector of CLPP scores for each result
clpp = all_results$clpp
problems = gsub("gwas_sumstats", "", all_results$base_gwas_file)
problems = as.numeric(gsub("_txt_gz", "", problems))

# Make ROC
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
answers = answer_key$is_coloc[match(problems, answer_key$test_case)]

require(pROC)
roc(answers, clpp, plot=TRUE)


# TODO now:
# - Replot with more test cases
# - Stratify by GWAS/eQTL effect sizes and sample sizes
# - Vary the percentage of all input sites that are controls vs. cases -- in a real application there will be many more controls than cases
# - Penalize sites where the eQTL/GWAS p-value aren't high enough (set them to 0 or lower or something)
# - Try it with other methods (coloc and RTC)
# - Compare each method when at its best

