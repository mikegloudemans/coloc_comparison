require(readr)
require(dplyr)

options(readr.num_columns = 0)

# Show real and fake coloc candidates side by side (random positioning).
# See if I can judge which one is which. (Might of course also depend on whether
# actually colocalized...)

interactive = TRUE

min_eqtl_pval = 1e-7
min_gwas_pval = 1e-7

real_padding = 500000	# How far to look on either side of the real site's lead SNP

real_sites_file = "/users/mgloud/projects/gwas/output/snps_to_test.txt"
sim_sites_dir = "/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15"

real_sites = suppressWarnings(read_delim(real_sites_file, delim="\t"))


# Examine examples until we say to stop
total = 0
total_correct = 0
features_correct = list(num_snps=0, strong_gwas=0, strong_eqtl=0, strong_either=0, complexity=0, likely_coloc=0)

limit = 1000
while (total < limit)
{
	# Load sim site and make sure it passes the cutoffs.
	# If not, keep trying
	while (TRUE)
	{
		sim_index = sample(1:5000, 1)
		sim_eqtl = suppressWarnings(read_delim(sprintf("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/eqtl/eqtl_sumstats%d.txt.gz", sim_index), delim="\t"))
		if (min(sim_eqtl$pvalue) > min_eqtl_pval)
		{
			next
		}
		sim_gwas = suppressWarnings(read_delim(sprintf("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas/gwas_sumstats%d.txt.gz", sim_index), delim="\t"))
		if (min(sim_gwas$pvalue) > min_gwas_pval)
		{
			next
		}

		# By current filters, we need at least one SNP that passes significance in both
		if (sum((min(sim_gwas$pvalue) < min_gwas_pval) & (min(sim_eqtl$pvalue) < min_eqtl_pval)) == 0)
		{
			next
		}
		# Pvalue shouldn't be reported as 0, how does this even happen?
		if (min(sim_gwas$pvalue, sim_eqtl$pvalue) == 0)
		{
			next
		}
		break
	}

	# Now load real site
	while (TRUE)
	{
		real_index = sample(1:(dim(real_sites)[1]), 1)
		real_site = real_sites[real_index,]

		if (is.na(real_site$chr))
		{
			next
		}

		# GWAS data
		header = strsplit(system(sprintf("zcat %s | head -n 1", real_site$gwas_file), intern=TRUE), "\t")[[1]]
		output = paste(system(sprintf("tabix %s %d:%d-%d", real_site$gwas_file, real_site$chr, real_site$snp_pos - real_padding, real_site$snp_pos + real_padding), intern=TRUE), collapse="\n")
		gwas = suppressWarnings(read_delim(output, delim="\t", col_names=header))
		# Filter down to only the trait of interest
		if ("trait" %in% colnames(gwas))
		{
			gwas = gwas[gwas$trait == real_site$trait,]
		}

		# eQTL data
		header = strsplit(system(sprintf("zcat %s | head -n 1", real_site$eqtl_file), intern=TRUE), "\t")[[1]]
		output = paste(system(sprintf("tabix %s %d:%d-%d", real_site$eqtl_file, real_site$chr, real_site$snp_pos - real_padding, real_site$snp_pos + real_padding), intern=TRUE), collapse="\n")
		eqtl = suppressWarnings(read_delim(output, delim="\t", col_names=header))
		# Filter down to only the gene of interest
		eqtl = eqtl[eqtl$gene == real_site$gene,]

		# Merge eQTL and GWAS
		combo = inner_join(eqtl, gwas, by=c("chr", "snp_pos"), suffix=c("_eqtl", "_gwas"))
		
		# Again, I don't know how this happens, but apparently it does
		if (is.na(min(c(combo$pvalue_eqtl, combo$pvalue_gwas))))
		{
			next
		}
		if (min(c(combo$pvalue_gwas, combo$pvalue_eqtl)) == 0)
		{
			next
		}
		break
	}

	# Show plots for the user to compare, unless we've set this to not happen
	if (interactive)
	{
		plot_layout = c(1,2,3,4)
		layout(matrix(c(1,2,3,4),nrow=2))

		# Give real and sim an equal chance of being on the left side of the plot grid
		swap_position = FALSE
		my_pch = 16
		if (swap_position)
		{
			plot(combo$snp_pos, -log10(combo$pvalue_eqtl), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch)
			plot(combo$snp_pos, -log10(combo$pvalue_gwas), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch)
			plot(sim_eqtl$snp_pos, -log10(sim_eqtl$pvalue), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch)
			plot(sim_gwas$snp_pos, -log10(sim_gwas$pvalue), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch)
			message = "The left plots were from the real datasets."
		}
		else
		{
			plot(sim_eqtl$snp_pos, -log10(sim_eqtl$pvalue), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch, ylim=c(0, max(c(10, -log10(sim_eqtl$pvalue), -log10(sim_gwas$pvalue)), na.rm=TRUE)))
			plot(sim_gwas$snp_pos, -log10(sim_gwas$pvalue), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch, ylim=c(0, max(c(10, -log10(sim_eqtl$pvalue), -log10(sim_gwas$pvalue)), na.rm=TRUE)))
			plot(combo$snp_pos, -log10(combo$pvalue_eqtl), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch, ylim=c(0, max(c(10, -log10(combo$pvalue_eqtl), -log10(combo$pvalue_gwas)), na.rm=TRUE)))
			plot(combo$snp_pos, -log10(combo$pvalue_gwas), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch, ylim=c(0, max(c(10, -log10(combo$pvalue_eqtl), -log10(combo$pvalue_gwas)), na.rm=TRUE)))
			message = "The right plots were from the real datasets."
		}
			

		input = toupper(readline(prompt="Which is real data? Enter (L)eft or (R)ight:"))[1]
		while ((input != "L") && (input != "R"))
		{	
			input = toupper(readline(prompt="Not a valid input. Enter (L)eft or (R)ight:"))
		}
		correct = ((input == "L") == swap_position)

		if (correct)
		{
			print("Correct!")
			total_correct = total_correct + 1
		}
		else
		{
			print("Nope.")
		}

		print(message)

		print("Your percentage is:")
		print(total_correct / total)
	}

	# Test a few other metrics to see how predictable it is...
	if (dim(combo)[1] < dim(sim_eqtl)[1])
	{
		features_correct[["num_snps"]] = features_correct[["num_snps"]] + 1
	}
	if (min(c(combo$pvalue_gwas, combo$pvalue_eqtl)) < min(c(sim_gwas$pvalue, sim_eqtl$pvalue)))
	{
		features_correct[["strong_either"]] = features_correct[["strong_either"]] + 1
	}
	if (min(combo$pvalue_gwas) < min(sim_gwas$pvalue))
	{
		features_correct[["strong_gwas"]] = features_correct[["strong_gwas"]] + 1
	}
	if (min(combo$pvalue_eqtl) < min(sim_eqtl$pvalue))
	{
		features_correct[["strong_eqtl"]] = features_correct[["strong_eqtl"]] + 1
	}

	# Complexity = proxy for LD; what fraction of SNPs have p-values not so far off
	# from the lead one?
	complexity_cutoff = 0.75
	real_eqtl_complexity = sum(log(combo$pvalue_eqtl) / log(min(combo$pvalue_eqtl)) > complexity_cutoff) / length(combo$pvalue_eqtl)
	real_gwas_complexity = sum(log(combo$pvalue_gwas) / log(min(combo$pvalue_gwas)) > complexity_cutoff) / length(combo$pvalue_gwas)
	real_avg_complexity = mean(real_eqtl_complexity, real_gwas_complexity)
	sim_eqtl_complexity = sum(log(sim_eqtl$pvalue) / log(min(sim_eqtl$pvalue)) > complexity_cutoff) / length(sim_eqtl$pvalue)
	sim_gwas_complexity = sum(log(sim_gwas$pvalue) / log(min(sim_gwas$pvalue)) > complexity_cutoff) / length(sim_gwas$pvalue)
	sim_avg_complexity = mean(sim_eqtl_complexity, sim_gwas_complexity)

	if (real_avg_complexity > sim_avg_complexity)

	{
		features_correct[["complexity"]] = features_correct[["complexity"]] + 1
	}

	# Do they have same lead variant, or very close?
	coloc_closeness_threshold = 0.95

	# Look up eQTL value for top GWAS variant -- is it near the top variant?
	real_gwas_min = which.min(combo$pvalue_gwas)
	real_coloc_likely = log(combo$pvalue_eqtl[real_gwas_min]) / log(min(combo$pvalue_eqtl)) > coloc_closeness_threshold
	sim_gwas_min = which.min(sim_gwas$pvalue)
	sim_coloc_likely = log(sim_eqtl$pvalue[sim_gwas_min]) / log(min(sim_eqtl$pvalue)) > coloc_closeness_threshold

	if (real_coloc_likely > sim_coloc_likely)
	{
		# Would have guessed real, so correct
		features_correct[["likely_coloc"]] = features_correct[["likely_coloc"]] + 1
	}
	else if (real_coloc_likely < sim_coloc_likely)
	{
		# Would have guessed sim, so incorrect
		# Do nothing
	}
	else
	{
		# Tie; can't choose.
		features_correct[["likely_coloc"]] = features_correct[["likely_coloc"]] + 0.5
	}
	
	total = total+1
	print(total)
	print(unlist(features_correct) / total)
}

