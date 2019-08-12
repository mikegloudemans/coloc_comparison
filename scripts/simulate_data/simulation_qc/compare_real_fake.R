require(readr)
require(dplyr)

options(readr.num_columns = 0)

# Show real and fake coloc candidates side by side (random positioning).
# See if I can judge which one is which. (Might of course also depend on whether
# actually colocalized...)

min_eqtl_pval = 1e-5
min_gwas_pval = 1e-5

real_padding = 500000	# How far to look on either side of the real site's lead SNP

real_sites_file = "/users/mgloud/projects/gwas/output/snps_to_test.txt"
sim_sites_dir = "/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15"

real_sites = suppressWarnings(read_delim(real_sites_file, delim="\t"))


# Examine examples until we say to stop
total = 0
total_correct = 0
while (TRUE)
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
		break
	}

	# Now load real site
	real_index = sample(1:(dim(real_sites)[1]), 1)
	real_site = real_sites[real_index,]

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

	plot_layout = c(1,2,3,4)
	layout(matrix(c(1,2,3,4),nrow=2))

	# Give real and sim an equal chance of being on the left side of the plot grid
	swap_position = runif(1) < 0.5
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
		plot(sim_eqtl$snp_pos, -log10(sim_eqtl$pvalue), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch)
		plot(sim_gwas$snp_pos, -log10(sim_gwas$pvalue), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch)
		plot(combo$snp_pos, -log10(combo$pvalue_eqtl), xlab = "Position", ylab="eQTL -log10(pvalue)", pch=my_pch)
		plot(combo$snp_pos, -log10(combo$pvalue_gwas), xlab = "Position", ylab="GWAS -log10(pvalue)", pch=my_pch)
		message = "The right plots were from the real datasets."
	}
		

	input = toupper(readline(prompt="Which is real data? Enter (L)eft or (R)ight:"))[1]
	while ((input != "L") && (input != "R"))
	{	
		input = toupper(readline(prompt="Not a valid input. Enter (L)eft or (R)ight:"))
	}
	correct = ((input == "L") == swap_position)

	total = total + 1
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

