require(dplyr)
require(readr)
require(reshape2)


data = read_delim("coloc_sample_size_comparison_results.txt", delim="\t")
data = data[data$n_snps >= 50,]

stats = data %>% group_by(ref_snp, eqtl_file, gwas_trait, feature, base_gwas_file) %>% summarize(range = max(clpp_h4) - min(clpp_h4), variance=var(clpp_h4))
hist(stats$range)

# See to what extent throwing away "N=10" runs reduces overall variability of results
sub_data = data[data$eqtl_N > 10 & data$gwas_N > 10,]
sub_stats = sub_data %>% group_by(ref_snp, eqtl_file, gwas_trait, feature, base_gwas_file) %>% summarize(range = max(clpp_h4) - min(clpp_h4), variance=var(clpp_h4))
hist(sub_stats$range)

sub_sub_data = data[data$eqtl_N > 100 & data$gwas_N > 100,]
sub_sub_stats = sub_sub_data %>% group_by(ref_snp, eqtl_file, gwas_trait, feature, base_gwas_file) %>% summarize(range = max(clpp_h4) - min(clpp_h4), variance=var(clpp_h4))
hist(sub_sub_stats$range)


# See to what extent the rankings agree for different specified sample sizes
# Use a gather for this...

casted = dcast(formula = ref_snp + eqtl_file + gwas_trait + feature + base_gwas_file ~ eqtl_N + gwas_N, data = data, value.var = "clpp_h4")

plot(rank(casted$"10_10"), rank(casted$"100_100"))
plot(rank(casted$"100_100"), rank(casted$"1000_1000"))
plot(rank(casted$"1000_1000"), rank(casted$"10000_10000"))

casted$diffs = abs(rank(casted$"1000_1000") - rank(casted$"10000_10000"))
sort_cast = casted[rev(order(casted$diffs)),]

