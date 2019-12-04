require(Rlab)

# Assuming 10,000 cases and 10,000 controls, what's our power to detect various odds
# ratios in a GWAS?

num_controls = 10000
num_cases = 10000

num_sims = 100

p_threshold = 5e-8

control_risk_allele_frequency = seq(0.01, 0.99, by=0.01)
gwas_odds_ratio = 1 + exp(seq(-7, 0, by = 0.1))

results = array(0, c(length(control_risk_allele_frequency), length(gwas_odds_ratio)))
colnames(results) = gwas_odds_ratio
rownames(results) = control_risk_allele_frequency

for (i in 1:length(control_risk_allele_frequency))
{
	control_raf = control_risk_allele_frequency[i]
	for (j in 1:length(gwas_odds_ratio))
	{
		gor = gwas_odds_ratio[j]
		success = 0
		for (iter in 1:num_sims)
		{
			# Simulate 10,000 controls
			control_genos = rbern(num_controls, control_raf) + rbern(num_controls, control_raf)
			
			case_raf = control_raf * gor
			case_genos = rbern(num_cases, case_raf) + rbern(num_cases, case_raf)
			
			# Just do a basic chi-squared test
			tbl = rbind(table(case_genos > 0), table(control_genos > 0))
			success = success + (chisq.test(tbl)$p.value < p_threshold)
		}
		print(control_raf)
		print(gor)
		print(success / num_sims)
		print("")
		results[i,j]
	}
}

filled.contour(x = control_risk_allele_frequency, y = log10(gwas_odds_ratio), results,	key.title=title(main = "%\ndetected"))
filled.contour(results, plot.axes = {},	key.title=title(main = "%\ndetected"))

axis(1, labels=control_risk_allele_frequency)
axis(2, labels=gwas_odds_ratio)


