conservation_window = 100


# Get sets of top colocalizations for each method
source("/users/mgloud/projects/coloc_comparisons/scripts/real_data/compare_coloc_finemap.R")

# Store single-position conservation and window conservation in a list
phylop_list = list()
phylop_window_list = list()
gerp_list = list()
gerp_window_list = list()
eigen_list = list()
eigen_window_list = list()
eigenpc_list = list()
eigenpc_window_list = list()
method_sets = list()

# An ensemble-ish method...just a test
#upset_matrix[,9] = (rowSums(upset_matrix[,c(1,3,4,6)]) >= 2) * 1
upset_matrix[,9] = (rowSums(upset_matrix[,c(1,3,4,6)]) >= 3) * 1
colnames(upset_matrix)[9] = "basic_ensemble"

# For each method tested
for (m in 1:dim(upset_matrix)[2]) 
{ 
	print(m)
	method = colnames(upset_matrix)[m]

	gwas_loci = upset_matrix_names[upset_matrix[,m] == 1]

	chr = sapply(gwas_loci, function(x) {strsplit(x ,"_")[[1]][1]})
	pos = sapply(gwas_loci, function(x) {strsplit(x ,"_")[[1]][2]})

	study_info = sapply(gwas_loci, function(x) {strsplit(x ,"_txt_gz_ENSG")[[1]][1]})
	study_info = sapply(study_info, function(x) {strsplit(x ,"_txt_gz_")[[1]][2]})

	trait = sapply(study_info, function(x) {strsplit(x ,"_")[[1]][1]})
	gwas_file = sapply(study_info, function(x) 
				{
					s = strsplit(x ,"_")[[1]]
					return(paste(s[2:length(s)], collapse="_"))
				})

	locus_mat = data.frame(list(chr=chr, pos=pos, trait=trait, gwas_file=gwas_file), stringsAsFactors=FALSE)
	method_sets[[m]] = locus_mat

	unique_hits = locus_mat[!duplicated(locus_mat),]
	
	phylop_window_list[[m]] = rep(0, dim(unique_hits)[1])
	phylop_list[[m]] = rep(0, dim(unique_hits)[1])
	gerp_window_list[[m]] = rep(0, dim(unique_hits)[1])
	gerp_list[[m]] = rep(0, dim(unique_hits)[1])
	eigen_window_list[[m]] = rep(0, dim(unique_hits)[1])
	eigen_list[[m]] = rep(0, dim(unique_hits)[1])
	eigenpc_window_list[[m]] = rep(0, dim(unique_hits)[1])
	eigenpc_list[[m]] = rep(0, dim(unique_hits)[1])

	for (l in 1:dim(unique_hits)[1])
	{
		print(l)
		locus = unique_hits[l,]

		# Get conservation in window of interest
		phylop_text_window = system(sprintf("tabix /mnt/lab_data/montgomery/shared/conservation/phylop/phyloP100way/hg19.100way.phyloP100way.bed.gz chr%s:%s-%s", locus$chr, as.numeric(locus$pos)-conservation_window, as.numeric(locus$pos)+conservation_window), intern=TRUE)
		phylop_window = as.numeric(sapply(phylop_text_window, function(x) {strsplit(x ,"\\t")[[1]][4]}))
		phylop_window_list[[m]][l] = mean(phylop_window, na.rm=TRUE)

		phylop_text = system(sprintf("tabix /mnt/lab_data/montgomery/shared/conservation/phylop/phyloP100way/hg19.100way.phyloP100way.bed.gz chr%s:%s-%s", locus$chr, as.numeric(locus$pos), as.numeric(locus$pos)), intern=TRUE)
		phylop_list[[m]][l] = max(as.numeric(sapply(phylop_text, function(x) {strsplit(x ,"\\t")[[1]][4]})), na.rm=TRUE)

		# Get GERP conservation in window of interest
		gerp_text_window = system(sprintf("tabix /mnt/lab_data/montgomery/shared/conservation/gerp/gerp.scores.bed.gz chr%s:%s-%s", locus$chr, as.numeric(locus$pos)-conservation_window, as.numeric(locus$pos)+conservation_window), intern=TRUE)
		gerp_window = as.numeric(sapply(gerp_text_window, function(x) {strsplit(x ,"\\t")[[1]][4]}))
		gerp_window_list[[m]][l] = mean(gerp_window, na.rm=TRUE)

		gerp_text = system(sprintf("tabix /mnt/lab_data/montgomery/shared/conservation/gerp/gerp.scores.bed.gz chr%s:%s-%s", locus$chr, as.numeric(locus$pos), as.numeric(locus$pos)), intern=TRUE)
		gerp_list[[m]][l] = max(as.numeric(sapply(gerp_text, function(x) {strsplit(x ,"\\t")[[1]][4]})), na.rm=TRUE)

		# Get EIGEN score in window of interest
		eigen_text_window = system(sprintf("tabix /mnt/lab_data/montgomery/shared/Eigen/EIGEN/Eigen_hg19_0916_chr%s.tab.bgz %s:%s-%s", locus$chr, locus$chr, as.numeric(locus$pos)-conservation_window, as.numeric(locus$pos)+conservation_window), intern=TRUE)
		eigen_window = as.numeric(sapply(eigen_text_window, function(x) {strsplit(x ,"\\t")[[1]][5]}))
		eigen_window_list[[m]][l] = mean(eigen_window, na.rm=TRUE)

		eigen_text = system(sprintf("tabix /mnt/lab_data/montgomery/shared/Eigen/EIGEN/Eigen_hg19_0916_chr%s.tab.bgz %s:%s-%s", locus$chr, locus$chr, as.numeric(locus$pos), as.numeric(locus$pos)), intern=TRUE)
		eigen_list[[m]][l] = max(as.numeric(sapply(eigen_text, function(x) {strsplit(x ,"\\t")[[1]][5]})), na.rm=TRUE)

		# Get EIGEN-PC score in window of interest
		eigenpc_text_window = system(sprintf("tabix /mnt/lab_data/montgomery/shared/Eigen/EIGEN-PC/Eigen-PC_hg19_0916_chr%s.tab.bgz %s:%s-%s", locus$chr, locus$chr, as.numeric(locus$pos)-conservation_window, as.numeric(locus$pos)+conservation_window), intern=TRUE)
		eigenpc_window = as.numeric(sapply(eigenpc_text_window, function(x) {strsplit(x ,"\\t")[[1]][5]}))
		eigenpc_window_list[[m]][l] = mean(eigenpc_window, na.rm=TRUE)

		eigenpc_text = system(sprintf("tabix /mnt/lab_data/montgomery/shared/Eigen/EIGEN-PC/Eigen-PC_hg19_0916_chr%s.tab.bgz %s:%s-%s", locus$chr, locus$chr, as.numeric(locus$pos), as.numeric(locus$pos)), intern=TRUE)
		eigenpc_list[[m]][l] = max(as.numeric(sapply(eigenpc_text, function(x) {strsplit(x ,"\\t")[[1]][5]})), na.rm=TRUE)
	}
}

# Also do a sort of "ensemble" method:

# TODO

sum(rowSums(upset_matrix[,c(1,3,4,6)]) >= 3)
gwas_loci = upset_matrix_names[]

# Compare distributions
# (maybe convert to box/violin plot eventually)
plot(density(phylop_list[[1]]))
lines(density(phylop_list[[2]]), col="red")
lines(density(phylop_list[[3]]), col="blue")
lines(density(phylop_list[[4]]), col="forestgreen")
lines(density(phylop_list[[5]]), col="darkgoldenrod4")
lines(density(phylop_list[[6]]), col="orange")
lines(density(phylop_list[[7]]), col="violet")
lines(density(phylop_list[[8]]), col="turquoise")
lines(density(phylop_list[[9]]), col="green")

sapply(phylop_list, mean)
sapply(phylop_list, function(x) {sum(x > 0.5) / length(x)})

plot(density(phylop_window_list[[1]]))
lines(density(phylop_window_list[[2]]), col="red")
lines(density(phylop_window_list[[3]]), col="blue")
lines(density(phylop_window_list[[4]]), col="forestgreen")
lines(density(phylop_window_list[[5]]), col="darkgoldenrod4")
lines(density(phylop_window_list[[6]]), col="orange")
lines(density(phylop_window_list[[7]]), col="violet")
lines(density(phylop_window_list[[8]]), col="turquoise")
lines(density(phylop_window_list[[9]]), col="green")

sapply(phylop_window_list, mean)
sapply(phylop_window_list, function(x) {sum(x > 0.5) / length(x)})

# Compare distributions
# (maybe convert to box/violin plot eventually)
plot(density(gerp_list[[1]]))
lines(density(gerp_list[[2]]), col="red")
lines(density(gerp_list[[3]]), col="blue")
lines(density(gerp_list[[4]]), col="forestgreen")
lines(density(gerp_list[[5]]), col="darkgoldenrod4")
lines(density(gerp_list[[6]]), col="orange")
lines(density(gerp_list[[7]]), col="violet")
lines(density(gerp_list[[8]]), col="turquoise")
lines(density(gerp_list[[9]]), col="green")

sapply(gerp_list, mean)
sapply(gerp_list, function(x) {sum(x > 0.5) / length(x)})

plot(density(gerp_window_list[[1]]))
lines(density(gerp_window_list[[2]]), col="red")
lines(density(gerp_window_list[[3]]), col="blue")
lines(density(gerp_window_list[[4]]), col="forestgreen")
lines(density(gerp_window_list[[5]]), col="darkgoldenrod4")
lines(density(gerp_window_list[[6]]), col="orange")
lines(density(gerp_window_list[[7]]), col="violet")
lines(density(gerp_window_list[[8]]), col="turquoise")
lines(density(gerp_window_list[[9]]), col="green")

sapply(gerp_window_list, mean)
sapply(gerp_window_list, function(x) {sum(x > 0.5) / length(x)})

# Compare distributions
# (maybe convert to box/violin plot eventually)
plot(density(eigen_list[[1]]))
lines(density(eigen_list[[2]]), col="red")
lines(density(eigen_list[[3]]), col="blue")
lines(density(eigen_list[[4]]), col="forestgreen")
lines(density(eigen_list[[5]]), col="darkgoldenrod4")
lines(density(eigen_list[[6]]), col="orange")
lines(density(eigen_list[[7]]), col="violet")
lines(density(eigen_list[[8]]), col="turquoise")
lines(density(eigen_list[[9]]), col="green")

sapply(eigen_list, mean)
sapply(eigen_list, function(x) {sum(x > 0.5) / length(x)})

plot(density(eigen_window_list[[1]], na.rm=TRUE))
lines(density(eigen_window_list[[2]], na.rm=TRUE), col="red")
lines(density(eigen_window_list[[3]], na.rm=TRUE), col="blue")
lines(density(eigen_window_list[[4]], na.rm=TRUE), col="forestgreen")
lines(density(eigen_window_list[[5]], na.rm=TRUE), col="darkgoldenrod4")
lines(density(eigen_window_list[[6]], na.rm=TRUE), col="orange")
lines(density(eigen_window_list[[7]], na.rm=TRUE), col="violet")
lines(density(eigen_window_list[[8]], na.rm=TRUE), col="turquoise")
lines(density(eigen_window_list[[9]], na.rm=TRUE), col="green")

sapply(eigen_window_list, mean, na.rm=TRUE)
sapply(eigen_window_list, function(x) {sum(x > 0.5, na.rm=TRUE) / sum(!is.na(x))})


eigenpc_list = lapply(eigenpc_list, function(x){x[is.finite(x)]})
# Compare distributions
# (maybe convert to box/violin plot eventually)
plot(density(eigenpc_list[[1]]))
lines(density(eigenpc_list[[2]]), col="red")
lines(density(eigenpc_list[[3]]), col="blue")
lines(density(eigenpc_list[[4]]), col="forestgreen")
lines(density(eigenpc_list[[5]]), col="darkgoldenrod4")
lines(density(eigenpc_list[[6]]), col="orange")
lines(density(eigenpc_list[[7]]), col="violet")
lines(density(eigenpc_list[[8]]), col="turquoise")
lines(density(eigenpc_list[[9]]), col="green")

sapply(eigenpc_list, mean)
sapply(eigenpc_list, function(x) {sum(x > 0.5) / length(x)})

eigenpc_window_list = lapply(eigenpc_window_list, function(x){x[is.finite(x)]})
plot(density(eigenpc_window_list[[1]]))
lines(density(eigenpc_window_list[[2]]), col="red")
lines(density(eigenpc_window_list[[3]]), col="blue")
lines(density(eigenpc_window_list[[4]]), col="forestgreen")
lines(density(eigenpc_window_list[[5]]), col="darkgoldenrod4")
lines(density(eigenpc_window_list[[6]]), col="orange")
lines(density(eigenpc_window_list[[7]]), col="violet")
lines(density(eigenpc_window_list[[8]]), col="turquoise")
lines(density(eigenpc_window_list[[9]]), col="green")

sapply(eigenpc_window_list, mean)
sapply(eigenpc_window_list, function(x) {sum(x > 0.5) / length(x)})





# Compare distributions

print(colnames(upset_matrix))

# How many unique SNPs are there in each method's set?
sapply(1:8, function(x) {dim(unique(method_sets[[x]][c("chr", "pos")]))[1] / dim(method_sets[[x]])[1]})

# How many unique SNP / GWAS file combos are there (even if there are
# multiple traits covered in that study?)
sapply(1:8, function(x) {dim(unique(method_sets[[x]][c("chr", "pos", "gwas_file")]))[1] / dim(method_sets[[x]])[1]})
