
main = function()
{
	answers = get_answer_key()

	answers = get_method_score(answers, "coloc", "COLOC_H4", "COLOC_h4")
	answers = get_method_score(answers, "finemap", "finemap_CLPP", "finemap_clpp")
	answers = get_method_score(answers, "finemap", "finemap_CLPP_mod", "finemap_clpp_mod")
	answers = get_method_score(answers, "baseline", "baseline1", "baseline1", invert=TRUE)
	answers = get_method_score(answers, "baseline", "baseline2", "baseline2", invert=TRUE)
	answers = get_method_score(answers, "baseline", "baseline3", "baseline3", invert=TRUE)
	answers = get_method_score(answers, "baseline", "baseline4", "baseline4", invert=TRUE)

	score_names = c("COLOC_H4", "finemap_CLPP", "finemap_CLPP_mod", "baseline1", "baseline2", "baseline3", "baseline4")

	# Remove rows in which none of the tests ran (probably meaning failed overlap)
	scores = answers[score_names]
	answers = answers[apply(scores, 1, function(x) {!(sum(x!=-1) == 0)}),]
	scores = scores[apply(scores, 1, function(x) {!(sum(x!=-1) == 0)}),]

	# Hopefully there aren't many other ones failing at this point, if any
	print(colSums(scores==-1))

	# But for those that do, fill them with the median (this may still be a little
	# biased for evaluation, but it's better than 0-filling or 1-filling)
	for (sn in score_names)
	{
		med = median(answers[sn][answers[sn] != -1,])
		answers[sn][answers[sn] == -1] = med
	}

	# Now rank scores for each method, for plotting purposes
	#
	# These rankings have no absolute meaning though; they just distinguish sites
	#	relative to one another
	for (sn in score_names)
	{
		vec = rank(answers[sn])
		vec = rank(vec + runif(length(vec), -0.4, 0.4))	# Break ties randomly to avoid annoying plotting artifacts
		answers[sprintf("%s_ranking", sn)] = vec 
	}

	write.table(answers, file="../../output/test_snakemake_on_sim_data/aggregated_scores.tsv", quote=FALSE, sep="\t", row.names=FALSE,col.names=TRUE)

	score_rank_names = paste(score_names, "_ranking", sep="")

	rank_table = answers[score_rank_names]

	print(cor(rank_table))

	colors = rep("black", length(answers$colocalization_status))
	colors[answers$colocalization_status == TRUE] = "red"
	print(colors)
	pairs(rank_table, labels=score_names, pch=20, cex=0.3, col=colors)

}

get_answer_key = function()
{

	answers = read.table("/oak/stanford/groups/smontgom/mgloud/projects/simulate_colocalizations/output/simulations/diverse_simulations/answer_key.txt", header=TRUE,stringsAsFactors=FALSE)

	answers$gwas_causal = sapply(answers$causal_variants, function(x){strsplit(strsplit(x,"gwas:")[[1]][2], ":")[[1]][1]})
	answers$eqtl_causal = sapply(answers$causal_variants, function(x){strsplit(strsplit(x,"eqtl:")[[1]][2], ":")[[1]][1]})

	answers$colocalization_status = (answers$gwas_causal == answers$eqtl_causal)

	answers = answers[c("test_case", "colocalization_status")] 


	return(answers)
}

get_method_score = function(data, method_name, new_score_name, old_score_name, invert=FALSE)
{
	print(sprintf("Getting %s", new_score_name))
	new_answers = data
	new_answers[new_score_name] = unlist(sapply(new_answers$test_case, function(x) 
		{
			if(!file.exists(sprintf("/oak/stanford/groups/smontgom/mgloud/projects/coloc_comparison/output/test_snakemake_on_sim_data/%s/TRAIT1-sim-gwas-%d.TRAIT2-sim-eqtl-%d.METHOD-%s.results.txt",method_name,x,x,method_name)))
			{
				return(-1)			
			}
		
			result = read.table(sprintf("/oak/stanford/groups/smontgom/mgloud/projects/coloc_comparison/output/test_snakemake_on_sim_data/%s/TRAIT1-sim-gwas-%d.TRAIT2-sim-eqtl-%d.METHOD-%s.results.txt",method_name,x,x,method_name), header=TRUE)

			if (dim(result)[1] == 0)
			{
				return(-1)	
			} else
			{
				if (!invert)
				{
					return(result[old_score_name])
				} else
				{
					return(1-result[old_score_name])
				}
			}
		}
	))

	return(new_answers)
}


add_finemap_scores = function(data)
{

}

add_baseline_scores = function(data)
{

}

main()

