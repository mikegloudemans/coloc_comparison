# Concatenate all results to tmp file

################################################
# COLOC and FINEMAP (1 causal variant max)
################################################

mkdir -p /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/plots

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc/*/*clpp* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_results.txt
done

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc/*/*h4* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/coloc_results.txt
done

for f in `ls /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc`; do
	finemap_file=`ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc/$f/*clpp* | head -n 1`
	coloc_file=`ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc/$f/*h4* | head -n 1`
	cat $finemap_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_results.txt
	cat $coloc_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/coloc_results.txt

	# Also copy the plots to a combined folder
	cp -r /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-coloc/$f/plots/* /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/plots 2> /dev/null
done

################################################
# FINEMAP (2 causal variants max)
################################################

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/finemap-multiple-causal/*/*clpp* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_c2_results.txt
done

mkdir -p /users/mgloud/projects/coloc_comparisons/output/finemap-multiple-causal/plots
for f in `ls /users/mgloud/projects/brain_gwas/output/finemap-multiple-causal`; do
	finemap_file=`ls -1 /users/mgloud/projects/brain_gwas/output/finemap-multiple-causal/$f/*clpp* | head -n 1`
	cat $finemap_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/coloc_vs_finemap/finemap_c2_results.txt

	# Also copy the plots to a combined folder
	cp -r /users/mgloud/projects/brain_gwas/output/finemap-multiple-causal/$f/plots/* /users/mgloud/projects/coloc_comparisons/output/finemap-multiple-causal/plots 2> /dev/null
done

###############################################
# Baseline
###############################################

mkdir -p /users/mgloud/projects/coloc_comparisons/output/baseline/plots

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/*/*baseline* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/baseline/baseline_results.txt
done

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/*/*ERROR* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/baseline/errors.txt
done

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/*/*skipped* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/baseline/skipped.txt
done

for f in `ls /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests`; do
	finemap_file=`ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/$f/*baseline* | head -n 1`
	cat $finemap_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/baseline/baseline_results.txt

	cat /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/$f/ERROR* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/baseline/errors.txt 2> /dev/null
	
	cat /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/$f/skipped* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/baseline/skipped.txt 2> /dev/null

	# Also copy the plots to a combined folder
	cp -r /users/mgloud/projects/brain_gwas/output/some-locus-compare-baseline-tests/$f/plots/* /users/mgloud/projects/coloc_comparisons/output/seline/plots 2> /dev/null
done


