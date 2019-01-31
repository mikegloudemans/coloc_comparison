# Concatenate all results to tmp file
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

