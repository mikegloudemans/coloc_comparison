###############################################
# SMR
###############################################

mkdir -p /users/mgloud/projects/coloc_comparisons/output/smr/plots

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/*/*smr* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/smr/smr_results.txt
done

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/*/*ERROR* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/smr/errors.txt
done

for f in `ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/*/*skipped* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/smr/skipped.txt
done

for f in `ls /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr`; do
	finemap_file=`ls -1 /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/$f/*smr* | head -n 1`
	cat $finemap_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/smr/smr_results.txt

	cat /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/$f/ERROR* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/smr/errors.txt 2> /dev/null
	
	cat /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/$f/skipped* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/smr/skipped.txt 2> /dev/null

	# Also copy the plots to a combined folder
	cp -r /users/mgloud/projects/brain_gwas/output/some-locus-compare-tests-with-smr/$f/plots/* /users/mgloud/projects/coloc_comparisons/output/smr/plots 2> /dev/null
done


