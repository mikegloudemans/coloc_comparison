# Concatenate all results to tmp file

base="/users/mgloud/projects/brain_gwas/output/coloc_comparisons/2019-07-23_16-15-08.365733"

mkdir -p /users/mgloud/projects/coloc_comparisons/output/full_comparison/plots

for f in `ls -1 $base/*/*twas* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/twas_results.txt
done
for f in `ls -1 $base/*/*finemap* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/finemap_results.txt
done
for f in `ls -1 $base/*/*coloc* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/coloc_results.txt
done
for f in `ls -1 $base/*/*rtc* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/rtc_results.txt
done
for f in `ls -1 $base/*/*_smr_* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/smr_results.txt
done
for f in `ls -1 $base/*/*gsmr* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/gsmr_results.txt
done
for f in `ls -1 $base/*/*baseline* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/baseline_results.txt
done

for f in `ls $base`; do
	finemap_file=`ls -1 $base/$f/*finemap* | head -n 1`
	coloc_file=`ls -1 $base/$f/*coloc* | head -n 1`
	rtc_file=`ls -1 $base/$f/*rtc* | head -n 1`
	smr_file=`ls -1 $base/$f/*_smr_* | head -n 1`
	gsmr_file=`ls -1 $base/$f/*gsmr* | head -n 1`
	twas_file=`ls -1 $base/$f/*twas* | head -n 1`
	baseline_file=`ls -1 $base/$f/*baseline* | head -n 1`
	cat $finemap_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/finemap_results.txt
	cat $coloc_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/coloc_results.txt
	cat $rtc_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/rtc_results.txt
	cat $smr_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/smr_results.txt
	cat $gsmr_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/gsmr_results.txt
	cat $twas_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/twas_results.txt
	cat $baseline_file | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/baseline_results.txt

	# Also copy the plots to a combined folder
	cp -r /users/mgloud/projects/brain_gwas/output/coloc_comparisons/$f/plots/* /users/mgloud/projects/coloc_comparisons/output/full_comparison/plots 2> /dev/null
done


# Save ERRORs and skipped variants too.
for f in `ls -1 $base/*/*ERROR* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/errors.txt
done

for f in `ls -1 $base/*/*skipped* | head -n 1`; do
	cat $f | head -n 1 > /users/mgloud/projects/coloc_comparisons/output/full_comparison/skipped.txt
done

for f in `ls $base`; do
	cat $base/$f/ERROR* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/errors.txt 2> /dev/null
	cat $base/$f/skipped* 2> /dev/null | tail -n +2 >> /users/mgloud/projects/coloc_comparisons/output/full_comparison/skipped.txt 2> /dev/null
done

