import json
import copy
import subprocess
import time

num_sims = 2500

with open("colocalization-config.json") as f:
	config = json.load(f)

for locus_num in range(1,num_sims+1):

	# Create specific config file
	
	cfg = copy.deepcopy(config)
	studies = list(cfg["studies"])
	for s in studies:
		cfg["studies"][s.format(locus_num)] = cfg["studies"][s]
		cfg["studies"][s.format(locus_num)]["file"] = cfg["studies"][s.format(locus_num)]["file"].format(locus_num)
		del cfg["studies"][s]

	cfg_file = f"tmp/config{locus_num}"
	with open(cfg_file, "w") as w:
		json.dump(cfg, w)

	# Run snakemake pipe then
	#subprocess.run(f"snakemake -np --snakefile /oak/stanford/groups/smontgom/mgloud/projects/snakemake_coloc_pipeline/brain_gwas/snakemake/Snakefile --configfile {cfg_file}".split())

	snakes = subprocess.run("ps -ef".split(), capture_output = True).stdout.decode('utf-8').split("\n")
	snakes = [s for s in snakes if "snakemake" in s and "mgloud" in s]	
	
	while len(snakes) > 15:
		time.sleep(3)	
		snakes = subprocess.run("ps -ef".split(), capture_output = True).stdout.decode('utf-8').split("\n")
		snakes = [s for s in snakes if "snakemake" in s and "mgloud" in s]	

	subprocess.Popen(f"snakemake --cores 1 --snakefile /oak/stanford/groups/smontgom/mgloud/projects/snakemake_coloc_pipeline/brain_gwas/snakemake/Snakefile --configfile {cfg_file}".split())
			
	

