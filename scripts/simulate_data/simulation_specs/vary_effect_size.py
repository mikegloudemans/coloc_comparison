#!/usr/bin/python

import subprocess

config_template_h0 = '''
{{
        "out_dir_group": "vary_effect_size/h0/es_{0}",

        "liftover": "True", 

	"gwas_reference_files": ["/mnt/lab_data/montgomery/mgloud/gwas/data/munged/GWAS_BMI-Height_Yengo_2018.txt.gz"],
	"gwas_threshold": 1e-7,
	"window": 500000,

	"total_test_sites": 100,
        
        "simulate_peer": "False",

        "no_genotypes": "True",

	"p_gwas_causal": 0,
	"p_eqtl_causal_given_gwas_causal": 0,
	"p_eqtl_causal_given_no_gwas_causal": 0,
	"p_same_causal_given_both_causal": 0,

	"gwas_min_effect_size": {0},
	"gwas_max_effect_size": {0},
	"eqtl_min_effect_size": 0.5,
        "eqtl_max_effect_size": 0.5,
	"gwas_min_control_log_sample_size": 4,
	"gwas_max_control_log_sample_size": 4,
	"gwas_min_case_log_sample_size": 4,
	"gwas_max_case_log_sample_size": 4,
	"eqtl_min_log_sample_size": 2.47,
	"eqtl_max_log_sample_size": 2.47
}}'''

config_template_h1 = '''
{{
        "out_dir_group": "vary_effect_size/h1/es_{0}",

        "liftover": "True", 

	"gwas_reference_files": ["/mnt/lab_data/montgomery/mgloud/gwas/data/munged/GWAS_BMI-Height_Yengo_2018.txt.gz"],
	"gwas_threshold": 1e-7,
	"window": 500000,

	"total_test_sites": 100,
        
        "simulate_peer": "False",

        "no_genotypes": "True",

	"p_gwas_causal": 1,
	"p_eqtl_causal_given_gwas_causal": 0,
	"p_eqtl_causal_given_no_gwas_causal": 0,
	"p_same_causal_given_both_causal": 0,

	"gwas_min_effect_size": {0},
	"gwas_max_effect_size": {0},
	"eqtl_min_effect_size": 0.5,
        "eqtl_max_effect_size": 0.5,
	"gwas_min_control_log_sample_size": 4,
	"gwas_max_control_log_sample_size": 4,
	"gwas_min_case_log_sample_size": 4,
	"gwas_max_case_log_sample_size": 4,
	"eqtl_min_log_sample_size": 2.47,
	"eqtl_max_log_sample_size": 2.47
}}'''

config_template_h2 = '''
{{
        "out_dir_group": "vary_effect_size/h2/es_{0}",

        "liftover": "True", 

	"gwas_reference_files": ["/mnt/lab_data/montgomery/mgloud/gwas/data/munged/GWAS_BMI-Height_Yengo_2018.txt.gz"],
	"gwas_threshold": 1e-7,
	"window": 500000,

	"total_test_sites": 100,
        
        "simulate_peer": "False",

        "no_genotypes": "True",

	"p_gwas_causal": 0,
	"p_eqtl_causal_given_gwas_causal": 0,
	"p_eqtl_causal_given_no_gwas_causal": 1,
	"p_same_causal_given_both_causal": 0,

	"gwas_min_effect_size": {0},
	"gwas_max_effect_size": {0},
	"eqtl_min_effect_size": 0.5,
        "eqtl_max_effect_size": 0.5,
	"gwas_min_control_log_sample_size": 4,
	"gwas_max_control_log_sample_size": 4,
	"gwas_min_case_log_sample_size": 4,
	"gwas_max_case_log_sample_size": 4,
	"eqtl_min_log_sample_size": 2.47,
	"eqtl_max_log_sample_size": 2.47
}}'''

config_template_h3 = '''
{{
        "out_dir_group": "vary_effect_size/h3/es_{0}",

        "liftover": "True", 

	"gwas_reference_files": ["/mnt/lab_data/montgomery/mgloud/gwas/data/munged/GWAS_BMI-Height_Yengo_2018.txt.gz"],
	"gwas_threshold": 1e-7,
	"window": 500000,

	"total_test_sites": 100,
        
        "simulate_peer": "False",

        "no_genotypes": "True",

	"p_gwas_causal": 1,
	"p_eqtl_causal_given_gwas_causal": 1,
	"p_eqtl_causal_given_no_gwas_causal": 0,
	"p_same_causal_given_both_causal": 0,

	"gwas_min_effect_size": {0},
	"gwas_max_effect_size": {0},
	"eqtl_min_effect_size": 0.5,
        "eqtl_max_effect_size": 0.5,
	"gwas_min_control_log_sample_size": 4,
	"gwas_max_control_log_sample_size": 4,
	"gwas_min_case_log_sample_size": 4,
	"gwas_max_case_log_sample_size": 4,
	"eqtl_min_log_sample_size": 2.47,
	"eqtl_max_log_sample_size": 2.47
}}'''

config_template_h4 = '''
{{
        "out_dir_group": "vary_effect_size/h4/es_{0}",

        "liftover": "True", 

	"gwas_reference_files": ["/mnt/lab_data/montgomery/mgloud/gwas/data/munged/GWAS_BMI-Height_Yengo_2018.txt.gz"],
	"gwas_threshold": 1e-7,
	"window": 500000,

	"total_test_sites": 100,
        
        "simulate_peer": "False",

        "no_genotypes": "True",

	"p_gwas_causal": 1,
	"p_eqtl_causal_given_gwas_causal": 1,
	"p_eqtl_causal_given_no_gwas_causal": 0,
	"p_same_causal_given_both_causal": 1,

	"gwas_min_effect_size": {0},
	"gwas_max_effect_size": {0},
	"eqtl_min_effect_size": 0.5,
        "eqtl_max_effect_size": 0.5,
	"gwas_min_control_log_sample_size": 4,
	"gwas_max_control_log_sample_size": 4,
	"gwas_min_case_log_sample_size": 4,
	"gwas_max_case_log_sample_size": 4,
	"eqtl_min_log_sample_size": 2.47,
	"eqtl_max_log_sample_size": 2.47
}}'''

for hypothesis in [config_template_h0, config_template_h1, config_template_h2, config_template_h3, config_template_h4]:

    for es in [0.02, 0.05, 0.06, 0.08, 0.18]:

        # Write template to file.
        with open("/users/mgloud/projects/coloc_comparisons/tmp/vary_effect_size.config", "w") as w:
            w.write(hypothesis.format(es))

        # Run simulation

        subprocess.check_call("python /users/mgloud/projects/coloc_comparisons/scripts/simulate_sum_stats.py /users/mgloud/projects/coloc_comparisons/tmp/vary_effect_size.config", shell=True)
