#!/usr/bin/python

sim_dir=$1
num_tests=$2

python run_rtc_tests.py $1 $2
python run_coloc_tests.py $1 $2
python run_finemap_tests.py $1 $2
python run_caviarbf_tests.py $1 $2
python run_ecaviar_tests.py $1 $2


