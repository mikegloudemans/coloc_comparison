# Location of config script
test_specs = "real_data_overlap.config"

# First, get all the overlaps of interest using the "overlap" script
# This way, if the exact set of tissues / traits changes, all we need
# to do is rerun this full wrapper script

# Then, given this list, run colocalization analysis
# (the file to use for pulling the SNP list will be clear from
# the config file)

# Then run post-analysis to create all plots for this part of the
# paper (Figure 2?)
