Author: Jeremy Tien

Key Files (standalone, runnable, produce plots):
- compare_methods_PAIRS.R --> produces the pairs plot (matrix of scatterplots) comparing all methods using REAL data. OUTPUT FILE: /users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/allmethods_pairs.pdf
- compare_methods_UPSET.R --> produces the upset plot comparing all methods using REAL data. OUTPUT FILE: /users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/allmethods_upset.pdf
- eval/eval_rand_forest/eval_rand_forest.py --> produces the correctly formatted ROC for the ensemble method. OUTPUT FILE: /users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/eval/eval_rand_forest/results/comp_ROC_NEW.png

Less Important Files:
- run-ensemble.py --> code was later copied and adapted into the pipeline (brain_gwas on Github); file right now just acts as a standalone that outputs ensemble score for a given input of real data
- preprocess-realdata.py --> not functional on its own; just a collection of useful shell commands for preparing real data before running through ensemble

NOTE: 
- decistree/, multilogreg/, rand_forest/, and svm/ are directories that contain rough implementations and evaluations of the decision tree, multinomial logistic regression, random forest, and support vector machine classification methods, respectively. Classification methods were evaluated on simulated data, metric for evaluation was model ACCURACY.
- Further evaluation for certain classification methods are in the eval/ directory -- plotting ROC, calculating area under ROC, etc.
