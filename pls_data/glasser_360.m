
clear;clc
addpath data
addpath pls
addpath results_glasser_lh

response_var_file = 'data\response_glasser_lh.csv'
predictor_var_file ='data\predictor_glasser_r0.2_lh.csv'
output_dir = 'results_glasser_lh\'

PLS_calculate_stats(response_var_file, predictor_var_file, output_dir)

PLS_bootstrap(response_var_file, predictor_var_file, output_dir)

fid1 = 'results_glasser_lh/PLS2_geneWeights.csv'
fid2 = 'data/Candidate_genes_schizophrenia.csv'
fid3 = 'results_glasser_lh/schizophrenia_pls2_stats.csv'
ABSOLUTE = false
PLS_candidate_genes(fid1,fid2,fid3,ABSOLUTE)

p_twoTailed = normcdf(PLS_Z);


