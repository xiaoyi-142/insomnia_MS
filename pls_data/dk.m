
clear;clc
addpath data
addpath pls
addpath results_dk_lh

response_var_file = 'data\response_dk_lh.csv'
predictor_var_file ='data\predictor_dk_r0.2_lh.csv'
output_dir = 'results_dk_lh\'

PLS_calculate_stats(response_var_file, predictor_var_file, output_dir)

PLS_bootstrap(response_var_file, predictor_var_file, output_dir)

fid1 = 'results_dk_lh/PLS1_geneWeights.csv'
fid2 = 'data/Candidate_genes_schizophrenia.csv'
fid3 = 'results_dk_lh/schizophrenia_pls2_stats.csv'
ABSOLUTE = false
PLS_candidate_genes(fid1,fid2,fid3,ABSOLUTE)

p_twoTailed = normcdf(PLS_Z);


