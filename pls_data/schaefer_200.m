
clear;clc
addpath data
addpath pls
addpath results_schaefer_200_lh

response_var_file = 'data\response_schaefer_200_lh.csv'
predictor_var_file ='data\predictor_schaefer_200_r0.2_lh.csv'
output_dir = 'results_schaefer_200_lh\'

PLS_calculate_stats(response_var_file, predictor_var_file, output_dir)

PLS_bootstrap(response_var_file, predictor_var_file, output_dir)

fid1 = 'results_schaefer_200_lh/PLS2_geneWeights.csv'
fid2 = 'data/Candidate_genes_schizophrenia.csv'
fid3 = 'results_schaefer_200_lh/schizophrenia_pls2_stats.csv'
ABSOLUTE = false
PLS_candidate_genes(fid1,fid2,fid3,ABSOLUTE)

p_twoTailed = normcdf(PLS_Z);


