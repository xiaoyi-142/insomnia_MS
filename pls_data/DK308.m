
clear;clc
addpath data
addpath pls
addpath results_dk308_lh

load("spin.mat");

response_var_file = 'data\response_dk308_lh.csv';
predictor_var_file ='data\predictor_dk308_r0.2_lh.csv';
output_dir = 'results_dk308_lh\';

PLS_calculate_stats(response_var_file, predictor_var_file, output_dir,spin);

PLS_bootstrap(response_var_file, predictor_var_file, output_dir);

threshold = 3; %权重值绝对值阈值
csv_file = 'PLS2_geneWeights.csv'; %要进行筛选的PLS成分
select_genes(threshold,csv_file,output_dir)

fid1 = 'results_dk308_lh/PLS1_geneWeights.csv'
fid2 = 'data/Candidate_genes_schizophrenia.csv'
fid3 = 'results_dk308_lh/schizophrenia_pls2_stats.csv'
ABSOLUTE = false;
PLS_candidate_genes(fid1,fid2,fid3,ABSOLUTE)




