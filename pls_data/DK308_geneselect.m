% Load PLS2 gene weights
pls_result = readtable('PLS2_geneWeights.csv');

% Get Z scores from the 3rd column
z_scores = pls_result{:, 3};

% Split into positive and negative loadings
positive_loadings_indices = find(z_scores > 0);
negative_loadings_indices = find(z_scores < 0);

% Get gene names and IDs
gene_names = pls_result{:, 1};
gene_ids = pls_result{:, 2};

% Process positive and negative loadings and save results
positive_selected_genes = process_loadings(positive_loadings_indices, z_scores, gene_names, gene_ids);
writetable(positive_selected_genes, fullfile(output_dir, 'positive_selected_genes.csv'));

negative_selected_genes = process_loadings(negative_loadings_indices, z_scores, gene_names, gene_ids);
writetable(negative_selected_genes, fullfile(output_dir, 'negative_selected_genes.csv'));

% Function to calculate FDR corrected p-values and select genes
function selected_genes = process_loadings(indices, z_scores, gene_names, gene_ids)
    p_values = 2 * (1 - normcdf(abs(z_scores(indices))));
    adj_p = mafdr(p_values, 'BHFDR', true);
    selected_indices = indices(find(adj_p < 0.05));
    selected_genes = table(gene_names(selected_indices), gene_ids(selected_indices), z_scores(selected_indices), 'VariableNames', {'GeneName', 'GeneID', 'Weight'});
end
