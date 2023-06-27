

clear;clc;
load("spin.mat");
m = 1;
response_var_file = 'data\response_dk308_lh.csv';
predictor_var_file ='data\predictor_dk308_r0.2_lh.csv';
output_dir = 'results_dk308_lh\';
% Import and process MRI response variables
MRI_data = importdata(response_var_file);
ROIname = MRI_data.textdata(2:end, 1);
ResponseVarNames = MRI_data.textdata(1, 2:end);
MRIdata = MRI_data.data;
%import predictor variables
indata=importdata(predictor_var_file);
GENEdata=indata.data;
% GENEdata(1,:)=[];
genes=indata.textdata(1,:)';
genes(1)=[]
geneindex=1:length(genes);
clear indata

%number of bootstrap iterations
bootnum=10;

%DO PLS in 2 dimensions (with 2 components)
X=GENEdata;
Y=zscore(MRIdata);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions IDs and weights in descending order of weight for both
%components
[R1,p1]=corr([XS(:,1),XS(:,2)],MRIdata);
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);



%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap
% disp('  Bootstrapping - could take a while')
% for i=1:bootnum
    Yp=spin(:,m);
%     [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run

    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run
% end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);


%print out results
PLS1_geneWeights = cell(length(genes), 3);
for i=1:length(genes)
    PLS1_geneWeights{i, 1} = PLS1{i};
    PLS1_geneWeights{i, 2} = geneindex1(i);
    PLS1_geneWeights{i, 3} = Z1(i);
end

PLS2_geneWeights = cell(length(genes), 3);
for i=1:length(genes)
    PLS2_geneWeights{i, 1} = PLS2{i};
    PLS2_geneWeights{i, 2} = geneindex2(i);
    PLS2_geneWeights{i, 3} = Z2(i);
end


% Define the threshold and file paths
threshold = 3; 



% Load PLS2 gene weights
pls_result = PLS2_geneWeights;

% Get Z scores from the 3rd column
z_scores = cell2mat(PLS2_geneWeights(:, 3));

% Split into positive and negative loadings
positive_loadings_indices = find(z_scores > 0);
negative_loadings_indices = find(z_scores < 0);


gene_names = PLS2_geneWeights(:, 1);
gene_ids = PLS2_geneWeights(:, 2);

% Process positive and negative loadings and save results
positive_selected_genes = process_loadings(positive_loadings_indices, z_scores, gene_names, gene_ids, threshold);


negative_selected_genes = process_loadings(negative_loadings_indices, z_scores, gene_names, gene_ids, threshold);



positive_selected_genes = table2cell(positive_selected_genes); % replace with your actual file name
negative_selected_genes = table2cell(negative_selected_genes); % replace with your actual file name

T = readtable('gene_entrez_ids.csv');

% Convert gene names to Entrez IDs and filter out genes without Entrez ID
entrezID_positive = zeros(size(positive_selected_genes, 1), 1);
loadings_positive = zeros(size(positive_selected_genes, 1), 1);
for i = 1:size(positive_selected_genes, 1)
    gene = positive_selected_genes{i, 1};
    row = find(strcmp(T{:, 1}, gene));
    if ~isempty(row)
        entrezID_positive(i) = T{row, 2};
        loadings_positive(i) = positive_selected_genes{i, 3};
    end
end
valid_positive = entrezID_positive ~= 0;
entrezID_positive = entrezID_positive(valid_positive);
loadings_positive = loadings_positive(valid_positive);

entrezID_negative = zeros(size(negative_selected_genes, 1), 1);
loadings_negative = zeros(size(negative_selected_genes, 1), 1);
for i = 1:size(negative_selected_genes, 1)
    gene = negative_selected_genes{i, 1};
    row = find(strcmp(T{:, 1}, gene));
    if ~isempty(row)
        entrezID_negative(i) = T{row, 2};
        loadings_negative(i) = negative_selected_genes{i, 3};
    end
end
valid_negative = entrezID_negative ~= 0;
entrezID_negative = entrezID_negative(valid_negative);
loadings_negative = loadings_negative(valid_negative);



% Function to calculate FDR corrected p-values and select genes
function selected_genes = process_loadings(indices, z_scores, gene_names, gene_ids, threshold)
    p_values = 2 * (1 - normcdf(abs(z_scores(indices))));
    adj_p = mafdr(p_values, 'BHFDR', true);
    selected_indices = indices(find(adj_p < 0.05));
    
    % Include only those genes for which absolute weight value is > threshold
    final_indices = selected_indices(abs(z_scores(selected_indices)) > threshold);
    
    selected_genes = table(gene_names(final_indices), gene_ids(final_indices), z_scores(final_indices), 'VariableNames', {'GeneName', 'GeneID', 'Weight'});
end


