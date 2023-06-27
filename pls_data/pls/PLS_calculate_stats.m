

function PLS_calculate_stats(response_var_file, predictor_var_file, output_dir,spin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS calculate stats function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% response_var_file ------ full path to the PLS_MRI_response_vars.csv file
%%%                           that is created by the NSPN_CorticalMyelination
%%%                           wrapper script
%%% predictor_var_file ----- full path to the PLS_gene_predictor_vars.csv file
%%%                           that is provided as raw data
%%% output_dir ------------- where to save the PLS_stats file (for PLS1 and PLS2 together)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Petra Vertes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Re-run PLS to get explained variance and associated stats')


% Import and process MRI response variables
MRI_data = importdata(response_var_file);
ROIname = MRI_data.textdata(2:end, 1);
ResponseVarNames = MRI_data.textdata(1, 2:end);
MRIdata = MRI_data.data;

% Import and process predictor variables
predictor_data = importdata(predictor_var_file);
genes = predictor_data.textdata(1, 2:end)';
GENEdata = predictor_data.data;
geneindex = 1:length(genes);


% %import response variables
% importdata(response_var_file);
% 
% %unwrap and tidy MRI response variable names
% ROIname=ans.textdata(:,1);
% ResponseVarNames=ans.textdata(1,:);
% ResponseVarNames=ans.textdata(1,2);
% ResponseVarNames=ans.textdata(1,:);
% ResponseVarNames(1)=[];
% ROIname(1)=[];
% %and store the response variables in matrix Y
% % MRIdata=ans.data(:,1);
% MRIdata=ans.data
% clear ans
% 
% %import predictor variables
% indata=importdata(predictor_var_file);
% GENEdata=indata.data;
% % GENEdata(1,:)=[];
% genes=indata.textdata(1,:)';
% genes(1)=[]
% geneindex=1:length(genes);
% clear indata


%DO PLS in 2 dimensions (with 2 components)
X=GENEdata;
Y=zscore(MRIdata);
X=zscore(X);
%挑选出自己需要的基因数据
% PLS1+
HTR1A = X(:,4695);
GRIA4 = X(:,4236);
IL4R = X(:,4841);
% PRKCD= X(:,8261);
% DKFZp779M0652= X(:,2705);
% % PLS1-
% UPP1 = X(:,11533);
% GPCPD1 = X(:,4157);
% EPB41 = X(:,3178);
% MIR29B2CHG = X(:,6391);
% URM1 = X(:,11550);


dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);



bar(1:dim, 100*PCTVAR(2,1:dim),'FaceColor',[140/255,0,0],'EdgeColor','none');
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);
grid on


[R1,p1]=corr([XS(:,1),XS(:,2)],MRIdata);
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end


%calculate correlations of PLS components with MRI variables
[R1,p1]=corr(XS(:,1),MRIdata);
[R2,p2]=corr(XS(:,2),MRIdata);
a=[R1',p1',R2',p2'];

% assess significance of PLS result
for j=1:5000
    order=randperm(size(Y,1));
%     order = spin;
    Yp=Y(order,:);
%     Yp=spin(:,j);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);
end
p=length(find(Rsq>=Rsquared))/j;


s1 = zscore(stats.W);

% convert the p-value to a string
p_str = num2str(p);

% create the title with the p-value
title_str = ['pspin = ', p_str];

% plot histogram
hist(Rsq,30)
hold on
plot(Rsquared,300,'.r','MarkerSize',15)
set(gca,'Fontsize',14)
xlabel('R squared','FontSize',14);
ylabel('Permuted runs','FontSize',14);
title(title_str)

%save stats
myStats=[PCTVAR; p, j];
csvwrite(fullfile(output_dir,'PLS_stats.csv'),myStats);
