clear;clc
addpath data
addpath results



importdata('schaefer100_ThickAvg.csv')
% importdata('schaefer100_GrayVol.csv')

importdata('sex.csv');
importdata('group.csv');
importdata('age.csv');
importdata('ICV.csv');

schaefer100_CT=readmatrix('schaefer100_ThickAvg.csv','Range',[2 2]);
% schaefer100_CT=readmatrix('schaefer100_GrayVol.csv','Range',[2 2]);



% 打开文件
fid = fopen('schaefer100_GrayVol.csv');

% 读取第一行
firstLine = fgetl(fid);

% 将字符串分割为单独的列名
columnNames = strsplit(firstLine, ',');
columnNames = columnNames(2:end)';
% 关闭文件
fclose(fid);
filepath_coord = 'F:\first_imag_trans\hansen_receptors\data\schaefer\coordinates\Schaefer_100_centres.txt';
coordinate_data = readtable(filepath_coord);

% 假设 coordinate_data 是一个表格
coordinate_data.Var1 = columnNames;



ICV=readmatrix('ICV.csv','Range',[2 2]);
sex=readmatrix('sex.csv','Range',[2 2]);
group=readmatrix('group.csv','Range',[2 2]);
age=readmatrix('age.csv','Range',[2 2]);

clear ans



% number of subjects- 151 for Maastricht, 115 for Dublin and 146 for Cobre
nregs=100;
% number of regions
nsubs=length(group); 
nregs_lh=50;
% ## Calculate the morphometric similarity matrices:

% z-score the inputs:
schaefer100_CT_zscore=zscore(transpose(schaefer100_CT));



% ## Global differences in morphometric similarity: 


x1=age;
x2=sex;
X = [ones(size(x1)) x1 x2 x1.*x2];
meanMS_regional=schaefer100_CT_zscore';
% Calculate regional residuals:
clear myresid_region
for region=1:nregs
    y = meanMS_regional(:,region);
    [b,bint,resid]=regress(y,X);
    YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
    myresid_region(region,:)=y-YFIT;
end


pats=find(group==2);
cons=find(group==1);

x = reshape(myresid_region(:,cons),[nregs*length(cons),1]);
y = reshape(myresid_region(:,pats),[nregs*length(pats),1]);

% Plot histograms of regional residuals for control subjects and patients:
figure
h1 = histogram(x);
hold on
h2 = histogram(y);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
legend('Controls','Patients')
xlim([-0.1 0.1])
xlabel('Cortical thickness regional residuals')
ylabel('Relative frequency')

% Calculate mean MS: 算出每个被试在每个脑区的总平均值
for subj=1:nsubs
    meanMS(subj)=mean(meanMS_regional(subj,:));
end

tbl = table(age,sex,group,transpose(meanMS));
tbl.sex = categorical(tbl.sex);
tbl.group = categorical(tbl.group);
lm = fitlm(tbl,'Var4~age*sex+group');
p_mean=lm.Coefficients{4,4} % p-value for the effect of group on mean MS

% Plot box plot:

for subj=1:nsubs
    myresid_region_mean(subj)=mean(myresid_region(:,subj));
end

figure
boxplot(myresid_region_mean,group,'notch','on','Labels',{'Controls','Patients'})
ylim([-3*10^(-3) 6*10^(-3)])
ylabel('Mean residual')



% ## Regional differences in morphometric similarity:

% The following code calculates a t-statistic for regional differences in morphometric similarity:


dummy=meanMS_regional;

clear mytstat mypval
for region=1:nregs
  tbl = table(age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~age*sex+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end



mytstat=transpose(mytstat);
mypval=transpose(mypval);
sigregs=find(mypval<0.05)
pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values 如何计算出总体的P值？
sigregs=find(pvalue_fdr<0.05)

dlmwrite('results/schaefer100_mytstat.dat',mytstat)
dlmwrite('results/schaefer100_mypval.dat',mypval)


meantstat = mytstat
meanMS_regional_con=meanMS_regional(find(group==2),:);
meanMS_con = mean(meanMS_regional_con,1)';

meanMS_regional_patient=meanMS_regional(find(group==1),:);
meanMS_patient= mean(meanMS_regional_patient,1)';

% Create a copy of mytstat for modification
mytstat_sig = mytstat;

% Find elements in mytstat_sig where corresponding mypval >= 0.05 and set them to 0
mytstat_sig(mypval >= 0.05) = 0;
% Create table
T = table(mytstat, mytstat_sig, meanMS_con, meanMS_patient);

% Change Variable Names
T.Properties.VariableNames{'mytstat'} = 'mytstat';
T.Properties.VariableNames{'mytstat_sig'} = 'mytstat_sig';
T.Properties.VariableNames{'meanMS_con'} = 'meanMS_con';
T.Properties.VariableNames{'meanMS_patient'} = 'meanMS_patient';

% Write table to file
% writetable(T, 'results/schaefer100_zq.csv');

writetable(T, 'results/schaefer100_zq_post.csv');


% ```
figure
scatter(meanMS_con,meantstat,'x')
refline(0,0)
hold on
plot([0 0], ylim)
plot([0 0], ylim,'-b')
xlabel('Mean control MS')
ylabel('Mean t-statistic')

% Calculate the percentage of scatter points in each quadrant:

a=0;
b=0;
c=0;
d=0;

xvalues=meanMS_con;
yvalues=meantstat;

% clear a b c d
for ind=1:nregs
    xval=xvalues(ind);
    yval=yvalues(ind);
    if ((xval<0)&&yval>0) % a is top left quadrant
        a=a+1;
        elseif ((xval>0)&&yval>0) % b is top right quadrant
        b=b+1;
        elseif ((xval<0)&&yval<0) % c is bottom left quadrant
        c=c+1;
        elseif ((xval>0)&&yval<0) % d is bottom right quadrant
        d=d+1;
    end
end

a=a/nregs % percentage of scatter points in the top left quadrant
b=b/nregs % percentage of scatter points in the top right quadrant
c=c/nregs % percentage of scatter points in the bottom left quadrant
d=d/nregs % percentage of scatter points in the bottom right quadrant
% ```
% 
% ## Regional differences in MS- lh only

% For the gene expression analyses, calculate the t-statistic for the lh only:
% 
% ```

for subj=1:nsubs
    meanMS_regional_lh(subj,:)=sum(subj_MSN_5{1,subj}(1:nregs_lh,1:nregs_lh))./(nregs_lh-1);
end

dummy=meanMS_regional_lh;
clear mytstat mypval
for region=1:nregs_lh
  tbl = table(age,sex,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var4~age*sex+group');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end
mytstat=transpose(mytstat);
mypval=transpose(mypval);
dlmwrite('results/schaefer100_mytstat_lh.dat',mytstat)
dlmwrite('results/schaefer100_mypval_lh.dat',mypval)