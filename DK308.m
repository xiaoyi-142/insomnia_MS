clear;clc;close all;
addpath data;
addpath results;


% importdata('DK308_ThickAvg.csv');
% importdata('DK308_SurfAvg.csv');
% importdata('DK308_MeanCurv.csv');
% importdata('DK308_GrayVol.csv');
% importdata('DK308_GausCurv.csv');
% importdata('sex.csv');
% importdata('group.csv');
% importdata('age.csv');
% importdata('ICV.csv');
% 
% DK308_CT=readmatrix('DK308_ThickAvg.csv','Range',[2 2]);
% DK308_SA=readmatrix('DK308_SurfAvg.csv','Range',[2 2]);
% DK308_MC=readmatrix('DK308_MeanCurv.csv','Range',[2 2]);
% DK308_GM=readmatrix('DK308_GrayVol.csv','Range',[2 2]);
% DK308_GC=readmatrix('DK308_GausCurv.csv','Range',[2 2]);
% ICV=readmatrix('ICV.csv','Range',[2 2]);
% sex=readmatrix('sex.csv','Range',[2 2]);
% group=readmatrix('group.csv','Range',[2 2]);
% age=readmatrix('age.csv','Range',[2 2]);


% 读取Excel文件中的不同工作表
DK308_CT = readmatrix('MS_orignal_data2.xlsx', 'Sheet', 'DK308_ThickAvg', 'Range', [2 2]);
DK308_SA = readmatrix('MS_orignal_data2.xlsx', 'Sheet', 'DK308_SurfAvg', 'Range', [2 2]);
DK308_MC = readmatrix('MS_orignal_data2.xlsx', 'Sheet', 'DK308_MeanCurv', 'Range', [2 2]);
DK308_GM = readmatrix('MS_orignal_data2.xlsx', 'Sheet', 'DK308_GrayVol', 'Range', [2 2]);
DK308_GC = readmatrix('MS_orignal_data2.xlsx', 'Sheet', 'DK308_GausCurv', 'Range', [2 2]);
% 只取前209行
% DK308_CT = DK308_CT(1:206,:);
% DK308_SA = DK308_SA(1:206,:);
% DK308_MC = DK308_MC(1:206,:);
% DK308_GM = DK308_GM(1:206,:);
% DK308_GC = DK308_GC(1:206,:);

% 从最后一张工作表中读取age、sex、education、ICV、group等变量
info_table = readtable('MS_orignal_data2.xlsx', 'Sheet', 'Demography'); % 这里你需要替换 'InfoSheet' 为你的实际工作表名字

% 从表格中提取单独的列作为变量
age = info_table.age;
sex = info_table.sex;
education = info_table.education;
ICV = info_table.ICV;
group = info_table.group;

% age = info_table.age(1:206);
% sex = info_table.sex(1:206);
% education = info_table.education(1:206);
% ICV = info_table.ICV(1:206);
% group = info_table.group(1:206);

clear ans



% number of subjects- 151 for Maastricht, 115 for Dublin and 146 for Cobre
nregs=308;
% number of regions
nsubs=length(group); 
nregs_lh=152;
% ## Calculate the morphometric similarity matrices:

% z-score the inputs:
DK308_CT_zscore=zscore(transpose(DK308_CT));
DK308_SA_zscore=zscore(transpose(DK308_SA));
DK308_GM_zscore=zscore(transpose(DK308_GM));
DK308_MC_zscore=zscore(transpose(DK308_MC));
DK308_GC_zscore=zscore(transpose(DK308_GC));


% Create a cell for each subject with all of the required inputs:
clear subj_features5
for subj=1:nsubs
    subj_features5{1,subj}(:,1)=DK308_CT_zscore(:,subj);
    subj_features5{1,subj}(:,2)=DK308_SA_zscore(:,subj);
    subj_features5{1,subj}(:,3)=DK308_GM_zscore(:,subj);
    subj_features5{1,subj}(:,4)=DK308_MC_zscore(:,subj);
    subj_features5{1,subj}(:,5)=DK308_GC_zscore(:,subj);
end

% Calculate the MS matrices by correlating all inputs and set the diagonal to zero:
for subj=1:nsubs
    subj_MSN_5{1,subj}=corr(transpose(subj_features5{1,subj}));
    subj_MSN_5{1,subj}(logical(eye(size(subj_MSN_5{1,subj})))) = 0;
end
%算出每个脑区和其他脑区相关的平均值
clear meanMS_regional
for subj=1:nsubs
    meanMS_regional(subj,:)=sum(subj_MSN_5{1,subj})./(nregs-1);
end


% ## Global differences in morphometric similarity: 
x1=age;
x2=sex;
x3=education;
x4= ICV;
X = [ones(size(x1)) x1 x2 x3];
% X = [ones(size(x1)) x1 x2 x3 x4];

% Calculate regional residuals:
clear myresid_region
for region=1:nregs
    y = meanMS_regional(:,region);
    [b,bint,resid]=regress(y,X);
    YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x3;
%     YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x3 + b(5)*x3;
    myresid_region(region,:)=y-YFIT;
end


pats=find(group==2);
cons=find(group==1);
% pats=find(group==1);
% cons=find(group==2);

x = reshape(myresid_region(:,cons),[nregs*length(cons),1]);
y = reshape(myresid_region(:,pats),[nregs*length(pats),1]);


color_matrix = [58/256,101/256,198/256
                56/256,200/256,145/256];   
% Plot histograms of regional residuals for control subjects and patients:
figure
h1 = histogram(x);
hold on
h2 = histogram(y);
h1.Normalization = 'probability';
h1.BinWidth = 0.01;
h2.Normalization = 'probability';
h2.BinWidth = 0.01;
legend('Controls','CID','FontName','Arial','FontSize',24)
xlim([-0.14 0.14])
set(gca,'FontName','Arial','FontSize',24);
% set(h1(1),'FaceColor',color_matrix(1,:));
% set(h2(1),'FaceColor',color_matrix(2,:));
xlabel('Mean MS residuals','FontName','Arial','FontSize',30);
ylabel('Relative frequency','FontName','Arial','FontSize',30);

% Calculate mean MS: 算出每个被试在每个脑区的总平均值
for subj=1:nsubs
    meanMS(subj)=mean(meanMS_regional(subj,:));
end

% 将 education 添加到表格
tbl = table(age, sex, group, education, transpose(meanMS));

% 将 sex 和 group 转换为分类变量
tbl.sex = categorical(tbl.sex);
tbl.group = categorical(tbl.group);

% 在线性模型中添加 education
lm = fitlm(tbl, 'Var5 ~ age + sex + group + education');

% 提取组别的 p 值
p_mean = lm.Coefficients{4, 4}; % 注意这里的下标可能需要根据你的模型的具体情况进行调整


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
  tbl = table(age,sex,group,education,dummy(:,region));
%   tbl = table(age,sex,group,education,ICV,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
%   lm = fitlm(tbl,'Var4~age*sex+group');
%   lm = fitlm(tbl, 'Var6 ~ age + sex + group + education + ICV');
  lm = fitlm(tbl, 'Var5 ~ age + sex + group + education');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

mytstat=transpose(mytstat);
mypval=transpose(mypval);
% sigregs=find(mypval<0.05);
pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values 如何计算出总体的P值？
sigregs=find(pvalue_fdr<0.03);


dlmwrite('results/DK308_mytstat.dat',mytstat)
dlmwrite('results/DK308_mypval.dat',mypval)


meantstat = mytstat;
meanMS_regional_con=meanMS_regional(find(group==1),:);
meanMS_con = mean(meanMS_regional_con,1)';
meanMS_regional_case=meanMS_regional(find(group==2),:);
meanMS_case = mean(meanMS_regional_case,1)';




% % z-sore2plot
% figure
% scatter(zscore(meanMS_con),zscore(meantstat),'x')
% refline(0,0)
% hold on
% plot([0 0], ylim)
% plot([0 0], ylim,'-b')
% xlabel('Mean control MS')
% ylabel('Mean t-statistic')


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

% Plot the correlation between the regional mean control MS and the mean t-statistic:为什么要做这一步？
% ```
figure
scatter(meanMS_con,meantstat,'x')
refline(0,0)
hold on
plot([0 0], ylim)
plot([0 0], ylim,'-b')
xlabel('Mean control MS')
ylabel('Mean t-statistic')


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
  tbl = table(age,sex,group,education,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var5~age+sex+group+education');
  mytstat(region)=lm.Coefficients{4,3};
  mypval(region)=lm.Coefficients{4,4};
end

mytstat=transpose(mytstat);
mypval=transpose(mypval);
% 测试另一半的相关性
% meanMS_regional_con_lh=meanMS_regional_lh(find(group==1),:);
% meanMS_con_lh = mean(meanMS_regional_con_lh,1)';
% figure
% scatter(meanMS_con_lh,mytstat,'x')
% refline(0,0)
% hold on
% plot([0 0], ylim)
% plot([0 0], ylim,'-b')
% xlabel('Mean control MS')
% ylabel('Mean t-statistic')

% 
% 
% dlmwrite('results/DK308_mytstat_lh.dat',mytstat)
% dlmwrite('results/DK308_mypval_lh.dat',mypval)

writematrix(mytstat, 'results/DK308_mytstat_lh.csv');
writematrix(mypval, 'results/DK308_mypval_lh.csv');
