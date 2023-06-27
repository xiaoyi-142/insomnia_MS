clear;clc
addpath data
importdata('CorticalMeasuresENIGMA_ThickAvg.csv')
importdata('CorticalMeasuresENIGMA_SurfAvg.csv')
importdata('CorticalMeasuresENIGMA_MeanCurv.csv')
importdata('CorticalMeasuresENIGMA_GrayVol.csv')
importdata('CorticalMeasuresENIGMA_GausCurv.csv')
importdata('sex.csv')
importdata('group.csv')
importdata('age.csv')
importdata('ICV.csv')

DK_CT=readmatrix('CorticalMeasuresENIGMA_ThickAvg.csv','Range',[2 2])
DK_SA=readmatrix('CorticalMeasuresENIGMA_SurfAvg.csv','Range',[2 2])
DK_MC=readmatrix('CorticalMeasuresENIGMA_MeanCurv.csv','Range',[2 2])
DK_GM=readmatrix('CorticalMeasuresENIGMA_GrayVol.csv','Range',[2 2])
DK_GC=readmatrix('CorticalMeasuresENIGMA_GausCurv.csv','Range',[2 2])
ICV=readmatrix('ICV.csv','Range',[2 2])
sex=readmatrix('sex.csv','Range',[2 2])
group=readmatrix('group.csv','Range',[2 2])
age=readmatrix('age.csv','Range',[2 2])

clear ans



nregs=68; % number of regions
nsubs=length(group); % number of subjects- 151 for Maastricht, 115 for Dublin and 146 for Cobre

% ## Calculate the morphometric similarity matrices:

% z-score the inputs:
DK_CT_zscore=zscore(transpose(DK_CT));
DK_SA_zscore=zscore(transpose(DK_SA));
DK_GM_zscore=zscore(transpose(DK_GM));
DK_MC_zscore=zscore(transpose(DK_MC));
DK_GC_zscore=zscore(transpose(DK_GC));


% Create a cell for each subject with all of the required inputs:
clear subj_features5
for subj=1:nsubs
    subj_features5{1,subj}(:,1)=DK_CT_zscore(:,subj);
    subj_features5{1,subj}(:,2)=DK_SA_zscore(:,subj);
    subj_features5{1,subj}(:,3)=DK_GM_zscore(:,subj);
    subj_features5{1,subj}(:,4)=DK_MC_zscore(:,subj);
    subj_features5{1,subj}(:,5)=DK_GC_zscore(:,subj);
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


x1=sex;
x2=age;
X = [ones(size(x1)) x1 x2 x1.*x2];

% Calculate regional residuals:
clear myresid_region
for region=1:nregs
    y = meanMS_regional(:,region);
    [b,bint,resid]=regress(y,X);
    YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
    myresid_region(region,:)=y-YFIT;
end


pats=find(group==1);
cons=find(group==2);

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
xlabel('MS- regional residuals')
ylabel('Relative frequency')

% Calculate mean MS: 算出每个被试在308个脑区的总平均值
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


dlmwrite('results/mytstat_Maast.dat',mytstat)
dlmwrite('results/mypval_Maast.dat',mypval)

% This code was run to calclulate t-statistics and p-values for all three datasets (mytstat_Maast, mytstat_Dublin, mytstat_Cobre, mypval_Maast, mypval_Dublin and mypval_Cobre). The p-values were then combined using Fisher's method:

clear pcomb
for region=1:nregs
  pcomb(region)=pfast([mypval_Maast(region),mypval_Dublin(region),mypval_Cobre(region)]);
end

pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values 如何计算出总体的P值？
sigregs=find(pvalue_fdr<0.05) % list of the statistically significant regions

meantstat = (mytstat_Maast + mytstat_Dublin + mytstat_Cobre)./3;
```

The code below calculates the mean regional MS for all controls, from all three datasets.

First you need to extract meanMS_regional (as calculated above) for the control subjects only from all three datasets, e.g. for Maastricht:
```
meanMS_regional_Maast_con=meanMS_regional(find(group==1),:);
```

Then simply average them:

```
meanMS_con=mean(vertcat(meanMS_regional_Maast_con,meanMS_regional_Dublin_con,meanMS_regional_Cobre_con),1);
```

Plot the correlation between the regional mean control MS and the mean t-statistic:为什么要做这一步？
```
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

clear a b c d
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
```

## Regional differences in MS- lh only

For the gene expression analyses, calculate the t-statistic for the lh only:

```
nregs_lh=152;
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

dlmwrite('mytstat_Maast_lh.dat',mytstat)
dlmwrite('mypval_Maast_lh.dat',mypval)
```

## von Economo/Yeo networks:

The code below assesses whether there are differences in MS within particular von Economo classes or Yeo networks. You will need to import the lists of which Yeo network/von Economo class each of the 308 cortical regions belongs to, which can be found in the files 'Yeo_500_overlap.txt' and 'vonEcon_500_overlap.txt'. The Yeo network mapping was performed by [Jakob Seidlitz](https://github.com/jms290) as part of the paper [Váša et al, Cereb Cortex. 2018](https://doi.org/10.1093/cercor/bhx249). The von Economo class mapping was performed by [Konrad Wagstyl](https://github.com/kwagstyl) and [Dr Kirstie Whitaker](https://github.com/kirstiejane) as part of the paper [Whitaker and Vértes, PNAS 2016](https://doi.org/10.1073/pnas.1601745113).

```
networks=Yeo_500_overlap.txt; % or vonEcon500overlap

% To calculate t-statistics and p-values (for the von Economo class and Yeo network Tables in the SI in the paper):

clear myt myp
for class=1:7
    myregions=find(networks==class);
    for subj=1:nsubs
        classreg(subj)=sum(meanMS_regional(subj,myregions));
    end

    tbl = table(age,sex,group,transpose(classreg));
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    lm = fitlm(tbl,'Var4~age*sex+group');
    myt(class)=lm.Coefficients{4,3};
    myp(class)=lm.Coefficients{4,4};
end

% To plot box plot for a specific network/class:

class=4; % insert which network/class you're interested in here
myregions=find(networks==class);
for subj=1:nsubs
    classreg(subj)=sum(meanMS_regional(subj,myregions));
end

x1=age;
x2=sex;
X = [ones(size(x1)) x1 x2 x1.*x2];
y=transpose(classreg);
[b,bint,resid]=regress(y,X);
YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x1.*x2;
myresid=y-YFIT;

figure
boxplot(myresid,group,'notch','on','labels',{'Controls','Patients'})
```
