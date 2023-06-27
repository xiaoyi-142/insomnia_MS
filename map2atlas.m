% ## von Economo/Yeo networks:
% 
% The code below assesses whether there are differences in MS within particular von Economo classes or Yeo networks. You will need to import the lists of which Yeo network/von Economo class each of the 308 cortical regions belongs to, which can be found in the files 'Yeo_500_overlap.txt' and 'vonEcon_500_overlap.txt'. The Yeo network mapping was performed by [Jakob Seidlitz](https://github.com/jms290) as part of the paper [Váša et al, Cereb Cortex. 2018](https://doi.org/10.1093/cercor/bhx249). The von Economo class mapping was performed by [Konrad Wagstyl](https://github.com/kwagstyl) and [Dr Kirstie Whitaker](https://github.com/kirstiejane) as part of the paper [Whitaker and Vértes, PNAS 2016](https://doi.org/10.1073/pnas.1601745113).
% 
% ```
% networks=load('vonEcon_500_overlap.txt')
networks=load('Yeo_500_overlap.txt'); % or vonEcon500overlap

% To calculate t-statistics and p-values (for the von Economo class and Yeo network Tables in the SI in the paper):

clear myt myp
for class=1:7
    myregions=find(networks==class);
    for subj=1:nsubs
        classreg(subj)=sum(meanMS_regional(subj,myregions));
    end
    tbl = table(age,sex,group,education,transpose(classreg));
    tbl.sex = categorical(tbl.sex);
    tbl.group = categorical(tbl.group);
    lm = fitlm(tbl, 'Var5 ~ age + sex + group + education');
    myt(class)=lm.Coefficients{4,3};
    myp(class)=lm.Coefficients{4,4};
end

myt=transpose(myt);
myp=transpose(myp);
myp_fdr = mafdr(myp,'BHFDR',1); % 
dlmwrite('results/DK308_myt_Yeo.csv',myt)
dlmwrite('results/DK308_myp_Yeo.csv',myp_fdr)
% To plot box plot for a specific network/class:

class=5; % insert which network/class you're interested in here
myregions=find(networks==class);
for subj=1:nsubs
    classreg(subj)=sum(meanMS_regional(subj,myregions));
end

% x1=age;
% x2=sex;
X = [ones(size(x1)) x1 x2 x3];
y=transpose(classreg);
[b,bint,resid]=regress(y,X);
YFIT = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x3;
myresid=y-YFIT;

figure
boxplot(myresid,group,'notch','on','labels',{'Controls','Patients'})