dummy=meanMS_regional;

clear mytstat mypval
for region=1:nregs
  tbl = table(age,sex,ICV,group,dummy(:,region));
  tbl.sex = categorical(tbl.sex);
  tbl.group = categorical(tbl.group);
  lm = fitlm(tbl,'Var5~age*sex+ICV+group');
  mytstat(region)=lm.Coefficients{5,3};
  mypval(region)=lm.Coefficients{5,4};
end

mytstat=transpose(mytstat);
mypval=transpose(mypval);
sigregs=find(mypval<0.05);
pvalue_fdr = mafdr(mypval,'BHFDR',1); % FDR corrected p-values 如何计算出总体的P值？
sigregs=find(pvalue_fdr<0.05);