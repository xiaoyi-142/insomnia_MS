
dim=15;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);


for j=1:1000
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);
end
p=length(find(Rsq>=Rsquared))/j;


