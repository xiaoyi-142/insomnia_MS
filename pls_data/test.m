rep=1000;
for dim=1:15
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);
    for j=1:rep
        %j
        order=randperm(size(Y,1));
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq(j) = temp(dim);
    end
dim
R(dim)=Rsquared
p(dim)=length(find(Rsq>=Rsquared))/rep
end
figure
plot(1:dim, p,'ok','MarkerSize',8,'MarkerFaceColor','r');
xlabel('Number of PLS components','FontSize',14);
ylabel('p-value','FontSize',14);
grid on

mypvalue_pls1 = [1.00000,0.99992,1.00000,0.00003,0.00061,1.00000,0.99177]
pvalue_fdr = mafdr(mypvalue_pls1,'BHFDR',1);
