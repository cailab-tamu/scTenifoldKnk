NCELLS=2000;
NGENES=400;

X0=nbinrnd(20,0.98,NGENES,NCELLS);
X0=X0(:,sum(X0)>120);

genelist=strings(NGENES,1);
for k=1:NGENES
    genelist(k)=sprintf("g%d",k);
end

T=sctenifoldknk(X0,genelist,"g1",...
                'qqplot',true);

