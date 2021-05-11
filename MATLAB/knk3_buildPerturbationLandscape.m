function [F]=knk3_buildPerturbationLandscape(A,genelist)

n=length(genelist);
assert(n==size(A,1));
F=zeros(n,n);
    for k=1:n
        k
        [t]=knk2_knockoutTargetGene(A,genelist,genelist(k),false);
        F(:,k)=t.drdist;
    end
end
