function [T]=knk2_knockoutTargetGene(A,genelist,targetgene,dosort)
    if nargin<4, dosort=true; end        
    idx=find(genelist==targetgene,1);
    if isempty(idx)
        error("TARGETGENE is not found in GENELIST");
    end
    A1=A;
    A1(idx,:)=0;
    [aln0,aln1]=i_ma(A,A1);
    T=i_dr(aln0,aln1,genelist,false,dosort);
end
