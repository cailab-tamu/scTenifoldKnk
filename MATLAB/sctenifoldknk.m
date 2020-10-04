function T=sctenifoldknk(X,genelist,kogene,varargin)
% T=sctenifoldknk(X,genelist,"Foxp3");
%
% X is a gene x cell matrix from the wild-type

    if nargin<3
        error(sprintf('USAGE: T=sctenifoldknk(X,genelist,kogene);\n       T=sctenifoldnet_m(X0,X1,genelist,''qqplot'',true);'));
    end
    if isscalar(kogene) && isnumeric(kogene)
        idx=kogene;
    else
        idx=find(genelist==kogene,1);
        if isempty(idx)
            error("KOGENE should be a member of GENELIST.");
        end
    end
   p = inputParser;
   addOptional(p,'qqplot',false,@islogical);
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',false,@islogical);   
   parse(p,varargin{:});
   doqqplot=p.Results.qqplot;
   tdmethod=p.Results.tdmethod;
   nsubsmpl=p.Results.nsubsmpl;
   csubsmpl=p.Results.csubsmpl;
   smplmethod=p.Results.smplmethod;
   savegrn=p.Results.savegrn;
   
   switch upper(tdmethod)
       case "CP"
           tdmethod=1;
       case "TUCKER"
           tdmethod=2;
   end
   switch lower(smplmethod)
       case "jackknife"
           usebootstrp=false;
       case "bootstrap"
           usebootstrp=true;
   end
   
    if size(X,1)~=length(genelist)
        error('Length of genelist should be the same as the number of rows of X0 or X1.');
    end
%     pw0=pwd;
%     pw1=fileparts(which(mfilename));
%     cd(pw1);
%     addpath('thirdparty/tensor_toolbox');
%     cd(pw0);
    if exist('i_td1.m','file')~=2
        error('Need i_td1.m in the scTendifoldNet https://github.com/cailab-tamu/scTenifoldNet/tree/master/MATLAB');
    end    
    
    if exist('sc_pcnet.m','file')~=2
        error('Need sc_pcnet.m in the scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    end    
    
    X=sc_norm(X,"type","libsize");    
    
    [XM]=i_nc(X,nsubsmpl,3,csubsmpl,usebootstrp);
    [A0]=i_td1(XM,tdmethod);
    if savegrn
        tstr=matlab.lang.makeValidName(datestr(datetime));
        save(sprintf('A0_%s',tstr),'A0','genelist','-v7.3');        
    end
    A1=A0;
    A1(idx,:)=0;
    [aln0,aln1]=i_ma(A0,A1);
    T=i_dr(aln0,aln1,genelist,doqqplot);
end
