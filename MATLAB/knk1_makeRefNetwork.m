function [A,outfile]=s1_makeRefNetwork(X,genelist)
A=sc_pcnetdenoised(X);
if nargout>1
    tstr=matlab.lang.makeValidName(datestr(datetime));
    outfile=sprintf('A_%s',tstr);
    outfile=fullfile(pwd,[outfile,'.mat']);
    save(outfile,'A','genelist','-v7.3');
end
end
