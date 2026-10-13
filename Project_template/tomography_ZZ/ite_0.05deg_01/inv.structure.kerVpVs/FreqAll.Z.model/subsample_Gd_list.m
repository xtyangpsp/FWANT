%This function get the subsets of the Gd_list file to: 1) avoid operating
%LSQR with large memory, 2) get statistical uncertainties on the inversion,
%3) Get bootstrap-alike reliability test.

nsubset=50; %number of subsets.
subsample_rate=0.7; %rate of subsample 0-1, excluding 0. The nearest integer 
%number of samples will be used.
outdir='subsamples';
if ~exist(outdir,'dir')
    mkdir(outdir)
end
infile='inv_Gd_list';
outdirlist=[outdir,'/subset.list'];

indata=readlines(infile); %read line by line, no need to seperate.
% ../G_spool/ZL.B06X_XO.KH28_BHZ_4_T1T2.P2.Kap.unformat ../G_spool/ZL.B06X_XO.KH28_BHZ_4_T1T2.P2.Kbp.unformat 
% -0.900 0.662252 544 461 461 0 923 0 1385 0 RL

nlines=length(indata);
if length(char(indata(end)))<1; nlines=nlines-1; end % incase the end line is empty.
disp(nlines)
fidout=fopen(outdirlist,'w');
for i=1:nsubset
    disp(['Subset -- ',num2str(i),'/',num2str(nsubset)])
    outsubdir=strcat(outdir,'/',num2str(i));
    if ~exist(outsubdir,'dir')
        mkdir(outsubdir)
    end
    fprintf(fidout,'%s\n',num2str(i));
    outfile=strcat(outsubdir,'/',infile);
    nidx=randsample(nlines,int32(nlines*subsample_rate));
    outdata=indata(nidx);

    writematrix(outdata,outfile)
end
fclose(fidout);
disp(['List of subset directories are saved in ',outdirlist]);
