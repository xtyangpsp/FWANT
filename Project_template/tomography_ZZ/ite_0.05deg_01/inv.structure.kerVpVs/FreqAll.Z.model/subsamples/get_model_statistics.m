subsetlistfile='subset.list';
%read all subset directories.
subsetlist=readlines(subsetlistfile);
nlines=length(subsetlist);
if length(char(subsetlist(end))) < 1; nlines = nlines - 1; end
disp(['working on [',num2str(nlines),'] subsets']);
outdir='model_statistics';

if ~exist(outdir,'dir')
	mkdir(outdir)
end

resultdir='result.1th';
smooth=8;
damp=8;
outfilebase=strcat('damp',num2str(damp),'.smot',num2str(smooth),'.st0.ev0.lo0.dat');
%try.damp8.smot8.st0.ev0.lo0.dat
mfile=strcat('try.',outfilebase);

%disp(mfile);
mfilesubpath=[resultdir,'/',mfile];
mdata_all=[];
for i = 1:nlines
	mfilepathfull=strcat(subsetlist(i),'/',mfilesubpath);
	if exist(mfilepathfull,'file')
		datain=load(mfilepathfull);
		nelement=length(datain);
		datain=reshape(datain,[1,nelement]);
		mdata_all = [mdata_all; datain];
	end
end
disp(['Data all with the size of: ',num2str(size(mdata_all))]);

mdata_mean=mean(mdata_all,1);

disp(['Data average with the size of: ', num2str(size(mdata_mean))]);

mdata_std=std(mdata_all,1);
mdata_median=median(mdata_all,1);

outfile_mean=strcat(outdir,'/model_mean.',outfilebase);
outfile_median=strcat(outdir,'/model_median.',outfilebase);
outfile_std=strcat(outdir,'/model_std.',outfilebase);

writematrix(mdata_mean,outfile_mean)
writematrix(mdata_median,outfile_median)
writematrix(mdata_std,outfile_std)


