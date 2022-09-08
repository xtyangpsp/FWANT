function kern_createfilterfile(fband,syndt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab script to create filter files which will be used during the kernel
% calculation  Modified from kern_createfilterfile.m
% Modified by Xiaotao Yang @ UMASS
% Changes: only need to specify frequency bands. 
% the program automatically derive other variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nfb, ~]=size(fband);
N=2; %number of poles in the filter
% -----------------------filters------------
for n=1:nfb
	Feff=1/syndt/2;
	Fc=fband(n,:);
	Wn=Fc/Feff;     
	[b,a]=butter(N,Wn); 
    	NOD1=length(b);
	
	fname=strcat('./butter.bp.',num2str(fband(n,1)),'_',...
			num2str(fband(n,2)),'.matlab.data');
	fid=fopen(fname,'w');
	fprintf(fid,'%i\n',NOD1);
	for m=1:NOD1
		fprintf(fid,'%20.16f %20.16f\n',a(m),b(m));
	end
	fclose(fid);
end	
end