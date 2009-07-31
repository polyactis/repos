% Makes all the dll for all the Mex*.c files in the folder

% This acts as a very simple makefile for building the C libraries necessary
% for GADA

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com

actualcd=pwd;
mypath=fileparts(mfilename('fullpath'));  % Put another path if sources in another path
cd(mypath);

commonfiles='./cgada/BaseGADA.c'; % Base source.
includes={'-I./cgada','-I./'}; % This is where the headers of the base reside. 
options={'-v','-argcheck'} % Set the compiler in high verbosity mode.
defines={'-D_MATLAB_'}; %Activates defines to make BaseGADA functions compile well for Matlab.


ListFiles=struct2cell(dir('Mex*.c')); %Obtain the list of all the files to compile.
ListFiles=ListFiles(1,:)';


for i=1:length(ListFiles)
    currentfile=ListFiles{i};
    fprintf('\n\n***********************************************************\n');
    fprintf('     Compiling %s \n',currentfile);
    fprintf('***********************************************************\n');    
    mexcall=cat(2,options,defines,includes,currentfile,commonfiles)
    mex(mexcall{:});
%     [fpath,fname,fext]=fileparts(currentfile);
%     fname=[fname '.m'];
%     if(~exist(fname,'file'));
%         fd=fopen(fname,'w');
%         fprintf(fd,'%% %s (Mex-file) documentation TODO or in %s \n%%See also: %s ',currentfile,fname(4:end-2),fname(4:end-2));
%         fclose(fd);
%     end        
    fprintf('\n');
    
end

%%
cd(actualcd)