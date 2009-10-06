function GADAJRNWrap(varargin)

%% 2009-9-27 wrapper to run GADA JRN.
%% Usage: GADAJRNWrap input_fname output_fname aAlpha TBackElim MinSegLen NumberOfProbes NumberOfSamples
%% 
%%  input_fname is the memory map file in binary format (matrix NumberOfProbes X NumberOfSamples).
%%  output format, tab-delimited:
%% 		sample-index, probe1-index, probe2-index, length, amplitude
%%
%% Example:
%% 
%% matlab -nodisplay -nojvm -r "GADAJRNWrap ~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr4_matrix_n100.memmap
%%	~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr4_n100_GADA_output.tsv 0.5 4 5 100 363"

% clear all;	% 2009-9-27 this will clear the varargin, which is not good.

for i=1:nargin
	fprintf(1, ' Argument %d is %s\n', i, varargin{i});
end;

no_of_required_args = 7;

if nargin<no_of_required_args
	disp(sprintf('need %d arguments. now is %d', no_of_required_args, nargin));
	exit;
end;


input_fname = varargin{1};
output_fname = varargin{2};
aAlpha = str2num(varargin{3});
TBackElim = str2num(varargin{4});
MinSegLen = int32(str2num(varargin{5}));
M = double(str2num(varargin{6}));	%% no. of probes. integer type doesn't work in the formatspec in memmapfile. only "double".
N = double(str2num(varargin{7}));	%% no. of samples

% M = 374112	% call method 43. chr1: 374112; chr2: 229253;
% N = 931;

disp('reading input ...');
addpath('~/script/variation/bin/GADA/mGADAr2');

%% Ymat = dlmread(input_fname, '\t');
%% [M,N]= size(Ymat);	%% M is the number of probes. N is the number of samples.

memmap_fname = input_fname; % 'mytemp.tmp';

disp('creating memmap file ...');

% fd=fopen(memmap_fname, 'wb');
% fwrite(fd, Ymat, 'double');
% fclose(fd);

Ymat2=memmapfile(memmap_fname, 'format', {'double', [M N], 'x'});	% the last 'x' means the matrix will be accessed as Ymat2.Data.x

% disp(size(Ymat2.Data.x) );

% aAlpha=0.5;   %Increase this value and the algorithm will be faster miss some breakpoints (default 0.5)
% TBackElim=4;  %Increase this value for less breakpoints ...
% MinSegLen=5;  %Minumum number of probes for calling a segment.
verbose=1;    %Show debugging messages.

CommonParameters={...
    'TBackElim',TBackElim,'MinSegLen',MinSegLen,...
    'verbose',verbose, 'aAlpha', aAlpha};

disp('running  GADA ...');
[GADA_JRN.Iext, GADA_JRN.Wext, GADA_JRN.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha, 'RefAdjust',1, 'RhoAdjust',1, 'TBackElim',3, 'aAlpha',0.2, 'MinSegLen',0,...
    CommonParameters{:});
    %% 2009-10-5 parameters set behind would override the eponymous ones set beforehand.
    %% CommonParameters{:} would override aAlpha, TBackElim, MinSegLen.

%% We don't need Ymat2 more, 
clear Ymat2;

outputGADA(GADA_JRN.Iext, GADA_JRN.Wext, output_fname);

delete *.tmp;

exit();