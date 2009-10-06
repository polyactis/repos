function outputGADA(Iext, Wext, output_fname)
% output GADAJRN's return data into a file
%
% Example:
%
% Iext={[0 33 66 99], [...], };
% Wext={[0 +10 -10], [...], };
%
% output format, tab-delimited:
% 	sample-index, probe1-index, probe2-index, length, amplitude

disp('outputting ...');

output_matrix = cell(1,5);	% sample-index, probe1-index, probe2-index, length, amplitude
line_of_output = 0;
N = length(Iext);
fid = fopen(output_fname, 'w');

for n=1:N;
	Iext_sample = Iext{n};
	Wext_sample = Wext{n};
	SegAmp = MexIextWextToSegAmp(Iext{n}, Wext{n});
	for k=1:length(SegAmp);
		output_matrix{1, 1} = n;
		output_matrix{1, 2} = Iext_sample(k)+1;
		output_matrix{1, 3} = Iext_sample(k+1);
		output_matrix{1, 4} = Iext_sample(k+1)-Iext_sample(k);
		output_matrix{1, 5} = SegAmp(k);
		fprintf(fid, '%d\t%d\t%d\t%d\t%f\n', output_matrix{1,:});
		line_of_output = line_of_output + 1;
	end;
	%K = length(Iext_sample)-2;	% Number of breakpoints
	%M = Iext_sample{K+1};	# last probe index
	%TotalMean=Wext_sample{1};
	%SegAmp = zeros(K+1,1);
	%SegAmp(1) = 0;
	%for k=2:K+1;
	%	SegAmp(k)=Wext_sample(k)/sqrt((M - Iext_sample{k})*Iext_sample{k}/M) + SegAmp(k-1);
	%AuxMean=0;
	%for k=1:K+1;
	%	AuxMean=AuxMean+SegAmp(k)*(Iext_sample{k+1}-Iext_sample{k});
	%AuxMean=AuxMean/M;
	%for k=1:K+1;
	%	SegAmp(k)=SegAmp(k)-AuxMean+TotalMean;
end;

fclose(fid);

%csvwrite(output_fname, cell2mat(output_matrix), line_of_output-1, 5);

%Xrec(n,:)=ReconstructFromIextWext(GADA_JRN.Iext{n},GADA_JRN.Wext{n});

%assert( length(Iext) == length(Wext)+1, 'Dimension mismatch');

%xRec=MexReconstructFromIextWext(int32(Iext),Wext);

disp(sprintf('%d samples. %d lines. %f segments per sample.', N, line_of_output, line_of_output/N));
