function xRec=ReconstructFromIextWext(Iext,Wext)
% ReconstructFromIextWext -- Reconstructs from breakpoint weights and positions using PWCnorm 
% xRec=ReconstructFromIextWext(Iext,Wext);
% Iext -- int32 with the breakpoint positions [0,i_1,..., i_K, M];
% Wext -- Wext(j+1) Breakpoint weights for i_j  for j=1...K, with Wext(1)
% being the overall mean.
% xRec -- Recontruction after concatenating the segments. 
%
% Example:
%
% Iext=[0 33 66 99];
% Wext=[0 +10 -10];
% SegLen=diff(Iext)
% SegAmp=IextWextToSegAmp(Iext,Wext)
% xRec=ReconstructFromIextWext(Iext,Wext);
%
% Calls the mex function  MexReconstructFromIextWext.c
%
% See also: IextWextToSegAmp

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


assert( length(Iext) == length(Wext)+1, 'Dimension mismatch');
xRec=MexReconstructFromIextWext(int32(Iext),Wext);

