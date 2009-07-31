function [h0,h1]=CompHs1(Ir,M)

% CompHs Computes Hs=inv(Fs'Fs) Tridiagonal symetric
%  [h0,h1]=CompHs1(Ir,M)
% Ir  Breakpoint Position Set
% M  length of the observed vector
% h0 Main diagonal
% h1 Lower/Upper diagonal
%
% This for the normalized basis case!
% 
% Extend this to the Iext notation.

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


Ms=length(Ir);
a=sqrt(Ir(1:Ms-1)./(M-Ir(1:Ms-1)).*(M-Ir(2:Ms))./Ir(2:Ms));
h0=1./(1-a.*a);
h1=-a.*h0;
h0(Ms)=1;
h0(2:Ms)=a.*a.*h0(1:Ms-1)+h0(2:Ms);


