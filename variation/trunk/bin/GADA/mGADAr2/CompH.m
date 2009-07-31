function [h0,h1]=CompH(dim)

% CompH - Computes H=inv(F'F) Tridiagonal symetric F unit norm.
% [h0,h1]=CompH(dim)
% h0 Main diagonal
% h1 Lower/Upper diagonal
% (Normalized case)

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


h0=(dim-1:-1:1).*(1:dim-1)/dim;
h1=-sqrt(h0(2:dim-1).*h0(1:dim-2));
h0=h0*2;

% Obs: F'F inverse is a tridiagonal matrix
% main diagonal calculation
% for i=(1:dim-1)
%     h1(i)=2*((dim-i)*i)/dim;
% end
% first diagonal calculation
% for i=(1:dim-2)
%     h1(i)=((-1)*(((dim-i-1)*(i+1)*(dim-i)*i)^(1/2)))/dim;
% end






