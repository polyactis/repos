function b=CompFdualXb(b)

% CompFdual -- Computes w=Fdual*b, F columns unit norm.
% b=CompFdualXb(b)
% It does this for normalized case F columns unit norm

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


dim=length(b);
b=diff(b);
b=b.*sqrt((dim-1:-1:1).*(1:dim-1)/dim);







