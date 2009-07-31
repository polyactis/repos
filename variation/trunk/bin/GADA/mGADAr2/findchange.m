function Ir=findchange(x,myeps)

% findchange - Finds breakpoint indices, Ir={i: |x(i+1)-x(i)|>eps}
% Ir=findchange(x,eps)
% Ir=find(abs(diff(x))>myeps);

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


if nargin<2
    myeps=0.00001;
end
Ir=find(abs(diff(x))>myeps);
