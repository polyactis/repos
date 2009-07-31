function y=TriSymGaxpy(t0,t1,x,y)

% TriSymGaxpy - y=T*x+y  T Symmetric Tridiagonal matrix
% y = TriGaxpy(t0,t1,x,y)
% t0 Central diagonal
% t1 Lower/Upper diagonal
% x vector
% y output vector
% See also: TriGaxpy TridSolve

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


if nargin<4
    y=zeros(size(x));
end

N=length(t0);
y=y+t0.*x;
y(1:N-1)=y(1:N-1)+t1.*x(2:N);
y(2:N)=y(2:N)+t1.*x(1:N-1);


