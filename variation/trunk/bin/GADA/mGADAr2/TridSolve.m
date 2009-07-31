function [x,err] = TridSolve(tu,tc,tl,b);
% TridSolve - Solves Tridiagonal system Tx=b by LDU and UDL factorization
%
% Synopsis:
%   TridSolve solves a tridiagonal system by combining two factorizations: 
%
% * The UDL factorization is obtained by Backward Gaussian elimination
% * The LDU factorization is obtained by Forward Gaussian elimination
% 
% The function can also return an error estimate that may be used to
% asses if the solution is stable. Additionally Fu,Bl can be used to
% compute the tridiagonal band of the inverse of T. 
%
% Syntax:  
%   [x,err] = TridSolve(tu,tc,tl,b);
%
% Inputs:
%   tu - Upper diagonal of the Tridiagonal matrix T   
%   tc - Main  diagonal of the Tridiagonal matrix T   
%   tl - Lower diagonal of the Tridiagonal matrix T   
%   b  - Right-side of the equation to solve
% 
% Outputs:
%   x  - Solution of the tridiagonal system.
%   err  - Error between forward and backward ||xB-xF||
%
% Example:
%   N=100;
%   tc=randn(1,N);
%   tu=randn(1,N-1); 
%   x=randn(1,N);
%   b=TriSymGaxpy(tc,tu,x);
%   [xFB,err] = TridSolve(tu,tc,tu,b);
%   xerr=sum(abs(xFB-x))
%   err
%     
% See also: ForwardTridSolve BackwardTridSolve TriSymGaxpy

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg
% March 2008; 

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com

N=numel(b);
assert(numel(tc)==N,'Dimensions do not match');
assert(numel(tu)==N-1,'Dimensions do not match');
assert(numel(tl)==N-1,'Dimensions do not match');

[xF,Fu,Fd,Fl] = MexForwardTridSolve(tu,tc,tl,b);
[xB,Bu,Bd,Bl] = MexBackwardTridSolve(tu,tc,tl,b);

if nargout>1
    err=sum(abs(xF-xB));
end

x=0.5*(xB+xF);


%   Fu  - Upper diagonal of U in T=LDU
%   Fl  - Lower diagonal of L in T=UDL
%   Fd  - Main diagonal of D in T=LDU



