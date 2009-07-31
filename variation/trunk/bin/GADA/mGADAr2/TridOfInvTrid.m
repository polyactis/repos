function [itu,itc,itl] = TridOfInvTrid(tu,tc,tl);
% TridOfInvTrid - Computes the Tridiagonal band of inv(T)
%
% From the UDL,LDU factorization this function computes the tridiagonal
% band of inv(T)
%
% Syntax:  
%   [itl,itc,itu] = TridOfInvTrid(tu,tc,tl);
%
% Inputs:
%   tu - Upper diagonal of the Tridiagonal matrix T   
%   tc - Main  diagonal of the Tridiagonal matrix T   
%   tl - Lower diagonal of the Tridiagonal matrix T   
%
% Outputs:
%   itu - Upper diagonal of the Tridiagonal matrix inv(T)   
%   itc - Main  diagonal of the Tridiagonal matrix inv(T)   
%   itl - Lower diagonal of the Tridiagonal matrix inv(T)
%
% Example:
%
%   N=100;
%   tc=randn(1,N);
%   tu=randn(1,N-1); 
%   [itl,itc,itu] = TridOfInvTrid(tu,tc,tu);
%   d=DiagOfTriXTri(tu,tc,tu,itl,itc,itu);
%   err=sum(abs(d))-N
%     
% See also: TridSolve

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg
% March 2008; 

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


N=length(tc);
assert(length(tu)==N-1,'Dimensions do not match');
assert(length(tl)==N-1,'Dimensions do not match');

b=zeros(size(tc));
[xF,Fu,Fd,Fl] = MexForwardTridSolve(tu,tc,tl,b);
[xB,Bu,Bd,Bl] = MexBackwardTridSolve(tu,tc,tl,b);
[itl,itc,itu] = MexTridOfInvTrid(Fu,Fd,Bl);

% debugging OK
% A=sparse(diag(tu,1)+diag(tc,0)+diag(tu,-1));
% sum(abs(diag(inv(A),1)-itu'))
% sum(abs(diag(inv(A),-1)-itl'))


