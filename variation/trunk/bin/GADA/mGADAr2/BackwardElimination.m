function [Iext,Wext]=BackwardElimination(Iext,Wext,sigma2,TBackElim,MinSegLen);
% BackwardElimination - Continues with the Backward elimination
% After using fitSingleGADA we can continue the Backward Elimination
% process with this function to continue removing the least likely
% breakpoints. 
%
% Syntax:  
%   [Iext,Wext]=BackwardElimination(Iext,Wext,sigma2,TBackElim,MinSegLen);
% Inputs:
%   Iext - Breakpoint positions
%   Wext - Breakpoint weights
%   sigma2 - Estimated Noise Variance
%   TBackElim - Backward Elimination Critical Value
%   MinSegLen - Minimum Number of Probes within a segment 
% 
% Outputs:
%   Iext - Breakpoint positions
%   Wext - Breakpoint weights
%
%   
% Implemented in mex built-in function MexBEwTandMinLen.c
%  
% See also: fitSingleGADA.m

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


[Iext,Wext]=MexBEwTandMinLen(Iext,Wext,sigma2,TBackElim,MinSegLen);
