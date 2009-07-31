function [Iext,Wext]=fitSingleGada(y,sigma2,aAlpha,TBackElim,MinSegLen);
% fitSingleGADA - GADA approach for a single sample.
% This function implements the Single GADA model as descrived by the
% original GADA paper [1]. fitSingleGADA encapsulates two algorithms.
% First a Sparse Bayesian Learning model is used to extract the most likely
% copy number breakpoints. Then, a BackwardElimination procedure is applied
% to recursively remove the breakpoints with less evidende. 
%
% Syntax:  
%   [Iext,Wext]=fitSingleGada(y,sigma2,aAlpha,TBackElim,MinSegLen);
% Inputs:
%   y - Observed array intensities ordered by chromosomal position
%   sigma2 - Estimated Noise Variance
%   aAlpha - SBL hyperprior
%   TBackElim - Backward Elimination Critical Value
%   MinSegLen - Minimum Number of Probes within a segment 
% 
% Outputs:
%   Iext - Breakpoint positions
%   Wext - Breakpoint weights
%   
% Implemented in mex built-in function MexSBL_BE.c
%
% References:
% [1] Pique-Regi R, Monso-Varona J,Ortega A, Seeger RC, Triche TJ,
%   Asgharzadeh S: "Sparse representation and Bayesian detection of the
%   genome copy number alterations from microarray data", Bioinformatics ,
%   Feb 2008
%  
% See also: fitGadaModel BackwardElimination

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


[Iext,Wext]=MexSBL_BE(y,sigma2,aAlpha,TBackElim,MinSegLen);

