% GADA -- Genome Alteration Detection Analysis Toolbox for Matlab
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com
%
%   In order to start using the toolbox you need to generate the mex files
%   using MakeMex. Then look the two test files for an example. 
%
%
% Base GADA functions:
%   fitGadaModel                  - fitGadaModel - Fits the original GADA model and its extensions to a large set of arrays. 
%   fitSingleGada                 - fitSingleGADA - GADA approach for a single sample.
%   BackwardElimination           - BackwardElimination - Continues with the Backward elimination
%   CollapseAmp                   - CollapseAmp -- Collapses the amplitudes of the non significantly altered segments into a base level. 
%
% Test functions:  
%   testFitGadaModel              - Test the GADA-JRN and GADA-SMN algorithms
%   testFitSingleGada             - Test the single sample GADA algorithm
%
% Functions for dealing with the PWC representation:
%   ReconstructFromIextWext       - ReconstructFromIextWext -- Reconstructs from breakpoint weights and positions using PWCnorm 
%   IextWextToSegAmp              - IextWextToSegAmp -- Obtains Segment amplitudes from breakpoint weights and positions using PWCnorm 
%   CompH                         - CompH - Computes H=inv(F'F) Tridiagonal symetric F unit norm.
%   CompHs1                       - CompHs Computes Hs=inv(Fs'Fs) Tridiagonal symetric
%   CompFdualXb                   - CompFdual -- Computes w=Fdual*b, F columns unit norm.
%
% Functions for operations dealing with Tridiagonal matrices:
%   TriSymGaxpy                   - TriSymGaxpy - y=T*x+y  T Symmetric Tridiagonal matrix
%   ForwardTridSolve              - ForwardTridSolve - Solves Tridiagonal system Tx=b by LDU factorization
%   BackwardTridSolve             - BackwardTridSolve - Solves Tridiagonal system Tx=b by UDL factorization
%   TridSolve                     - TridSolve - Solves Tridiagonal system Tx=b by LDU and UDL factorization
%   TridOfInvTrid                 - TridOfInvTrid - Computes the Tridiagonal band of inv(T)
%
% Other files:
%   MakeMex                       - Makes all the dll for all the Mex*.c files in the folder
%   parsepvpairsparstruct         - Validate parameter name/value pairs and throw errors if necessary.
%   verbose                       - just a fprintf, that is printed when verbosity is >0
%   comparedisc                   - comparedisc -- Compares two discontinuity sets using Willenbrock measure
%   findchange                    - findchange - Finds breakpoint indices, Ir={i: |x(i+1)-x(i)|>eps}
%