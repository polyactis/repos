function verbose(verbosity,varargin)

%VERBOSE just a fprintf, that is printed when verbosity is >0
% verbose(verbosity,varargin)
% This is useful when we desire different levels of debugging messages. 

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


if verbosity>0
    fprintf(varargin{:})
end
