function parstruct = parsepvpairsparstruct(parstruct, varargin)
%PARSEPVPAIRSPARSTRUCT Validate parameter name/value pairs and throw errors if necessary.
%   Given the struct of valid parameter names as fields, a corresponding
%   values as parameter default values, and a variable length list of parameter 
%   name/value pairs, validate the specified name/value pairs and assign values to output
%   parameters struct.
%
%   parstruct = parsepvpairsparstruct(parstruct, 'Name1', Value1, 'Name2', Value2, ...)
%
%   Inputs:
%   parstruct - struct with fields as parameter names, and values the
%               default parameter
%
%     Field# - Character strings of parameter names to be validated and
%              assigned the corresponding value that immediately follows each
%              in the input argument list. Parameter name validation is
%              case-insensitive and partial string matches are allowed provided
%              no ambiguities exist.
%
%     Value# - The values assigned to the corresponding parameter that
%              immediately precede each in the input argument list.
%
%   Outputs:
%         P# - Parameters assigned the parameter values Value1, Value2, ...
%              in the same order as the names listed in Names. Parameters
%              corresponding to entries in Names that are not specified in the
%              name/value pairs are set to the corresponding value listed in
%              Defaults.

% Initialize some variables.
nInputs   = length(varargin);  % # of input arguments

% Ensure parameter/value pairs.
if mod(nInputs, 2) ~= 0
   error('parsepvpairs:incorrectNumberOfInputs', ...
      'Input parameters must be in name/value pair format.');

else
   names = fields(parstruct);
   %names = lower(names);
   defaults=struct2cell(parstruct);
   
   % Process p/v pairs.
   for j = 1:2:nInputs
      pName = varargin{j};

      if ~ischar(pName)
         error('finance:parsepvpairs:nonTextString', ...
            'Parameter names must be character strings.');
      end

      i = strmatch(lower(pName), lower(names),'exact');

      if isempty(i)
         error('finance:parsepvpairs:invalidParameter', ...
            'Invalid parameter name:  %s.', pName);

      elseif length(i) > 1
         % If ambiguities exist, check for exact match to narrow search.
         i = strmatch(lower(pName),  lower(names), 'exact');
         if length(i) == 1
            %varargout{i} = varargin{j+1};
            parstruct.(names{i})=varargin{j+1};
         else
            error('finance:parsepvpairs:ambiguousParameter', ...
               'Ambiguous parameter name:  %s.', pName);
         end

      else
         %varargout{i} = varargin{j+1};
         parstruct.(names{i})=varargin{j+1};
      end
   end
end
