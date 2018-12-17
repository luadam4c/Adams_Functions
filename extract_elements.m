function elements = extract_elements (vecs, extractMode, varargin)
%% Extracts elements from vectors using a certain mode ('first', 'last', 'min', 'max')
% Usage: elements = extract_elements (vecs, extractMode, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       elements    - element(s) from each vector extracted
%                   specified as a numeric vector 
%                       or a cell array of column numeric vectors
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array of numeric arrays
%       extractMode - mode of extraction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - first element of each vector
%                       'last'  - last element of each vector
%                       'min'   - minimum-valued element of each vector
%                       'max'   - maximum-valued element of each vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/create_average_time_vector.m
%       cd/plot_protocols.m

% File History:
% 2018-12-15 Created by Adam Lu
% 

%% Hard-coded parameters
validExtractModes = {'first', 'last', 'min', 'max'};

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'extractMode', ...
    @(x) any(validatestring(x, validExtractModes)));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vecs, extractMode, varargin{:});
% param1 = iP.Results.param1;

% Validate extractMode
extractMode = validatestring(extractMode, validExtractModes);

%% Do the job
if iscell(vecs)
    elements = cellfun(@(x) extract_element(x, extractMode), vecs);
else
    elements = arrayfun(@(x) extract_element(vecs(:, x), extractMode), ...
                        transpose(1:size(vecs, 2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function element = extract_element (x, extractMode)

switch extractMode
    case 'first'
        element = x(1);
    case 'last'
        element = x(end);
    case 'min'
        element = min(x);
    case 'max'
        element = max(x);
    otherwise
        error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%