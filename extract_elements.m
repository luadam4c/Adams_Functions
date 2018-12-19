function [elements, idxElement] = extract_elements (vecs, extractMode, varargin)
%% Extracts elements from vectors using a certain mode ('first', 'last', 'min', 'max')
% Usage: [elements, idxElement] = extract_elements (vecs, extractMode, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       elements    - element(s) from each vector extracted
%                   specified as a numeric vector 
%                       or a cell array of column numeric vectors
%       idxElement  - indices(s) of elements extracted
%                   specified as a numeric vector 
%                       or a cell array of column numeric vectors
% Arguments:
%       vecs        - vector(s)
%                   must be a numeric array or a cell array
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
%       cd/create_error_for_nargin.m
%       cd/create_indices.m
%       cd/parse_pulse_response.m
%       cd/plot_protocols.m

% File History:
% 2018-12-15 Created by Adam Lu
% 2018-12-17 Now returns idxElement as well
% TODO: Add 'MaxNum' as an optional argument with default Inf
% TODO: Add 'Indices', 'Endpoints' and 'Windows' as optional arguments
%           and use extract_subvectors.m
% 

%% Hard-coded parameters
validExtractModes = {'first', 'last', 'min', 'max'};

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs', ...                  % vectors to extract
    @(x) assert(isnumeric(x) || iscell(x), ...
                ['vecs must be either a numeric array', ...
                    'or a cell array of vectors!']));
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
switch extractMode
case {'first', 'last', 'min', 'max'}
    if iscell(vecs)
        % Do for all elements
        [elements, idxElement] = ...
            cellfun(@(x) extract_single_element(x, extractMode), vecs);
    else
        % Do for all columns
        [elements, idxElement] = ...
            arrayfun(@(x) extract_single_element(vecs(:, x), extractMode), ...
                    transpose(1:size(vecs, 2)));
    end
otherwise
    error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [element, idxElement] = extract_single_element (x, extractMode)

switch extractMode
    case 'first'
        element = x(1);
        idxElement = 1;
    case 'last'
        element = x(end);
        idxElement = length(x);
    case 'min'
        [element, idxElement] = min(x);
    case 'max'
        [element, idxElement] = max(x);
    otherwise
        error('Code logic error!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

@(x) assert(isnumeric(x) || iscellnumeric(x), ...
            ['vecs must be either a numeric array', ...
                'or a cell array of numeric arrays!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%