function varargout = extract_columns (arrays, varargin)
%% Extracts columns from arrays
% Usage: varargout = extract_columns (arrays, colNumbers (opt), varargin)
% Explanation:
%       TODO
% Example(s):
%       [a, b] = extract_columns({magic(3); ones(4)}, [1:3])
%       [a, b] = extract_columns({magic(3); ones(3); zeros(4)}, ...
%                               {[2, 3], [1:3], [1, 3]})
%       a = extract_columns({magic(3); ones(4)}, [1:3], 'OutputMode', 'single')
%
% Outputs:
%       varargout   - extracted column #1s, column #2s, etc.
%                       or extracted columns for each array 
%                           if 'OutputMode' is 'single'
%                   specified as a numeric column vector 
%                       or a cell array of numeric column vectors
%                       or a cell array of cell arrays of numeric column vectors
% Arguments:    
%       arrays      - arrays to extract columns from
%                   must be a numeric array
%                       or a cell array of numeric arrays
%       colNumbers  - (opt) column number(s) to extract instead of 1, 2, 3, ...
%                   must be either 'all' or a positive integer vector
%                       or a cell array of positive integer vectors
%                   default == 'all'
%       varargin    - 'OutputMode': mode of outputs
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'multiple' - all columns as separate outputs
%                       'single'   - all columns as one output
%                   default == 'multiple'
%
% Requires:
%       cd/force_column_numeric.m
%       cd/match_dimensions.m
%       cd/ispositiveintegervector.m
%       cd/iscellnumeric.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-24 Created by Adam Lu
% 2018-10-25 Updated usage of match_dimnensions()
% 

%% Hard-coded parameters
validOutputModes = {'multiple', 'single'};

%% Default values for optional arguments
colNumberDefault  = 'all';      % extract all columns by default
outputModeDefault = 'multiple'; % separate outputs by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'arrays', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['arrays must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'colNumbers', colNumberDefault, ...
    @(x) assert((ischar(x) || isstring(x)) && strcmpi(x, 'all') || ...
                ispositiveintegervector(x) || ...
                iscell(x) && iscellnumeric(x), ...
                ['colNumbers must be either ''all'', "all" ', ...
                    'or a positive integer vector ', ...
                    'or a cell array of positive integer vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputMode', outputModeDefault, ...
    @(x) any(validatestring(x, validOutputModes)));

% Read from the Input Parser
parse(iP, arrays, varargin{:});
colNumbers = iP.Results.colNumbers;
outputMode = validatestring(iP.Results.OutputMode, validOutputModes);

% Check if colNumbers is compatible with arrays
if iscell(colNumbers)
    if ~iscell(arrays)
        error(['colNumbers cannot be a cell array ', ...
                'if arrays is not a cell array!']);
    elseif numel(arrays) ~= numel(colNumbers)
        error(['If colNumbers is a cell array, ', ...
                'it must have the same number of elements as arrays!'])
    end
end

%% Preparation
% If arrays is empty, return empty outputs
if isempty(arrays)
    varargout = cell(1, nargout);
    return
end

% Count the number of arrays
if isnumeric(arrays)
    nArrays = 1;
else
    nArrays = numel(arrays);
end

% Count the number of columns
if isnumeric(arrays)
    nCols = size(arrays, 2);
else
    nCols = cellfun(@(x) size(x, 2), arrays);
end

% Make sure provided column numbers are within range
if iscell(colNumbers)
    % Get the maximum column numbers
    maxColNumbers = cellfun(@max, colNumbers);

    % If any maximum column number is greater than 
    %   the corresponding number of columns, return error
    if any(maxColNumbers > nCols)
        indProblematic = find(maxColNumbers > nCols);
        for idx = 1:length(indProblematic)
            error(['The maximum column number requested exceeded ', ...
                    'the minimum number of columns for array number %d!'], ...
                    indProblematic(idx));
        end
    end    
elseif isnumeric(colNumbers)
    % Get the maximum column number
    maxColNumber = max(colNumbers);

    % Get the minimum number of columns
    minNCols = min(nCols);

    % If the maximum column number is greater than 
    %   the minimum number of columns, return error
    if maxColNumber > minNCols
        error(['The maximum column number cannot exceed ', ...
                'the minimum number of columns!']);
    end
end

% Modify colNumbers to be compatible with arrays
if iscell(arrays)
    if iscell(colNumbers)
        % Match the dimensions of the cell array
        colNumbers = match_dimensions(colNumbers, size(arrays));
    elseif isnumeric(colNumbers)
        % Place colNumbers in a cell and repmat it as many times as 
        %   the number of arrays
        colNumbers = repmat({colNumbers}, size(arrays));
    end
elseif ischar(colNumbers) || isstring(colNumbers)
    % Generate column numbers
    if isnumeric(arrays)
        colNumbers = transpose(1:nCols);
    else
        colNumbers = cellfun(@(x) transpose(1:x), nCols, ...
                            'UniformOutput', false);
    end
end

% Count the number of columns requested
if isnumeric(colNumbers)
    nColNumbers = length(colNumbers);
elseif iscell(colNumbers)
    nColNumbers = cellfun(@length, colNumbers);
end

% Count the number of output arguments requested
switch outputMode
    case 'multiple'
        % Return as many outputs as needed
        nOutputs = nargout;

        % If the number of outputs requested is greater than 
        %   the number of columns requested, return error
        if nOutputs > max(nColNumbers)
            error(['The number of outputs requested ', ...
                    'cannot exceed the number of columns requested!']);
        end
    case 'single'
        % Return a single output
        nOutputs = 1;

        % If more than one output requested, return error
        if nargout > 1
            error('There can only be one output under ''single'' mode!');
        end
    otherwise
        error('outputMode unrecognized!');
end

%% Extract columns
switch outputMode
    case 'multiple'
        % Extract columns
        varargout = cell(1, nOutputs);
        for iOutput = 1:nOutputs
            if isnumeric(arrays)
                varargout{iOutput} = arrays(:, colNumbers(iOutput));
            else
                varargout{iOutput} = ...
                    cellfun(@(x, y) x(:, y(iOutput)), arrays, colNumbers, ...
                            'UniformOutput', false);
            end
        end
    case 'single'
        % Transform arrays into cell arrays of column vectors
        if isnumeric(arrays)
            varargout{1} = force_column_numeric(arrays(:, colNumbers));
        else
            varargout{1} = ...
                cellfun(@(x, y) force_column_numeric(x(:, y)), ...
                            arrays, colNumbers, ...
                            'UniformOutput', false);
        end
    otherwise
        error('outputMode unrecognized!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%