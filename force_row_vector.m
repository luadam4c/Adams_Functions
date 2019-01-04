function vectors = force_row_vector (vectors, varargin)
%% Transform column vector(s) or array(s) to row vector(s)
% Usage: vectors = force_row_vector (vectors, varargin)
% Explanation:
%       Starting with a cell array of mixed vectors, some row and some column,
%           this function makes sure each vector is a row vector.
%       If a single vector is provided, the function makes sure 
%           it's a row vector.
% Example(s):
%       vector = force_row_vector(vector);
%       vectors = force_row_vector(vectors);
%       force_row_vector({[3, 4], [5; 6], magic(3)})
% Outputs:
%       vectors     - vectors transformed
%                   specified as a numeric array 
%                       or a cell array of numeric vectors
%
% Arguments:
%       vectors     - original vectors
%       varargin    - 'IgnoreNonVectors': whether to ignore non-vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/force_column_cell.m
%
% Used by:    

% File History:
% 2018-10-25 Modified from force_column_vector.m
% 2018-10-27 Added 'IgnoreNonVectors' as an optional argument
% 

%% Default values for optional arguments
ignoreNonVectorsDefault = false;    % don't ignore non-vectors by default

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
addRequired(iP, 'vectors');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNonVectors', ignoreNonVectorsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vectors, varargin{:});
ignoreNonVectors = iP.Results.IgnoreNonVectors;

%% Do the job
if ~iscell(vectors) && ~isrow(vectors)
    if isempty(vectors)
        % Do nothing
    elseif isvector(vectors)
        % Must be a column vector, so transpose it
        vectors = transpose(vectors);
    else
        % Must be a non-vectors
        if ~ignoreNonVectors
            % TODO: Make this more efficient by modifying force_column_cell.m
            % Reassign as a column cell array of column vectors
            vectors = force_column_cell(vectors);

            % Change the column vectors to row vectors
            vectors = force_row_vector(vectors);
        end
    end
elseif iscell(vectors)
    % Extract as a cell array
    %   Note: this will have a recursive effect
    vectors = cellfun(@force_row_vector, vectors, 'UniformOutput', false);
else
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

vectors = reshape(vectors, 1, numel(vectors));

%                   must be a numeric vector or a cell array of numeric arrays
@(x) assert(isnum(x) || iscellnumeric(x), ...
            ['vectors must be either a numeric array', ...
                'or a cell array of numeric arrays!']));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
