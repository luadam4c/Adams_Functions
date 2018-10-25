function vectors = force_row_numeric (vectors)
%% Transform column numeric vector(s) or numeric array(s) to row numeric vector(s)
% Usage: vectors = force_row_numeric (vectors)
% Explanation:
%       Starting with a cell array of mixed vectors, some row and some column,
%           this function makes sure each vector is a row vector.
%       If a single vector is provided, the function makes sure 
%           it's a row vector.
% Example(s):
%       vector = force_row_numeric(vector);
%       vectors = force_row_numeric(vectors);
% Outputs:
%       vectors     - vectors transformed
%                   specified as a numeric array 
%                       or a cell array of numeric vectors
%
% Arguments:
%       vectors     - original vectors
%                   must be a numeric vector or a cell array of numeric arrays
%
% Requires:
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%
% Used by:    
%       cd/m3ha_compute_single_neuron_errors.m TODO

% File History:
% 2018-10-25 Modified from force_column_numeric.m
% 

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
addRequired(iP, 'vectors', ...                   % vectors
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vectors must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Read from the Input Parser
parse(iP, vectors);

%% Do the job
if isnumeric(vectors) && ~iscolumn(vectors)
    if isvector(vectors)
        % Reassign as a row
        vectors = reshape(vectors, 1, numel(vectors));
    else
        % TODO: Make this more efficient by modifying force_column_cell.m
        % Reassign as a column cell array of column vectors
        vectors = force_column_cell(vectors);

        % Change the column vectors to row vectors
        vectors = force_row_numeric(vectors);
    end
elseif iscell(vectors)
    % Extract as a cell array
    %   Note: this will have a recursive effect
    vectors = cellfun(@force_row_numeric, vectors, 'UniformOutput', false);
else
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%