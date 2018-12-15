function vectorsCell = force_column_cell (vectorsOrig)
%% Transforms a row cell array or a numeric array to a column cell array
% Usage: vectorsCell = force_column_cell (vectorsOrig)
% Explanation:
%       This is an attempt to standardize the way multiple vectors are stored
%           -- always as column cell arrays
%       1. Row cell arrays are converted to column cell arrays
%       2. Empty numeric arrays and character arrays are placed in a cell array
%       3. Numeric vector arrays are forced as a column vector
%           (force_column_numeric.m is used), then placed in a cell array
%       4. Numeric non-vector arrays are transformed to a column cell array
%           of column numeric vectors based on the original columns
%
% Example(s):
%       vectors = force_column_cell(data);
%       vectors = force_column_cell(vectors);
%
% Outputs:
%       vectorsCell - vectors as a column cell array
%                   specified as a column cell array
%
% Arguments:
%       vectorsOrig - original vectors
%                   Note: If an array, each column is a vector 
%                           to be placed in a cell
%                   must be a numeric array or a cell array 
%                       or a character vector
% Requires:
%       cd/count_vectors.m
%       cd/force_column_numeric.m
%
% Used by:
%       cd/compute_rms_error.m
%       cd/construct_fullpath.m
%       cd/filter_and_extract_pulse_response.m
%       cd/find_pulse_endpoints.m
%       cd/force_column_numeric.m
%       cd/force_row_numeric.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_individual_traces.m
%       cd/match_format_vectors.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m
%       cd/run_neuron.m

% File History:
% 2018-10-10 Created by Adam Lu
% 2018-10-19 Now accepts character vectors
% 2018-10-27 Now places empty numeric arrays in a cell array
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
addRequired(iP, 'vectorsOrig', ...
    @(x) isnumeric(x) || iscell(x) || ischar(x));

% Read from the Input Parser
parse(iP, vectorsOrig);

%% Do the job
if iscell(vectorsOrig) && ~iscolumn(vectorsOrig)
    % Reassign as a column
    vectorsCell = vectorsOrig(:);
elseif ischar(vectorsOrig)
    % Place in a cell array
    vectorsCell = {vectorsOrig};
elseif isnumeric(vectorsOrig)
    if isempty(vectorsOrig)
        % Place in a cell array
        vectorsCell = {vectorsOrig};
    else
        % Force any numeric vector as a column vector
        vectorsOrig = force_column_numeric(vectorsOrig, ...
                                            'IgnoreNonVectors', true);

        % Count the number of vectors
        nVectors = count_vectors(vectorsOrig);

        % Extract as a cell array
        vectorsCell = arrayfun(@(x) vectorsOrig(:, x), ...
                                transpose(1:nVectors), ...
                                'UniformOutput', false);
    end
else
    % Do nothing
    vectorsCell = vectorsOrig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

nVectors = size(vectorsOrig, 2);

% Reassign as a column
vectorsCell = vectorsCell(:);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%