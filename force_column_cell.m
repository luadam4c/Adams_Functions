function vectorsCell = force_column_cell (vectorsOrig)
%% Transforms a row cell array or a numeric array to a column cell array
% Usage: vectorsCell = force_column_cell (vectorsOrig)
% Explanation:
%       This is an attempt to standardize the way multiple vectors are stored
%           -- always as column cell arrays
% Example(s):
%       vectors = force_column_cell(data);
%       vectors = force_column_cell(vectors);
% Outputs:
%       vectorsCell - vectors as a column cell array
%                   specified as a column cell array
%
% Arguments:
%       vectorsOrig - original vectors
%                   Note: If an array, each column is a vector 
%                           to be placed in a cell
%                   must be a numeric array or a cell array
%
% Used by:
%       cd/count_samples.m
%       cd/match_vector_numbers.m
%       cd/parse_pulse.m
%       cd/parse_pulse_response.m

% File History:
% 2018-10-10 Created by Adam Lu
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
    @(x) isnumeric(x) || iscell(x));

% Read from the Input Parser
parse(iP, vectorsOrig);

%% Do the job
if iscell(vectorsOrig) && ~iscolumn(vectorsOrig)
    % Reassign as a column
    vectorsCell = vectorsOrig(:);
elseif isnumeric(vectorsOrig)
    % Count the number of vectors
    nVectors = size(vectorsOrig, 2);

    % Extract as a cell array
    vectorsCell = arrayfun(@(x) vectorsOrig(:, x), 1:nVectors, ...
                            'UniformOutput', false);

    % Reassign as a column
    vectorsCell = vectorsCell(:);
else
    % Do nothing
    vectorsCell = vectorsOrig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%