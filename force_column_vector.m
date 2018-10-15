function vectors = force_column_vector (vectors)
%% Transform row numeric vector(s) to column numeric vector(s)
% Usage: vectors = force_column_vector (vectors)
% Explanation:
%       Starting with a cell array of mixed vectors, some row and some column,
%           this function makes sure each vector is a column vector.
%       If a single vector is provided, the function makes sure 
%           it's a column vector.
% Example(s):
%       vector = force_column_vector(vector);
%       vectors = force_column_vector(vectors);
% Outputs:
%       vectors     - vectors transformed
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
%
% Arguments:
%       vectors     - original vectors
%                   must be a numeric vector or a cell array of numeric vectors
%
% Used by:    
%       cd/compute_average_trace.m
%       cd/fit_2exp.m
%       cd/plot_cfit_pulse_response.m

% File History:
% 2018-10-12 Created by Adam Lu
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
addRequired(iP, 'vectors', ...
    @(x) isnumeric(x) && isvector(x) || ...
        iscell(x) && all(cellfun(@(x) isnumeric(x) && isvector(x), x)) );

% Read from the Input Parser
parse(iP, vectors);

%% Do the job
if isnumeric(vectors) && ~iscolumn(vectors)
    % Reassign as a column
    vectors = vectors(:);
elseif iscell(vectors)
    % Extract as a cell array
    vectors = cellfun(@force_column_vector, vectors, 'UniformOutput', false);
else
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%