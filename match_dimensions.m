function arrayNew = match_dimensions (arrayOld, dimNew, varargin)
%% Reshapes or expands an array to match given dimensions
% Usage: arrayNew = match_dimensions (arrayOld, dimNew, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       arrayNew    - array matched
%                   specified as a numeric, logical, cell or struct array
% Arguments:    
%       arrayOld    - array to match
%                       if a character array, put in a cell
%                   must be a numeric, logical, char, cell or struct array
%       dimNew      - new dimensions
%                   must be a positive integer vector
%
% Requires:
%       cd/ispositiveintegervector.m
%
% Used by:    
%       cd/compute_weighted_average.m
%       cd/create_indices.m
%       cd/create_time_vectors.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/match_format_vector_sets.m
%       cd/match_array_counts.m
%       cd/normalize_by_initial_value.m
%       cd/parse_assyst_swd.m
%       cd/parse_atf_swd.m
%       cd/parse_pulse_response.m
%       cd/plot_traces_abf.m

% File History:
% 2018-10-24 Created by Adam Lu
% 2018-10-25 Changed the second argument to dimNew
% 

%% Hard-coded parameters

%% Default values for optional arguments

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
addRequired(iP, 'arrayOld', ...
    @(x) validateattributes(x, ...
        {'numeric', 'logical', 'char', 'cell', 'struct'}, {'3d'}));
addRequired(iP, 'dimNew', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'vector'}));

% Read from the Input Parser
parse(iP, arrayOld, dimNew, varargin{:});

%% Preparation
% Place character arrays in a cell array
if ischar(arrayOld)
    arrayOld = {arrayOld};
end

% Query the old dimensions
dimOld = size(arrayOld);

% If the old array is empty 
%   or if the new dimensions are the same as the old ones, 
%   just return the old array
if isempty(arrayOld) || isequal(dimNew, dimOld)
    arrayNew = arrayOld;
    return
end

% Query the number of dimensions
nDimsOld = length(dimOld);
nDimsNew = length(dimNew);

% Query the number of elements
nElementsOld = prod(dimOld);
nElementsNew = prod(dimNew);

%% Do the job
% Decide based on the relative number of elements
if nElementsOld == nElementsNew
    % Match dimensions by reshaping if there are equal number of elements
    if nDimsOld == nDimsNew
        % Reshape arrayOld to match dimNew
        arrayNew = reshape(arrayOld, dimNew);
    elseif nDimsOld > nDimsNew
        % Squeeze arrayOld, then reshape to match dimNew
        arrayNew = reshape(squeeze(arrayOld), dimNew);    
    elseif nDimsOld < nDimsNew
        % Look for the minimum in dimNew
        [~, indTemp] = min(dimNew);

        % Choose the first minimum dimension
        idxMinDim = indTemp(1);

        % Expand arrayOld
        if idxMinDim == 1
            arrayNew(1, :, :) = arrayOld;
        elseif idxMinDim == 2
            arrayNew(:, 1, :) = arrayOld;
        elseif idxMinDim == 3
            arrayNew(:, :, 1) = arrayOld;
        else
            error('Code logic error!');
        end
    else
        error('Code logic error!');
    end
elseif nElementsOld < nElementsNew
    % If there are fewer elements in the array than required, 
    %   try expanding it
    if nDimsOld == nDimsNew
        % Get the factor to expand for each dimension
        factorToExpand = dimNew ./ dimOld;

        if ispositiveintegervector(factorToExpand)
            % If all factors are positive integers, use repmat
            arrayNew = repmat(arrayOld, factorToExpand);
        else
            % Otherwise, return error
            error(['Factor to expand not all integers.\n', ...
                    'Array ''%s'' cannot be expanded to match ', ...
                    'the requested dimensions!'], inputname(1));
        end
    else
        error('Not implemented yet!');
    end
else
    error(['There are more elements in the array ', ...
            'than possible for the requested dimensions!']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function arrayOld = match_dimensions(arrayOld, array2, varargin)
%       array2      - array to be matched
%                   must be a cell array                   
addRequired(iP, 'array2', ...
    @(x) validateattributes(x, {'cell'}, {'3d'}));
parse(iP, arrayOld, array2, varargin{:});

dimNew = size(array2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
