function arrayNew = match_row_count (arrayOld, nRowsNew, varargin)
%% Expands or truncates an array to match a given number of rows (dimension #1)
% Usage: arrayNew = match_row_count (arrayOld, nRowsNew, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       match_row_count([1, 2, 3], 6)
%       match_row_count([1; 2; 3], 6)
%       match_row_count([1, 2; 3, 4], 7)
%       match_row_count([1, 2; 3, 4], 7, 'ExpansionMethod', 'repeat')
%       match_row_count([1, 2; 3, 4], 7, 'ExpansionMethod', 'patchNaNs')
%
% Outputs:
%       arrayNew    - array matched
%                   specified as a numeric, cell or struct array
%
% Arguments:    
%       arrayOld    - array to match
%                   must be a numeric, cell or struct array
%       nRowsNew    - new number of rows
%                   must be a positive integer scalar
%       varargin    - 'ExpansionMethod': method for expanding vector
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'repeat'        - repeat the vector
%                       'patchNaNs'     - patch with NaNs
%                       'patchZeros'    - patch with zeros
%                       'patchOnes'    - patch with ones
%                   default == 'repeat'
%                   - 'DimToMatch': dimension to match
%                   must be empty or a positive integer scalar
%                   default == 1 (match rows)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_lts_errors.m
%       cd/compute_single_neuron_errors.m
%       cd/compute_sampsizepwr.m
%       cd/compute_sweep_errors.m
%       cd/count_samples.m
%       cd/decide_on_colormap.m
%       cd/match_format_vector_sets.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_create_sim_params.m
%       cd/m3ha_plot_figure07.m
%       cd/match_column_count.m
%       cd/match_format_vectors.m
%       cd/match_time_info.m
%       cd/parse_current_family.m
%       cd/parse_ipsc.m
%       cd/parse_repetitive_pulses.m
%       cd/parse_spike2_mat.m
%       cd/plot_fitted_traces.m
%       cd/plot_selected.m
%       cd/plot_struct.m
%       cd/xolotl_add_current_injection.m
%       cd/xolotl_add_voltage_clamp.m
%       cd/xolotl_set_simparams.m

% File History:
% 2018-10-25 Created by Adam Lu
% 2019-02-20 Now uses create_error_for_nargin
% 2019-08-15 Added 'ExpansionMethod' as an optional argument
% 2019-08-15 Implemented 'patchNaNs', 'patchZeros', 'patchOnes'
% 2020-01-02 Added 'DimToMatch' as an optional argument
% 2020-01-02 Now accepts any array type

%% Hard-coded parameters
validExpansionMethods = {'repeat', 'patchNaNs', 'patchZeros', 'patchOnes'};

%% Default values for optional arguments
expansionMethodDefault = 'repeat';
dimToMatchDefault = [];       % use the sum() function default by default

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
addRequired(iP, 'arrayOld');
addRequired(iP, 'nRowsNew', ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ExpansionMethod', expansionMethodDefault, ...
    @(x) any(validatestring(x, validExpansionMethods)));
addParameter(iP, 'DimToMatch', dimToMatchDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['DimToMatch must be either empty ', ...
                    'or a positive integer scalar!']));

% Read from the Input Parser
parse(iP, arrayOld, nRowsNew, varargin{:});
expansionMethod = validatestring(iP.Results.ExpansionMethod, validExpansionMethods);
dimToMatch = iP.Results.DimToMatch;

%% Preparation
% Set default dimension to match
if isempty(dimToMatch)
    dimToMatch = 1;
end

% Query the old number of rows
nRowsOld = size(arrayOld, dimToMatch);

% If arrayOld is empty 
%   or if the new number of rows are the same as the old ones, 
%   just return the old array
if isempty(arrayOld) || nRowsNew == nRowsOld
    arrayNew = arrayOld;
    return
end

% Query the number of dimensions
nDims = ndims(arrayOld);

%% Expand or truncate array
if nRowsNew > nRowsOld
    switch expansionMethod
    case 'repeat'
        % Compute the factor to expand by
        factorToExpand = floor(nRowsNew / nRowsOld);

        % Compute the remaining number of rows to fill after expansion
        remainingNRows = mod(nRowsNew, nRowsOld);

        % Expand array by repetition
        if nDims == 2
            switch dimToMatch
            case 1
                % First expand rows by factorToExpand
                arrayNew = repmat(arrayOld, [factorToExpand, 1]);

                % Fill in remaining rows by the first rows
                arrayNew = cat(1, arrayNew, arrayOld(1:remainingNRows, :));
            case 2
                % First expand columns by factorToExpand
                arrayNew = repmat(arrayOld, [1, factorToExpand]);

                % Fill in remaining columns by the first columns
                arrayNew = cat(2, arrayNew, arrayOld(:, 1:remainingNRows));
            case 3
                % First expand columns by factorToExpand
                arrayNew = repmat(arrayOld, [1, 1, factorToExpand]);
            otherwise
                error('dimToMatch unrecognized!');
            end
        elseif nDims == 3
            switch dimToMatch
            case 1
                % First expand by factorToExpand
                arrayNew = repmat(arrayOld, [factorToExpand, 1, 1]);

                % Fill in remaining rows by the first rows
                arrayNew = cat(1, arrayNew, arrayOld(1:remainingNRows, :, :));
            case 2
                % First expand columns by factorToExpand
                arrayNew = repmat(arrayOld, [1, factorToExpand, 1]);

                % Fill in remaining columns by the first columns
                arrayNew = cat(2, arrayNew, arrayOld(:, 1:remainingNRows, :));
            case 3
                % First expand by factorToExpand
                arrayNew = repmat(arrayOld, [1, 1, factorToExpand]);

                % Fill in remaining rows by the first rows
                arrayNew = cat(3, arrayNew, arrayOld(:, :, 1:remainingNRows));
            otherwise
                error('dimToMatch unrecognized!');
            end
        end
    case {'patchNaNs', 'patchZeros'}
        % Get the old dimensions
        dimOld = size(arrayOld);

        % Set new dimensions
        if dimToMatch > nDims
            dimNew = [dimOld, 1];
        else
            dimNew = dimOld;
            dimNew(dimToMatch) = nRowsNew;
        end

        % Initialize as NaNs or zeros
        switch expansionMethod
        case 'patchNaNs'
            arrayNew = nan(dimNew);
        case 'patchZeros'
            arrayNew = zeros(dimNew);
        case 'patchOnes'
            arrayNew = ones(dimNew);
        end

        % Expand array by patching with NaNs
        if nDims == 2
            switch dimToMatch
            case 1
                arrayNew(1:nRowsOld, :) = arrayOld;
            case 2
                arrayNew(:, 1:nRowsOld) = arrayOld;
            case 3
                arrayNew(:, :, 1) = arrayOld;
            otherwise
                error('dimToMatch unrecognized!');
            end
        elseif nDims == 3
            switch dimToMatch
            case 1
                arrayNew(1:nRowsOld, :, :) = arrayOld;
            case 2
                arrayNew(:, 1:nRowsOld, :) = arrayOld;
            case 3
                arrayNew(:, :, 1:nRowsOld) = arrayOld;
            otherwise
                error('dimToMatch unrecognized!');
            end
        end
    otherwise
        error('ExpansionMethod unrecognized!');
    end
elseif nRowsNew < nRowsOld
    % Truncate array
    if nDims == 2
        switch dimToMatch
        case 1
            arrayNew = arrayOld(1:nRowsNew, :);
        case 2
            arrayNew = arrayOld(:, 1:nRowsNew);
        case 3
            arrayNew = arrayOld(:, :, 1:nRowsNew);
        otherwise
            error('dimToMatch unrecognized!');
        end
    elseif nDims == 3
        switch dimToMatch
        case 1
            arrayNew = arrayOld(1:nRowsNew, :, :);
        case 2
            arrayNew = arrayOld(:, 1:nRowsNew, :);
        case 3
            arrayNew = arrayOld(:, :, 1:nRowsNew);
        otherwise
            error('dimToMatch unrecognized!');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Query the dimensions
dims = size(arrayOld);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
