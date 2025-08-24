function [y, ia, ic] = unique_custom (x, varargin)
%% Returns the unique values in x, optionally without NaN
% Usage: [y, ia, ic] = unique_custom (x, optArg (opt), varargin)
% Explanation:
%       This function is the same as the default unique() function
%           but treats all NaN as equal by default (so there would
%           be only one NaN if any in the output)
%
% Example(s):
%       unique_custom([5, NaN, 3, 5, NaN])
%       unique_custom([5, NaN, 3, 5, NaN], 'IgnoreNaN', true)
%       unique_custom([5, NaN, 3, 5, NaN], 'stable', 'IgnoreNaN', true)
%       unique_custom([5, NaN, 3, 5, NaN], 'TreatNanAsEqual', false)
%       unique_custom({'a', 'a', 'b', ''})
%       unique_custom({'a', 'a', 'b', ''}, 'IgnoreEmpty', true)
%       
% Outputs:
%       y           - All unique values in x
%                   specified as a array
% Arguments:
%       x           - Matrix to check unique values
%                   must be a array
%       optArg      - (opt) either setOrder or occurrence 
%                           for the unique() function
%                   must be an unambiguous, case-insensitive match to one of: 
%                       - ''
%                       - 'sorted'
%                       - 'stable'
%                       - 'first'
%                       - 'last'
%                   default == ''
%       varargin    - 'IgnoreNaN': whether to include NaN as distinct elements
%                                   for numeric arrays
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'IgnoreEmpty': whether to include empty strings 
%                                   as distinct elements for cell arrays
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatNanAsEqual': whether to treat all NaN values
%                                           as the same
%                       Note: If 'IgnoreNaN' == true, 
%                           'TreatNanAsEqual' has no effect
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for the unique() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/remove_empty.m
%
% Used by:
%       cd/adjust_edges.m
%       cd/compute_combined_trace.m
%       cd/compute_confidence_ellipse.m
%       cd/combine_swd_resp_data.m
%       cd/convert_to_rank.m
%       cd/create_default_grouping.m
%       cd/extract_subvectors.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_phase_info.m
%       cd/plot_measures.m
%       cd/plot_tuning_curve.m
%       cd/unique_groups.m

% File History:
% 2019-04-01 BT - Adapted from https://www.mathworks.com/matlabcentral/
%                         answers/42561-treating-nan-as-a-
%                         unique-value-instead-of-as-a-distinct#answer_52371
% 2019-07-24 Added optArg
% 2019-07-24 Improved algorithm
% 2019-10-09 Added 'IgnoreEmpty' as an optional argument

%% Hard-coded parameters
validOptArgs = {'', 'sorted', 'stable', 'first', 'last'};

%% Default values for optional arguments
optArgDefault = '';
ignoreNanDefault = false;  	    % do not ignore NaN by default
ignoreEmptyDefault = false;     % do not ignore empty by default
treatNanAsEqualDefault = true;  % treat all NaN values equal by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'x');

% Add optional inputs to the Input Parser
addOptional(iP, 'optArg', optArgDefault, ...
    @(x) any(validatestring(x, validOptArgs)));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreNaN', ignoreNanDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'IgnoreEmpty', ignoreEmptyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'TreatNanAsEqual', treatNanAsEqualDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));

% Read from the Input Parser
parse(iP, x, varargin{:});
optArg = validatestring(iP.Results.optArg, validOptArgs);
ignoreNan = iP.Results.IgnoreNaN;
ignoreEmpty = iP.Results.IgnoreEmpty;
treatNanAsEqual = iP.Results.TreatNanAsEqual;

% Keep unmatched arguments for the unique_custom() function
otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments

%% Preparation
% Deal with NaNs of empty entries
if isnum(x)
    if ignoreNan
        % Remove NaNs from the data
        x = remove_empty(x);
    elseif treatNanAsEqual
        % Decide if there is an NaN
        if any(isnan(x))
            % Remove NaNs from the data
            x = remove_empty(x);

            % Place one at the end
            if iscolumn(x)
                x = [x; NaN];
            else
                x = [x, NaN];
            end
        end
    end
elseif iscell(x)
    if ignoreEmpty
        x = remove_empty(x);
    end
end

%% Do the job
% Initial unique matrix y of x
if isempty(optArg)
    [y, ia, ic] = unique(x, otherArguments{:});
else
    [y, ia, ic] = unique(x, optArg, otherArguments{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% All NaN indices
indNaN = isnan(y(1:end));

% Delete all NaN elements in y, ia, and ic
y(indNaN) = [];
ia(indNaN) = [];
ic(indNaN) = [];

% NaN indices, does not include last NaN if present
indNaN = isnan(y(1:end-1));

% Delete all NaN elements in y, ia, and ic
y(indNaN) = [];
ia(indNaN) = [];
ic(indNaN) = [];
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
