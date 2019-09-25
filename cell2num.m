function numArray = cell2num (cellArray, varargin)
%% This is the reverse of num2cell, replacing empty entries with NaNs
% Usage: numArray = cell2num (cellArray, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       cell2num({[], 2, [3, 4, 5], []})
%       cell2num({[], 2; [3, 4, 5], []})
%       cell2num({[], 2; [3, 4, 5], []}, 'CombineMethod', 'nan')
%       cell2num({[], 2; [3, 4, 5], []}, 'CombineMethod', 'mean')
%
% Outputs:
%       numArray    - numeric array
%                   specified as a numeric array
% Arguments:
%       cellArray   - cell array with at most one number per cell
%                   must be a cell array of numeric vectors
%       varargin    - 'CombineMethod': method for combining numbers 
%                                       if there are more in the same cell
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - take the first value
%                       'nan'   - return nan if more than one number
%                       'average' or 'mean' - take the average value
%                   default == 'first'
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/isnum.m
%
% Used by:
%       cd/compute_combined_trace.m
%       cd/compute_sampsizepwr.m
%       cd/plot_tuning_curve.m

% File History:
% 2019-08-20 Created by Adam Lu
% 2019-09-25 Added 'CombineMethod' as an optional argument
% 

%% Hard-coded parameters
validCombineMethods = {'first', 'nan', 'average', 'mean'};

%% Default values for optional arguments
combineMethodDefault = 'first';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'cellArray');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'CombineMethod', combineMethodDefault, ...
    @(x) any(validatestring(x, validCombineMethods)));

% Read from the Input Parser
parse(iP, cellArray, varargin{:});
combineMethod = validatestring(iP.Results.CombineMethod, validCombineMethods);

%% Do the job
numArray = cellfun(@(x) force_numeric_scalar(x, combineMethod), cellArray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = force_numeric_scalar(value, combineMethod)

if isempty(value) || ~isnum(value)
    out = NaN;
elseif numel(value) > 1
    switch combineMethod
        case 'first'
            out = value(1);
        case 'nan'
            out = NaN;
        case {'mean', 'average'}
            out = nanmean(value);
        otherwise
            error('combineMethod unrecognized!');
    end
else
    out = value;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
