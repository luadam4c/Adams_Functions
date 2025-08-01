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
%       cell2num({[], 2, [3, 4, 5], []}, 'CombineMethod', 'padNaN')
%       cell2num({[], 2; [3, 4, 5], []}, 'CombineMethod', 'padNaN')
%
% Outputs:
%       numArray    - numeric array
%                   specified as a numeric array
% Arguments:
%       cellArray   - cell array of vectors (no more than 1 dimension)
%                   must be a cell array of numeric vectors
%       varargin    - 'CombineMethod': method for combining numbers 
%                                       if there are more in the same cell
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'first' - take the first value
%                       'nan'   - return nan if more than one number
%                       'average' or 'mean' - take the average value
%                       'padNaN'   - preserve all values and pad with NaN
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
%       vIRt\virt_moore.m

% File History:
% 2019-08-20 Created by Adam Lu
% 2019-09-25 Added 'CombineMethod' as an optional argument
% 2025-08-01 Added 'padNaN' as an optional argument
% 

%% Hard-coded parameters
validCombineMethods = {'first', 'nan', 'average', 'mean', 'padNaN'};

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

% Preparation
if strcmp(combineMethod, 'padNaN')
    % Count each array
    counts = cellfun(@numel, cellArray);

    % Obtain maximum length of array
    maxLength = max(max(counts));

    % Pad each cell content
    paddedArray = cellfun(@(x) pad_vector(x, maxLength), ...
                        cellArray, 'UniformOutput', false);

    % Combine vectors
    numArray = horzcat(paddedArray{:});
else
    % Force each cell content into a numeric scalar
    numArray = cellfun(@(x) force_numeric_scalar(x, combineMethod), cellArray);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = pad_vector(value, maxLength)

% Initialize out vector
out = nan(maxLength, 1);

% Replace first part of vector with value
out(1:numel(value)) = value;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = force_numeric_scalar(value, combineMethod)

if isempty(value) || ~isnum(value)
    out = NaN;
else
    if numel(value) > 1
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
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
