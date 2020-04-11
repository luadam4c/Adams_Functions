function scaleStr = islog2scale (isLog)
%% Converts whether is log-scaled (a Boolean) to a string ('log' or 'linear')
% Usage: scaleStr = islog2scale (isLog)
% Explanation:
%       TODO
%
% Example(s):
%       scaleStr = islog2scale(true)
%       scaleStrs = islog2scale([true; false])
%       scaleStrs = islog2scale([true, false])
%
% Outputs:
%       scaleStr    - a string ('log' or 'linear')
%                   specified as a character array 
%                       or a cell array of character arrays
%
% Arguments:
%       isLog       - whether is log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/plot_tuning_curve.m

% File History:
% 2020-04-11 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
if numel(isLog) < 1
    scaleStr = '';
elseif numel(isLog) == 1
    scaleStr = islog2scale_helper(isLog);
else
    scaleStr = arrayfun(@islog2scale_helper, isLog, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaleStr = islog2scale_helper (isLog)

if isLog
    scaleStr = 'log';
else
    scaleStr = 'linear';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%