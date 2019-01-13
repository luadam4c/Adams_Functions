function isNum = isnum (x, varargin)
%% Returns whether the input is numeric in the general sense (numeric, logical, datetime or duration)
% Usage: isNum = isnum (x, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       isNum       - whether the input is numeric in the general sense
%                   specified as a logical scalar
% Arguments:
%       x           - an input to check
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_axis_limits.m
%       cd/compute_combined_trace.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/isnumericvector.m
%       cd/compute_combined_data.m
%       cd/plot_horizontal_line.m
%       cd/plot_vertical_line.m

% File History:
% 2018-12-28 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isNum = isnumeric(x) || islogical(x) || isdatetime(x) || isduration(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
