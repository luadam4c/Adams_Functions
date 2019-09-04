function frame = create_empty_frame (height, width)
%% Creates an empty MATLAB movie frame
% Usage: frame = create_empty_frame (height, width)
% Explanation:
%       TODO
%
% Example(s):
%       create_empty_frame(3, 2)
%
% Outputs:
%       frame       - empty MATLAB frame, with fields:
%                       cdata    - RGB intensity data
%                       colormap - color map used
%                   specified as a structure
%
% Arguments:
%       height      - height of the frame
%                   must be a positive scalar
%       width       - width of the frame
%                   must be a positive scalar
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/read_all_frames.m

% File History:
% 2019-09-04 Adapted from https://www.mathworks.com/help/matlab/import_export/read-video-files.html
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'height', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addRequired(iP, 'width', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, height, width);

%% Do the job
frame = struct('cdata', zeros(height, width, 3, 'uint8'), 'colormap', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%