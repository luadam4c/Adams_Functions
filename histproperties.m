function [harea, edges] = histproperties (X)
%% Computes the area, edges of the histogram for given data array
% Usage: [harea, edges] = histproperties (X)
% Example(s):
%       TODO
% Outputs:
%       harea       - histogram area
%                   specified as a numeric scalar
%       edges       - edges of the histogram bins
%                   specified as a numeric vector
% Arguments:
%       X           - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%
% Used by:    
%       cd/plot_pdf.m
%
% File History:
% 2018-06-12 Created by Adam Lu

%% Hard-coded parameters

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
addRequired(iP, 'X', ...                    % data to distribute among bins
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, X);

%% Do the job
% Compute the histogram edges
[~, edges] = histcounts(X);

% Compute the binWidth
binWidth = edges(2) - edges(1);

% Compute the total number of data points
npts = length(X);

% Compute the total histogram area
harea = binWidth * npts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       npts        - number of data points
%                   specified as a positive integer scalar

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
