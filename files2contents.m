function newCellStr = files2contents (oldCellStr, varargin)
%% Replaces file names with file contents in a cell array of strings
% Usage: newCellStr = files2contents (oldCellStr, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       newCellStr  - TODO: Description of output1
%                   specified as a TODO
% Arguments:    
%       oldCellStr  - TODO: Description of oldCellStr
%                   must be a TODO
%
% Used by:    
%       cd/run_neuron.m

% File History:
% 2018-10-21 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments

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
addRequired(iP, 'oldCellStr', ...
    @iscellstr);

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, oldCellStr, varargin{:});

%% Do the job
% Count the number of strings
nStrings = numel(oldCellStr);

% Convert all file names to contents
newCellStr = cell(size(oldCellStr));
%parfor iString = 1:nStrings
for iString = 1:nStrings
    % Get this string
    thisStr = oldCellStr{iString};

    % If this string is a file name, read and place in newCellStr,
    %   otherwise, just copy over this string
    if isfile(thisStr)
        newCellStr{iString} = fileread(thisStr);
    else
        newCellStr{iString} = thisStr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%