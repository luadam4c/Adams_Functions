function [structArray, structNames] = structofstruct2structarray (structOfStruct, varargin)
%% Converts a structure of structures to a structure array
% Usage: [structArray, structNames] = structofstruct2structarray (structOfStruct, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       [a, b] = structofstruct2structarray(myStructOfStruct)
%
% Outputs:
%       structArray - a structure array
%                   specified as a structure array
%       structNames - a cell array of array names
%
% Arguments:
%       structOfStruct  - a structure of structures
%                       must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:

% File History:
% 2019-09-02 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'structOfStruct', @isstruct);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, structOfStruct, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Get all the structure names
structNames = fieldnames(structOfStruct);

% Extract all the structures
structArray = cellfun(@(x) structOfStruct.(x), structNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
