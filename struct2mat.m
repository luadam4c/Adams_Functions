function [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
%% Saves each variable in a structure as a variable in a MAT-file and create a logHeader and a logVariables
% Usage: [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
% Outputs: TODO
% Arguments: TODO
%       structure   - must be a structure
%       varargin    - 'FileName' - file name to save as, can include .mat or not    
%                   must be a TODO
%                   default == 'mystruct.mat'
%                   - 'LogHeader' - a cell array of labels 
%                                   for the fields of the structure 
%                   must be a cell array with the same length as 
%                       the number of variables in the structure
%                   default == ''
%
% Used by:
%
% File History:
% 2016-11-07 Envisioned
% 2017-05-21 Created
% TODO: Add a parameter FileName to allow user to provide a custom file name
%       fileName = iP.FileName;
% TODO: Add a parameter LogHeader to allow user to provide a custom log header
%       logHeader = iP.LogHeader;
% 

%% Constants
DEFAULT_FILENAME = 'mystruct.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments 
%% TODO: Use input Parser instead
if nargin < 1
    error('Not enough input arguments, type ''help struct2mat'' for usage');
elseif ~isstruct(structure)
    error('First argument must be a structure array!');
end

%% Perform task
% If no file name provided, create a file name
if isempty(fileName)
    % Find original structure name, if any
    structName = inputname(1);

    % Create file name
    if ~isempty(structName)             % original structure name provided
        % Use original structure name to for file name
        fileName = [structName, '.mat'];
    else                                % no original structure name provided
        % Use default file name
        fileName = DEFAULT_FILENAME;
    end
end

% Replace any original extension in fileName with 'mat'
[filePath, fileBase, fileExt] = fileparts(fileName);
fileName = fullfile(filePath, [fileBase, '.mat']);

% Save fields as individual variables in a MAT-file named by fileName
save(fileName, '-struct', 'structure', '-v7.3');    

% Return MAT-file object
m = matfile(fileName, 'Writable', true);    % MAT-file object for fileName

% Create a log of variables (the field names of the structure)
logVariables = fieldnames(structure);

% Create a header for the variable log
if isempty(logHeader)
    % Use field names as header
    logHeader = logVariables;
else
    % TODO: check length of logHeader: must be the same as number of variables
    %       otherwise, return error message
end

% Store in MAT-file
m.logHeader = logHeader;
m.logVariables = logVariables;

