function [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
%% Saves each variable in a structure as a variable in a MAT-file and create a logHeader and a logVariables
% Usage: [m, fileName, logHeader, logVariables] = struct2mat (structure, varargin)
% Outputs: TODO
% Arguments: TODO
%       structure   - must be a structure
%       varargin    - 'FileName' - file name to save as, can include .mat or not    
%                   must be a string scalar or a character vector
%                   default == 'mystruct.mat'
%                   - 'LogHeader' - a cell array of labels 
%                                   for the fields of the structure 
%                   must be a cell array with the same length as 
%                       the number of variables in the structure TODO
%                   default == ''
%
% Used by:
%
% File History:
% 2016-11-07 Envisioned
% 2017-05-21 Created
% 2018-07-30 Implemented input parser 
% 2018-07-30 Add a parameter FileName to allow user to provide 
%               a custom file name
% 2018-07-30 Add a parameter LogHeader to allow user to provide 
%               a custom log header
% 

%% Hard-coded parameters
DEFAULT_FILENAME = 'my_struct.mat';

%% Default values for optional arguments
fileNameDefault = '';           % default file name
logHeaderDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments 
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'fit_logIEI';

% Add required inputs to the Input Parser
addRequired(iP, 'structure', ...              % a vector of IEIs
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileName', fileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'LogHeader', logHeaderDefault, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, structure, varargin{:});
fileName = iP.Results.FileName;
logHeader = iP.Results.LogHeader;

%% Perform task
% If no file name provided, create a file name
if isempty(fileName)
    % Find original structure name, if any
    structName = inputname(1);

    % Create file name
    if ~isempty(structName)             % original structure name provided
        % Use original structure name appended with '_struct' for file name
        fileName = [structName, '_struct.mat'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%