function create_input_file (dataDir, defParams, varargin)
%% Create an input Excel file from data file names in a directory based on default parameters
% Usage: create_input_file (dataDir, defParams, varargin)
% Explanation:
%       TODO
% Side Effects:
%       Generates an Excel file in the present working directory
% Arguments:    
%       dataDir     - a data directory
%                   must be a directory
%       defParams   - default parameters
%                   must be a structure with fields containing scalar values
%       varargin    - 'DataType': data type
%                   must be a scalartext
%                   default == 'mat'
%
% Requires:
%
% Used by:    
%
% File History:
% 2018-04-30 Created by Adam Lu
% 

%% Hard-coded parameters
dataFileHeader = 'File Names';      % header for the data file names

% Default parameters
% defParams.WL = 0.5;         % Window size for standard deviation computations (in seconds)
% defParams.SNR = sqrt(2);    % Signal to noise ratio
% defParams.minD = 2;         % Minimum duration (in seconds)
% defParams.minF = 6.0;       % Minimum frequency (in Hz)
% defParams.maxF = 11;        % Maximum frequency (in Hz)
% defParams.FP = 2;           % Frequency proximity (in Hz) 
%                             %   i.e., how close the harmonics are to each other

%% Default values for optional arguments
dataTypeDefault = 'mat';            % detect matfiles by default

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
addRequired(iP, 'dataDir', ...                  % a data directory
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'defParams', ...                % default parameters
    @(x) validateattributes(x, {'struct'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'DataType', dataTypeDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, dataDir, defParams, varargin{:});
dataType = iP.Results.DataType;

%% Deal with files
% Check if the data directory exists and exit the program
if ~isdir(dataDir)
    fprintf('The data directory %s does not exist!!\n', dataDir);
    return;
end

% Look for files of the given data type in the data directory
%   This returns a row structure array
dataFiles = dir(['*.', dataType]);

% Extract the file names as a row cell array, sorted alphabetically
dataFileNames = sort({dataFiles.name});

% Display a warning if no files of the given data type are found
if isempty(dataFiles)
    fprintf('Warning: There are no .%s files in %s!!\n', dataType, dataDir);
end

% Count the number of files
nFiles = numel(dataFiles);

%% Deal with parameters
% Extract the parameter names
%   This returns a column cell array
paramNames = fieldnames(defParams);

% Extract the parameter values
% paramValues = struct2cell(defParams);
paramValues = cellfun(@num2str, struct2cell(defParams), ...
                        'UniformOutput', false);

% Count the number of parameters
nParams = numel(paramNames);

% Check if all fields in the default parameters structure are scalar values
if nParams > 0 && ~all(structfun(@isscalar, defParams))
    % Determine which fields are not scalar
    idxfields = find(~structfun(@isscalar, defParams));

    % Print all field names that are not scalars and exit the program
    for idx = idxfields
        fprintf('The parameter %s is not a scalar!!\n', paramNames{idx});
    end
    return;    
end

%% Create cell array
% Allocate a cell array with files as rows and parameters as columns
cellArray = cell(nFiles + 1, nParams + 1);

% Write the header for the data files
cellArray{1, 1} = dataFileHeader;

% Place the parameter names in the cell array as the header
cellArray(1, 2:(nParams + 1)) = paramNames';

% Place the file names in the first column of the cell array
cellArray(2:(nFiles + 1), 1) = dataFileNames';

% Place the values in each row of the cell array
cellArray(2:(nFiles + 1), 2:(nParams + 1)) = repmat(paramValues', nFiles, 1);

%% Write cell array
% Split the data directory path into parts
dataDirParts = strsplit(dataDir, '/');

% The last non-empty part is the data directory name
dataDirName = dataDirParts{find(cellfun(@numel, dataDirParts), 1, 'last')};

% Create input Excel file based on data directory name and data type
% inputFileName = [dataDirName, '_all', dataType, 'files.xls'];

% Write cell array to an Excel file in the present working directory
% xlswrite(fullfile(pwd, inputFileName), cellArray);

% Write cell array to a csv file in the present working directory
inputFileName = [dataDirName, '_all', dataType, 'files.csv'];
fid = fopen(inputFileName, 'w');
for iRow = 1:size(cellArray, 1)
    % Print the row
    fprintf(fid, [strjoin(cellArray(iRow, :), ', '), '\n']);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
TODO: Place older versions of the code that you want to save here, 
        in case you need it back in the future

%}
