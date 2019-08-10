function varargout = all_dependent_functions (mFileName, varargin)
%% Prints all dependent files used by a given MATLAB script/function
% Usage: [fTableOrList, pTableOrList] = all_dependent_functions (mFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       all_dependent_functions('parse_pulse')
%       all_dependent_functions('parse_pulse.m')
%       [fTable, pTable] = all_dependent_functions('parse_pulse')
%       TODO:
%       [fList, pList] = all_dependent_functions('parse_pulse', 'OriginalOutput', true)
%
% Outputs:
%       fTableOrList    - a table or cell array of function paths
%                       specified as a table or cell
%       pTableOrList    - a table or structure of MATLAB products
%                       specified as a table or struct
% Arguments:
%       mFileName   - .m file name
%                   must be a string scalar or a character vector
%       varargin    - 'TopOnly': display only the functions used directly 
%                                   by the given script/function
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-09 Created by Adam Lu
% TODO: Add 'SaveFlag' as an optional argument
% TODO: Add 'PrintFlag' as an optional argument
% TODO: Add 'OutFolder' as an optional argument
% TODO: Add 'FunctionListPath' as an optional argument
% TODO: Add 'MatlabProductListPath' as an optional argument
% TODO: Add 'OriginalOutput' as an optional argument
% TODO: Add 'SheetType' as an optional argument

%% Hard-coded parameters
functionListSuffix = 'function_list';
matlabProductListSuffix = 'matlab_product_list';

% TODO: Make these optional arguments
saveFlag = true;
printFlag = true;
outFolder = pwd;
functionListPath = [];
matlabProductListPath = [];
originalOutput = false;
sheetType = 'csv';

%% Default values for optional arguments
topOnlyDefault = [];

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
addRequired(iP, 'mFileName');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TopOnly', topOnlyDefault);

% Read from the Input Parser
parse(iP, mFileName, varargin{:});
topOnly = iP.Results.TopOnly;

%% Preparation 
% Extract just the file name without '.m'
if regexp(mFileName, '.m$')
    mFileBase = extractBefore(mFileName, '.m');
else
    mFileBase = mFileName;
end

% Set default function list path
if isempty(functionListPath)
    functionListPath = [mFileBase, '_', functionListSuffix, '.', sheetType];
end

% Set default MATLAB product list path
if isempty(matlabProductListPath)
    matlabProductListPath = [mFileBase, '_', matlabProductListSuffix, '.', sheetType];
end

%% Do the job
% Retrieve dependent functions and MATLAB products
%   Note: Introduced in R2014a
if topOnly
    [functionList, matlabProductList] = ...
        matlab.codetools.requiredFilesAndProducts(mFileName, 'toponly');
else
    [functionList, matlabProductList] = ...
        matlab.codetools.requiredFilesAndProducts(mFileName);
end

if originalOutput
    varargout{1} = functionList;
    varargout{2} = matlabProductList;
    return
end

% Rename as fullPath and force as a column cell array
fullPath = force_column_cell(functionList);

% Extract just the function name
functionName = extract_fileparts(fullPath, 'base');

% Extract the containing directory
directory = extract_fileparts(fullPath, 'directory');

% Extract the common directory for all dependent functions
commonDirectory = extract_common_directory(functionList);

% Add the file separater
parentDirectoryWithFileSep = [commonDirectory, filesep];

% Extract relative paths
relativePath = extractAfter(fullPath, parentDirectoryWithFileSep);

%% Convert to tables
% Expand commonDirectory as cell array
commonDirectory = repmat({commonDirectory}, size(fullPath));

% Convert to tables
functionListTable = table(functionName, directory, commonDirectory, relativePath, fullPath);
matlabProductListTable = struct2table(matlabProductList);

%% Save tables
if saveFlag
    writetable(functionListTable, functionListPath);
    writetable(matlabProductListTable, matlabProductListPath);
end

%% Display results
if printFlag
    functionListTable
    matlabProductListTable
end

%% Output results
varargout{1} = functionListTable;
varargout{2} = matlabProductListTable;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%