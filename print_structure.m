function print_structure (structure, varargin)
%% Display all fields of a structure recursively
% Usage: print_structure (structure, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Side Effects:
%       Prints to standard output
% Arguments:    
%       structure   - the structure to print
%                   must be a structure
%       varargin    - 'FileID': file ID returned by fopen()
%                   must be an integer
%                   default == 1 (standard output)
%                   - 'NTabs': number of tabs before fields
%                   must be a positive integer scalar
%                   default == 1
%                   - 'StructName': name of structure 
%                                   if different from first argument
%                   must be a string scalar or a character vector
%                   default == inputname(1)
%
% Requires:
%       cd/print_cellstr.m
%
% Used by:
%       cd/find_passive_params.m
%
% 2016-11-02 Created
% 2018-06-21 Added comments and input parser
% 2018-06-21 Added cell string and now uses print_cellstr.m
% 2018-06-21 Now uses mat2str
% 2018-06-21 Added FileID, NTabs, StructName
% TODO: Deal with non-scalar structures
% TODO: Change disp() to fprintf for it to work when fileID is not 1

%% Default values for optional arguments
fileIdDefault = 1;                      % default: print to standard output
nTabsDefault = 1;                       % default number of tabs before fields
structNameDefault = inputname(1);       % default structure name

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
addRequired(iP, 'structure', @isstruct);        % the structure to print

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileId', fileIdDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'integer'}));
addParameter(iP, 'NTabs', nTabsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'StructName', structNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, structure, varargin{:});
fileId = iP.Results.FileId;
nTabs = iP.Results.NTabs;
structName = iP.Results.StructName;

%% Preparation
% Display warning if a nonscalar structure array is passed
if numel(structure) > 1
    fprintf(['Warning: Only the first entry of a non-scalar ', ...
            'structure array will be printed.\n\n']);
end

%% Perform task
% Print structure name
fprintf(fileId, [repmat('\t', 1, nTabs - 1), '''%s'':\n'], structName);

% Prepare tabs
tabs = repmat('\t', 1, nTabs);

% Get all field names
fields = fieldnames(structure);

% Count the number of fields
nFields = numel(fields);

% Print all fields
for iField = 1:nFields
    % Get the current field name
    fieldName = fields{iField};

    % Get the current field value
    fieldValue = structure.(fieldName);

    % Print according to data type
    if isnumeric(fieldValue) || islogical(fieldValue)
        fprintf(fileId, [tabs, '''%s'': %s\n'], ...
                        fieldName, mat2str(fieldValue));
    elseif ischar(fieldValue)
        fprintf(fileId, [tabs, '''%s'': ''%s''\n'], ...
                        fieldName, fieldValue);
    elseif iscellstr(fieldValue)
        fprintf(fileId, [tabs, '''%s'': '], fieldName);
        print_cellstr(fieldValue, 'FileID', fileId);
    elseif isstruct(fieldValue)
        print_structure(fieldValue, 'FileID', fileId, ...
                                    'StructName', fieldName, ...
                                    'NTabs', nTabs + 1);
    else
        fprintf(fileId, [tabs, '''%s'':\n'], fieldName);
        disp(fieldValue);   % TODO: Will not work if fileId is not 1
    end

    % Place a spacing line between each field
    fprintf(fileId, '\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

elseif ~isstruct(structure)
    error('First argument must be a structure array!');
end

% This will print the variable name
display(fieldValue);

fprintf('''%s'' == %g\n', fieldName, fieldValue);

elseif islogical(fieldValue) && isscalar(fieldValue)
    if fieldValue
        logicString = 'true';
    else
        logicString = 'false';
    end
    fprintf('''%s'' == %s\n', fieldName, logicString);        

fprintf('''%s'' is a cell array with the contents:\n', fieldName);
print_cellstr(fieldValue);

fprintf(fileId, 'Printing structure ''%s'':\n\n', inputname(1));
fprintf(fileId, '''%s'' == %s\n', fieldName, mat2str(fieldValue));
fprintf(fileId, '''%s'' == ''%s''\n', fieldName, fieldValue);
fprintf(fileId, '''%s'' is a cell array with the contents:\n', fieldName);
fprintf(fileId, '''%s'' has value:\n', fieldName);

fprintf(fileId, [repmat('\t', 1, nTabs - 1), '''%s'':\n'], inputname(1));

%% Add directories to search path for required functions across servers
if ~isdeployed
    if exist(fullfile(pwd, 'Miras_Functions'), 'dir') == 7
        functionsdirectory = pwd;
    elseif exist('/home/Matlab/', 'dir') == 7
        functionsDirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsDirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsDirectory does not exist!');
    end
    addpath(fullfile(functionsDirectory, 'Miras_Functions')); 
                                            % for print_cellstr.m
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%