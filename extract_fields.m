function varargout = extract_fields (structs, varargin)
%% Extracts field(s) from an array of structures or a cell array of structures
% Usage: varargout = extract_fields (structs, fieldNames (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       isMarried = extract_fields({blab, blab}, 'isMarried')
%       [myLogicalScalars, myNumericScalars] = extract_fields(myStructArray, {'myLogicalScalar', 'myNumericScalar'}, 'UniformOutput', true)
%
% Outputs:
%       varargout   - extracted field #1s, field #2s, etc.
%                       or extracted fields for each structure 
%                           if 'OutputMode' is 'single' TODO
%                   specified as TODO
%
% Arguments:
%       structs     - structures to extract from
%                   must be a struct array or a cell array
%       fieldNames  - (opt) name(s) of field(s) to extract
%                   must be a character vector, a string array 
%                       or a cell array of character vectors
%                   default == all fields
%       varargin    - 'UniformOutput': whether the output is not a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/spike2Mat2Text.m

% File History:
% 2019-09-03 Created by Adam Lu
% 2019-09-06 Updated so that [] is returned instead of NaN 
%               if UniformOutput is false
% TODO: accept substrings
% TODO: OutputMode
% 

%% Hard-coded parameters

%% Default values for optional arguments
fieldNamesDefault = '';             % set later
uniformOutputDefault = false;       % default TODO: Description of uniformOutput

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
addRequired(iP, 'structs', ...
    @(x) validateattributes(x, {'struct', 'cell'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'fieldNames', fieldNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fieldNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'uniformOutput', uniformOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, structs, varargin{:});
fieldNames = iP.Results.fieldNames;
uniformOutput = iP.Results.uniformOutput;

%% Preparation
if isempty(fieldNames)
    if iscell(structs)
        fieldNames = fieldnames(structs{1});
    else
        fieldNames = fieldnames(structs(1));
    end
else
    fieldNames = force_column_cell(fieldNames);
end

%% Do the job
if iscell(structs)
    varargout = cellfun(@(y) ...
                            cellfun(@(x) get_field(x, y, uniformOutput), ...
                                structs, 'UniformOutput', uniformOutput), ...
                        fieldNames, 'UniformOutput', false);
else
    varargout = cellfun(@(y) ...
                            arrayfun(@(x) get_field(x, y, uniformOutput), ...
                                structs, 'UniformOutput', uniformOutput), ...
                        fieldNames, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = get_field (myStruct, fieldName, uniformOutput)

% Get the field or return NaN
if isfield(myStruct, fieldName)
    field = myStruct.(fieldName);
else
    if uniformOutput
        field = NaN;
    else
        field = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%