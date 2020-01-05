function varargout = extract_fields (structs, varargin)
%% Extracts field(s) from an array of structures/tables/objects or a cell array of structures/tables/objects
% Usage: varargout = extract_fields (structs, fieldNames (opt), varargin)
% Explanation:
%       TODO
%       See also:
%           cd/first_matching_field.m
%
% Example(s):
%       load_examples;
%       isMarried = extract_fields({blab, blab}, 'isMarried')
%       [myLogicalScalars, myNumericScalars] = extract_fields(myStructArray, {'myLogicalScalar', 'myNumericScalar'})
%       [myLogicalScalars, myNumericScalars] = extract_fields(myStructArray, {'myLogicalScalar', 'myNumericScalar'}, 'UniformOutput', false)
%
% Outputs:
%       varargout   - extracted field #1s, field #2s, etc.
%                       or extracted fields for each structure 
%                           if 'OutputMode' is 'single' TODO
%                   specified as TODO
%
% Arguments:
%       genStructs  - structures in the general sense
%                       (structures/tables/objects) to extract from
%                   must be a struct/table/object array or a cell array
%       fieldNames  - (opt) name(s) of field(s) to extract
%                   must be a character vector, a string array 
%                       or a cell array of character vectors
%                   default == all fields
%       varargin    - 'UniformOutput': whether the output is not a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true whenever possible
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/is_field.m
%
% Used by:
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_plot_violin.m
%       cd/parse_spike2_mat.m

% File History:
% 2019-09-03 Created by Adam Lu
% 2019-09-06 Updated so that [] is returned instead of NaN 
%               if UniformOutput is false
% 2019-12-30 Now allows the first argument to be objects or tables
% 2020-01-02 Changed the default UniformOutput to true whenever possible
% TODO: accept substrings
% TODO: OutputMode
% 

%% Hard-coded parameters

%% Default values for optional arguments
fieldNamesDefault = '';             % set later
uniformOutputDefault = [];          % set later

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
addRequired(iP, 'structs');

% Add optional inputs to the Input Parser
addOptional(iP, 'fieldNames', fieldNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fieldNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'UniformOutput', uniformOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, structs, varargin{:});
fieldNames = iP.Results.fieldNames;
uniformOutput = iP.Results.UniformOutput;

%% Preparation
if isempty(fieldNames)
    if iscell(structs)
        fieldNames = field_names(structs{1});
    else
        fieldNames = field_names(structs(1));
    end
else
    fieldNames = force_column_cell(fieldNames);
end

%% Do the job
if isempty(uniformOutput)
    try
        [varargout{1:nargout}] = ...
            extract_fields(structs, fieldNames, 'UniformOutput', true);
    catch
        [varargout{1:nargout}] = ...
            extract_fields(structs, fieldNames, 'UniformOutput', false);
    end
else
    varargout = cellfun(@(y) array_fun(@(x) get_field(x, y, uniformOutput), ...
                                structs, 'UniformOutput', uniformOutput), ...
                        fieldNames, 'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = get_field (myStruct, fieldName, uniformOutput)

% Get the field or return NaN
if is_field(myStruct, fieldName)
    field = myStruct.(fieldName);
else
    if uniformOutput
        field = NaN;
    else
        field = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fieldNames = field_names (genStruct)
% TODO: Pull out as its own function

if isstruct(genStruct) || isobject(genStruct)
    fieldNames = fieldnames(genStruct);
elseif istable(genStruct)
    fieldNames = genStruct.Properties.VariableNames;
else
    error('myStruct type unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%