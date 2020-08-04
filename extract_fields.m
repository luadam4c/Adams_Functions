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
%       [students, isMarried] = extract_fields({blab, blab}, {'students', 'isMarried'})
%       [a, b] = extract_fields(myStructArray)
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
%       cd/create_plot_movie.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/m3ha_find_decision_point.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_plot_violin.m
%       cd/parse_spike2_mat.m
%       cd/plot_correlation_coefficient.m
%       cd/update_figure_for_corel.m

% File History:
% 2019-09-03 Created by Adam Lu
% 2019-09-06 Updated so that [] is returned instead of NaN 
%               if UniformOutput is false
% 2019-12-30 Now allows the first argument to be objects or tables
% 2020-01-02 Changed the default UniformOutput to true whenever possible
% 2020-04-20 Now applies default UniformOutput to each field separately
% 2020-07-07 Fixed the case when structs are in cell arrays
% 2020-08-04 Now returns a field from a scalar struct 
%               without cell arrays by default
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
varargout = cellfun(@(f) extract_one_field(structs, f, uniformOutput), ...
                    fieldNames, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fieldValues = extract_one_field (structs, fieldName, uniformOutput)

if isempty(uniformOutput)
    try
        fieldValues = extract_one_field(structs, fieldName, true);
    catch
        fieldValues = extract_one_field(structs, fieldName, false);
        if numel(structs) == 1
            fieldValues = fieldValues{1};
        end
    end
    return
else
    % Extract field values from each structure
    %   Note: Don't use array_fun.m here
    fieldValues = arrayfun(@(x) get_field(x, fieldName, uniformOutput), ...
                            structs, 'UniformOutput', uniformOutput);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function field = get_field (myStruct, fieldName, returnNan)

% Remove from cell array
if iscell(myStruct)
    myStruct = myStruct{1};
end

% Get the field or return NaN
if is_field(myStruct, fieldName)
    field = myStruct.(fieldName);
else
    if returnNan
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