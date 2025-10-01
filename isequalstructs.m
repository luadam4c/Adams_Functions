function areEqual = isequalstructs (varargin)
%% Compares two structures or a cell array of structures for equality.
% Usage: areEqual = isequalstructs (struct1, struct2, 'IgnoreFields', {'field1'}, ...)
%        areEqual = isequalstructs (structsCellArray, 'IgnoreFields', {'field1'}, ...)
%
% Explanation:
%       This function determines if two structures are equal or if all
%       structures within a cell array are equal. It allows for specified
%       fields to be ignored during the comparison. Any parameter-value
%       pairs not recognized by this function are passed directly to the
%       built-in isequal() function, allowing for advanced comparisons
%       (e.g., isequal(s1, s2, 'nan')).
%
% Example(s):
%       s1 = struct('a', 1, 'b', [1 2]);
%       s2 = struct('a', 1, 'b', [1 2]);
%       s3 = struct('a', 1, 'b', [3 4]);
%       s4 = struct('a', 99, 'b', [3 4]);
%
%       areEqual = isequalstructs(s1, s2)
%       % returns true
%
%       areEqual = isequalstructs({s1, s2, s3})
%       % returns false
%
%       areEqual = isequalstructs({s3, s4}, 'IgnoreFields', 'a')
%       % returns true
%
% Outputs:
%       areEqual    - A logical scalar (true or false) indicating whether
%                   the structures are considered equal.
%
% Arguments:
%       structs1    - The first structure to compare, or a cell array
%                   containing all structures to be compared against each other.
%                   must be a struct or a cell array of structs
%       structs2    - (opt) The second structure to compare. This argument is
%                   required if 'structs1' is a structure. It is ignored if
%                   'structs1' is a cell array.
%                   must be a struct
%       varargin    - 'IgnoreFields': A field name or a cell array of field
%                   names to ignore during the comparison.
%                   must be a character vector, a string, or a cell array
%                       of character vectors/strings
%                   default == {}
%                   - Any other parameter-value pair for isequal()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       \Shared\Code\vIRt-Moore\virt_moore.m

% File History:
% 2025-09-30 Created by Gemini from Adam's prompt
%

%% Default values for optional arguments
ignoreFieldsDefault = {};               % default fields to ignore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Get the first argument
structs1 = varargin{1};
pvArgs = {}; % Initialize parameter-value arguments

% Determine the usage mode based on input types
if iscell(structs1)
    % Mode 1: isequalstructs(cellArray, 'param', value, ...)
    if ~all(cellfun(@isstruct, structs1))
        error('When the first argument is a cell array, it must only contain structures.');
    end
    structs2 = []; % Not used in this mode
    pvArgs = varargin(2:end);
elseif isstruct(structs1)
    % Mode 2: isequalstructs(struct1, struct2, 'param', value, ...)
    if nargin < 2 || ~isstruct(varargin{2})
        error('When the first argument is a structure, the second argument must also be a structure.');
    end
    structs2 = varargin{2};
    pvArgs = varargin(3:end);
else
    error('The first argument must be a structure or a cell array of structures.');
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options for isequal()

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IgnoreFields', ignoreFieldsDefault, ...
    @(x) assert(ischar(x) || isstring(x) || iscellstr(x), ...
        ['IgnoreFields must be a character vector, string, ', ...
         'or a cell array of character vectors.']));

% Read from the Input Parser
parse(iP, pvArgs{:});
ignoreFields = iP.Results.IgnoreFields;

% Keep unmatched arguments for the isequal() function
otherArguments = struct2arglist(iP.Unmatched);

% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs will be passed to isequal(): \n');
    disp(iP.Unmatched);
end

% Standardize ignoreFields to be a cell array
if ~iscell(ignoreFields)
    ignoreFields = {ignoreFields};
end

%% Do the job
areEqual = true; % Assume equality until proven otherwise

if iscell(structs1)
    % Mode 1: Compare all structs in a cell array
    if numel(structs1) < 2
        % A single struct or an empty cell array is trivially "equal"
        return;
    end

    % Get the first struct as the reference
    refStruct = structs1{1};
    
    % Remove ignored fields from the reference struct, if any
    if ~isempty(ignoreFields)
        validIgnoreFields = intersect(ignoreFields, fieldnames(refStruct));
        refStruct = rmfield(refStruct, validIgnoreFields);
    end

    % Loop through the rest of the structs and compare against the reference
    for i = 2:numel(structs1)
        currentStruct = structs1{i};
        
        % Remove ignored fields from the current struct
        if ~isempty(ignoreFields)
            validIgnoreFields = intersect(ignoreFields, fieldnames(currentStruct));
            currentStruct = rmfield(currentStruct, validIgnoreFields);
        end

        % Compare the structs, passing through any extra arguments
        if ~isequal(refStruct, currentStruct, otherArguments{:})
            areEqual = false;
            return; % Exit early on first mismatch
        end
    end
else
    % Mode 2: Compare two provided structs
    s1_mod = structs1;
    s2_mod = structs2;

    % Remove ignored fields from both structs
    if ~isempty(ignoreFields)
        validIgnoreFields1 = intersect(ignoreFields, fieldnames(s1_mod));
        if ~isempty(validIgnoreFields1)
            s1_mod = rmfield(s1_mod, validIgnoreFields1);
        end

        validIgnoreFields2 = intersect(ignoreFields, fieldnames(s2_mod));
        if ~isempty(validIgnoreFields2)
            s2_mod = rmfield(s2_mod, validIgnoreFields2);
        end
    end

    % Compare the modified structs, passing through any extra arguments
    areEqual = isequal(s1_mod, s2_mod, otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%