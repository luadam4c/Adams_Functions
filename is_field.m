function isField = is_field (genStruct, candName, varargin)
%% Tests whether a name is a field (in the general sense) in a structure/table/object
% Usage: isField = is_field (genStruct, candName, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       isField     - whether candName is a field in 
%                       the general sense in genStruct
%                   specified as a logical scalar
%
% Arguments:
%       genStruct   - structure in the general sense
%                       (a structure, a table or a MATLAB object)
%                   must be a struct/table/object array
%       candName    - candidate field name
%                   must be a character vector or a string scalar 
%       varargin    - Any other parameter-value pair for is_var_in_table()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/is_var_in_table.m
%
% Used by:
%       cd/extract_fields.m
%       cd/first_matching_field.m

% File History:
% 2019-12-30 Moved from first_matching_field.m
% 

%% Hard-coded parameters

%% Default values for optional arguments

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'genStruct');
addRequired(iP, 'candName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, genStruct, candName, varargin{:});

% Keep unmatched arguments for the is_var_in_table() function
otherArguments = iP.Unmatched;

%% Do the job
if isstruct(genStruct)
    isField = isfield(genStruct, candName);
elseif istable(genStruct)
    isField = is_var_in_table(candName, genStruct, otherArguments);
elseif isobject(myStruct)
    isField = isprop(genStruct, candName);
else
    error('myStruct type unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%