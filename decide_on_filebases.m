function fileBases = decide_on_filebases (fileStrings, varargin)
%% Create filebases if empty, or extract file bases
% Usage: fileBases = decide_on_filebases (fileStrings, nFiles (opt), varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       fileBases   - base of file name(s) (without extension)
%                   specified as a character vector, a string vector 
%                       or a cell array of character vectors
%                   default == 'unnamed_1', 'unnamed_2', ...
%
% Arguments:
%       fileStrings - strings containing file bases
%                   must be empty
%                       or a character vector, a string vector 
%                       or a cell array of character vectors
%       nFiles      - number of files
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fileparts.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m

% File History:
% 2019-11-15 Created by Adam Lu
% 

%% Hard-coded parameters
filePrefix = 'unnamed_';

%% Default values for optional arguments
nFilesDefault = 1;
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'fileStrings', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fileStrings must be a character array or a string array ', ...
            'or cell array of character arrays!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'nFiles', nFilesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, fileStrings, varargin{:});
nFiles = iP.Results.nFiles;
% param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% Do the job
if isempty(fileStrings)
    fileBases = create_labels_from_numbers(transpose(1:nFiles), ...
                                            'Prefix', filePrefix);
else
    pathBases = extract_fileparts(fileStrings, 'pathbase');
    fileBases = extract_fileparts(pathBases, 'dirbase');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%