function [fullfilename] = construct_fullfilename (filename, varargin)
%% Constructs full file name based on filename and an optional full directory path with optional suffices and/or Name-Value pairs
% Usage: [fullfilename] = construct_fullfilename (filename, varargin)
% Outputs:
%       fullfilename    - the full path to file constructed
% Arguments:
%       filename    - filename with file extension, 
%                       with or without directory, 
%                       e.g. 'A100110_0008_18.mat'
%                   must be a string scalar or a character vector
%       varargin    - 'Directory': a full directory path, 
%                       e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'Suffices': suffix(ces) to add to fileBase
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%                   default == 'nosuffices'
%                   - 'NameValuePairs': Name-Value pairs that are changed
%                   must be a 2-element cell array whose first element 
%                       is a string/char array or cell array 
%                       and whose second element is a numeric array
%                   default == {'nopairs', NaN}
%        
%
% Requires:
%       cd/construct_suffix.m
%
% Used by:
%       cd/check_dir.m
%       /media/adamX/RTCl/neuronlaunch.m
% 
% 2017-03-27 Created
% 2017-05-04 Removed ntrials
% 2017-05-04 Added input Parser scheme
% 2017-05-04 Changed paramname, paramvalue -> 'NameValuePairs' 
%                and made it a parameter 
% 2017-05-04 Made 'suffices' a parameter 
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-10-03 Now accepts file names that are full paths

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

% Add required inputs to an input Parser
addRequired(iP, 'filename', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B

% Add parameter-value pairs to the input Parser
addParameter(iP, 'Directory', '', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'Suffices', '', ...
    @(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || ...
                min(cellfun(@isstring, x))) || isstring(x), ...
                ['Suffices must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));
addParameter(iP, 'NameValuePairs', {'', NaN}, ...
    @(x) assert(iscell(x) && numel(x) == 2 ...
            && (ischar(x{1}) || iscell(x{1}) || isstring(x{1})) ...
            && isnumeric(x{2}), ...
        ['NameValuePairs must be a 2-element cell array whose ', ...
            'first element is a string/char array or cell array ', ...
            'and whose second element is a numeric array!']));

% Read from the input Parser
parse(iP, filename, varargin{:});
directory = iP.Results.Directory;
suffices = iP.Results.Suffices;
namevaluepairs = iP.Results.NameValuePairs;

%% Do the job
% Separate fileDir, fileBase and fileExt from filename
[fileDir, fileBase, fileExt] = fileparts(filename);

% Decide on the directory if not provided
if isempty(directory)
    if ~isempty(fileDir)
        % Use the original file directory
        directory = fileDir;
    else
        % Use the present working directory
        directory = pwd;
    end
end

% Construct final suffix
finalSuffix = construct_suffix('Suffices', suffices, ...
                                'NameValuePairs', namevaluepairs);

                            
% Construct full file name
if isempty(finalSuffix)                % if nothing provided
    % Construct path based on directory, fileBase and fileExt
    fullfilename = fullfile(directory, ...
                            [fileBase, fileExt]);    
elseif ~isempty(finalSuffix)              % suffix(ces) is(are) provided
    % Construct path based on directory, fileBase, final suffix and fileExt
    fullfilename = fullfile(directory, ...
                            [fileBase, '_', finalSuffix, fileExt]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
