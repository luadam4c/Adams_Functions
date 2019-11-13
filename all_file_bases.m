function fileBases = all_file_bases (varargin)
%% Returns all the file bases in a given directory (optionally recursive) that matches a prefix, keyword, suffix or extension
% Usage: fileBases = all_file_bases (varargin)
% Explanation:
%       This function is the same as all_files.m but returns just the 
%           file base names
%
% Example(s):
%       fileBases = all_file_bases;
%       fileBases = all_file_bases('Recursive', true);
%
% Outputs:
%       fileBases   - base names of files found
%                   specified as a cell array of character arrays
%
% Arguments:
%       varargin    - see all_files.m
%
% Requires:
%       cd/all_files.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%
% Used by:    
%       cd/find_LTS.m

% File History:
% 2017-02-16 Created
% 2018-11-27 Renamed find_filebases -> all_file_bases
% 2019-01-14 Now uses all_files.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Read from the Input Parser
parse(iP, varargin{:});

% Keep unmatched arguments for all_files.m
otherArguments = iP.Unmatched;

%% Do the job
% Find all files
files = all_files(varargin{:});

% Extract file base names
fileNames = transpose({files.name});

% Remove the extensions
fileBases = extract_fileparts(fileNames, 'base');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
