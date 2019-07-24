function fileBases = all_file_bases (infolder, subdirs, fileext)
%% Finds base names for files in infolder/subdir and return as a cell array or cell arrays of strings
% Usage: fileBases = all_file_bases (infolder, subdirs, fileext)
% Arguments:
%       infolder    - the directory that contains the subdirectories which contains the special case files
%       subdirs     - a cell array of subdirectory names within infolder
%       fileext     - (opt) file extension to be restricted to, a character array, e.g., 'png'
%                   default == ''
%
% Used by:    
%       /home/Matlab/Adams_Functions/find_LTS.m
%
% 2017-02-16 Created
% 2018-11-27 Renamed find_filebases -> all_file_bases
% TODO: Use all_files.m?
% TODO: Input parser

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
    error('infolder & subdirs are required, type ''help %s'' for usage', mfilename);
elseif isempty(infolder) || ~ischar(infolder)
    error('infolder must be a character array!, type ''help %s'' for usage', mfilename);
elseif isempty(subdirs) || ~iscell(subdirs)
    error('subdirs must be a cell array!, type ''help %s'' for usage', mfilename);
elseif nargin >= 3 && ~ischar(fileext)
    error('fileext must be a character array!, type ''help %s'' for usage', mfilename);
end

%% Set default arguments
if nargin < 3
    fileext = '';
end

%% Find all filebases to override
ndirs = numel(subdirs);
fileBases = cell(1, ndirs);
for k = 1:ndirs
    if ~isempty(fileext);
        fileList = dir(fullfile(infolder, subdirs{k}, ['*.', fileext]));
    else
        fileList = dir(fullfile(infolder, subdirs{k}));
    end
    ct = 0;
    for file = fileList'
        % Extract base name
        ct = ct + 1;
        [~, fileBases{k}{ct}, ~] = fileparts(file.name);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
