function varargout = all_subdirs(varargin)
%% Returns all the subdirectories in a given directory
% Usage: [subDirs, fullPaths] = all_subdirs(varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [subDirs, fullPaths] = all_subdirs;
%       [subDirs, fullPaths] = all_subdirs('SortBy', 'date');
%       [subDirs, fullPaths] = all_subdirs('Recursive', true);
%
% Outputs:
%       subDirs     - a files structure for the subdirectories
%                   specified as a structure array with fields:
%                       name
%                       folder
%                       date
%                       bytes
%                       isdir
%                       datenum
%       fullPaths   - the full paths to the subdirectories
%                   specified as a column cell array of character vectors
%
% Arguments:    
%       varargin    - 'Level': level of subdirectories under directory
%                   must be a positive integer scalar
%                   default == 1
%                   - Any other parameter-value pair for all_files()
%
% Requires:
%       cd/all_files.m
%       cd/apply_over_cells.m
%       cd/rmfield_custom.m
%
% Used by:
%       cd/apply_to_all_subdirs.m
%       cd/m3ha_plot_figure07.m.m
%       cd/m3ha_simulate_population.m
%       /media/ashleyX/Recordings/analyze_recordings.m TODO: Update this
%       /home/zhongxiao/SCIDmiceLTP/Code/analyze_SCIDmiceLTP.m

% File History:
% 2018-09-27 Adapted from the web by Adam Lu
%               https://www.mathworks.com/matlabcentral/answers/
%                   166629-is-there-any-way-to-list-all-folders-only-in-the-level-directly-below-a-selected-directory
% 2018-10-04 Renamed subdirs() -> all_subdirs()
% 2019-07-22 Added 'ForceCellOutput' as an optional argument
% 2019-07-22 Now uses varargout
% 2019-12-13 Now uses all_files.m
% 2020-03-06 Added 'Level' as an optional argument

%% Default values for optional arguments
levelDefault = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Level', levelDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
level = iP.Results.Level;

% Keep unmatched arguments for the all_files() function
otherArguments = iP.Unmatched;

%% Find subdirectories
if level == 1
    % Use all_files.m with the 'SubDirInstead' option
    if nargout >= 1
        [varargout{1:nargout}] = ...
            all_files('SubDirInstead', true, otherArguments);
    end
elseif level == 2
    % Remove user's input 'Recursive' field if exists
    otherArguments = rmfield_custom(otherArguments, 'Recursive');

    % Find subdirectories one level down
    [~, subDirPathsOneLevelDown] = ...
        all_subdirs('Recursive', false, otherArguments);

    % Remove user's input 'Directory' field if exists
    otherArguments = rmfield_custom(otherArguments, 'Directory');

    % Find network subdirectories within the seed number subdirectories
    [subDirsCell, subDirPathsCell] = ...
        cellfun(@(x) all_subdirs('Directory', x, 'Recursive', false, ...
                                    otherArguments), ...
                subDirPathsOneLevelDown, 'UniformOutput', false);

    % Concatenate all the network subdirectories
    varargout{1} = apply_over_cells(@vertcat, subDirsCell);
    varargout{2} = apply_over_cells(@vertcat, subDirPathsCell);
else
    error('Not implemented yet!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%