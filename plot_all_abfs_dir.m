function plot_all_abfs_dir (varargin)
%% Plots all abf files in a directory
% Usage: plot_all_abfs_dir (varargin)
%
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the abf files, e.g. '20161216'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'ExpMode': experiment mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EEG'   - EEG data; x axis in seconds; y-axis in uV
%                       'patch' - patch data; x axis in ms; y-axis in mV
%                   must be consistent with plot_traces_abf.m
%                   default == 'patch'
%                   - 'Individually': whether sweeps are plotted individually
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == a subdirectory named by {fileName}_traces in pwd
%                   - 'TimeUnits': units for time
%                   must be a string scalar or a character vector
%                   default == 's' for 2-data data and 'ms' for 3-data data
%                   - 'TimeStart': the start of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == 0
%                   - 'TimeEnd': the end of the time interval of interest 
%                                   (in units set by TimeUnits)
%                   must be a numeric nonnegative scalar
%                   default == tVec(end)
%
% Requires:
%       cd/plot_traces_abf.m
%       cd/plot_FI.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%       /home/Matlab/Downloaded_Functions/dirr.m
%       /home/Matlab/Brians_Functions/identify_CI.m
%
% File history: 
% 2016-09-22 - Created
% 2017-04-11 - Added expmode as arguments
% 2017-04-11 - Now uses dirr.m to find abf files in subdirectories too
% 2017-04-17 - BT - Creates F-I plot for current injection protocols
% 2017-04-19 - BT - Changed detection method to difference of sweep averages
% 2018-01-24 - Added isdeployed
% 2018-07-24 - Now uses a try catch statement
% 2018-09-17 - Added the input parser

%% Hard-coded parameters
validExpModes = {'EEG', 'patch'};

%% Default values for optional arguments
directoryDefault = '';          % set later
expModeDefault = 'patch';       % assume traces are patching data by default
individuallyDefault = false;    % plot all sweeps together by default
outFolderDefault = '';          % set later
timeUnitsDefault = '';          % set later
timeStartDefault = [];          % set later
timeEndDefault = [];            % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                            % for dirr.m, abf2load.m or abfload.m
    addpath(fullfile(functionsdirectory, '/Brians_Functions/'));        
                                            % for identify_CI.m
end

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) any(validatestring(x, validExpModes)));
addParameter(iP, 'Individually', individuallyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeUnits', timeUnitsDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TimeStart', timeStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
expMode = validatestring(iP.Results.ExpMode, validExpModes);
individually = iP.Results.Individually;
outFolder = iP.Results.OutFolder;
timeUnits = iP.Results.TimeUnits;
timeStart = iP.Results.TimeStart;
timeEnd = iP.Results.TimeEnd;

% Set dependent argument defaults
if isempty(directory)
    directory = pwd;
end

%% Find all .abf files
[~, ~, filenames] = dirr(directory, '.abf', 'name');
if isempty(filenames)
    fprintf('No abf files in current directory!\n');
    fprintf('Type ''help plot_all_abfs_dir'' for usage\n');
end
nfiles = numel(filenames);

%% Plot traces from each file using plot_traces_abf.m
parfor k = 1:nfiles
    % Plot traces from this abf file
    try 
        [data, siUs] =...
            plot_traces_abf(filenames{k}, ...
                            'ExpMode', expMode, 'Individually', individually, ...
                            'OutFolder', outFolder, 'TimeUnits', timeUnits, ...
                            'TimeStart', timeStart, 'TimeEnd', timeEnd);

        % Identify whether this is a current injection protocol
        isCI = identify_CI(data, siUs);

        % If it's a current injection protocol, 
        %   detect spikes for each sweep and make an F-I plot
        if isCI
            plot_FI(filenames{k}, data, siUs);
        end
    catch ME
        fprintf('Traces for %s cannot be plotted!\n', filenames{k});
        fprintf([ME.identifier, ': ', ME.message]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

files = dir(directory);
if strfind(files.name, '.abf')

files = dirr(directory, '.abf');
for file = files'

    if nargin >= 3
        plot_traces_abf(fullfile(directory, filenames{k}), expmode, recmode);
    else
        plot_traces_abf(fullfile(directory, filenames{k}), expmode);
    end

    % Load abf file
    abffilename_full = construct_abffilename(filenames{k});    % creates full path to abf file robustly
    if exist('abf2load', 'file') == 2
        [data, siUs] = abf2load(abffilename_full);
    elseif exist('abfload', 'file') == 2
        [data, siUs] = abfload(abffilename_full);
    end

        injection_data = current_data(12000:20000,:,:);
        if std(injection_data,0,1) < 4 & abs(max(max(injection_data)) - min(min(injection_data))) > 100
            plot_FI_bt(filenames{k}, data, siUs);
        end
    if length(current_data) > 20000            % tests sweeps if values within injection time range are consistent
        
        injection_data = squeeze(current_data(12000:20000, 1, :));    % current values within typical injection time range
        avgs_byswp = mean(injection_data, 1);                % average of sweeps
        reduction = abs(diff(diff(avgs_byswp)));            % reduces sweeps into differences between successive sweep averages        %%% TODO: this is not differences but rather differences of differences
        [~, max_swp] = max(avgs_byswp);                    % highest sweep by average
        max_swp_peaks_avg = mean(findpeaks(injection_data(:, max_swp)));    % average peak value of greatest sweep
        if reduction < maxSwpSpacing & max_swp_peaks_avg > 100 & size(data, 3) > 1    % sweep avgs should be separated by constant
            plot_FI(filenames{k}, data, siUs);
        end
    end

maxSwpSpacing = 2;

function plot_all_abfs_dir (directory, expmode)

%}
