function [output1] = plot_traces_spike2_mat (spike2Path, varargin)
%% Plots traces from a Spike2-exported .mat file
% Usage: [output1] = plot_traces_spike2_mat (spike2Path, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       spike2Path     - TODO: Description of spike2Path
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for plot_traces() TODO
%
% Requires:
%       cd/create_error_for_nargin.m
% TODO:
% cd/create_time_vectors.m
% cd/extract_fileparts.m
% cd/parse_spike2_mat.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-30 Moved from plethRO1_analyze.m
% 

%% Hard-coded constants
S_PER_MIN = 60;

%% Hard-coded parameters
channelNamesUser = {'O2'; 'Pleth 2'; 'WIC#2'; 'WIC#1'};
parseGas = false;

%% Default values for optional arguments
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'spike2Path', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, spike2Path, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% TODO

%% Load trace data
% Load data and parse gas trace if necessary
data = parse_spike2_mat(spike2Path, 'ParseGas', parseGas, ...
                        'ChannelNames', channelNamesUser);

% Extract the channel names and values
channelNames = data.channelNames;
channelValues = data.channelValues;
siSeconds = nanmean(data.siSeconds);
nSamples = nanmean(data.nSamples);

% Create a time vector
timeVec = create_time_vectors(nSamples, 'TimeUnits', 'min', ...
            'SamplingIntervalSec', siSeconds);

%% Load gas time data
% Find corresponding gas table
% TODO: Make a function that finds corresponding files easily
%       Consult code in load_matching_sheets.m
spike2PathBase = extract_fileparts(spike2Path, 'pathbase');
gasTablePath = [spike2PathBase, '_gas_pulses.csv'];

% Load the gas table
gasTable = readtable(gasTablePath);
startTimesSec = gasTable.startTime;
startTimesMin = startTimesSec / S_PER_MIN;
exampleCenter = startTimesMin(1);

% Create time window
exampleRelWindow = [-10, 20];
exampleNoSwdRelWindow = [-6.3, -6];
exampleSwdRelWindow = [6, 6.3];
exampleWindow = exampleCenter + exampleRelWindow;

% Restrict data
data = force_matrix(exampleChannelValues);
endPoints = find_window_endpoints(exampleWindow, timeVec);
dataShifted = extract_subvectors(data, 'EndPoints', endPoints);
timeVecShifted = extract_subvectors(timeVec, 'EndPoints', endPoints);
timeVecRelative = timeVecShifted - exampleCenter;

% Decide on plot channel names
plotChannelNames = {'O2'; 'Pleth'; 'EEG'; 'EMG'; 'Spect'};

% Compute the spectrogram
eegToShow = dataShifted(:, 3);
spectWinSeconds = 1;             % 1 second windows
spectWinSamples = round(spectWinSeconds / siSeconds);
spectNOverLap = spectWinSamples / 2;
spectSamplingFreqHz = 1 / siSeconds;
[spectData, freqHz, timeInstantsRelSeconds] = ...
    spectrogram(eegToShow, spectWinSamples, ...
                spectNOverLap, [], spectSamplingFreqHz);
timeInstantsMin = exampleRelWindow(1) + timeInstantsRelSeconds / S_PER_MIN;

% Create figure paths
[figPathBase1, figPathBase2, figPathBase3] = ...
    argfun(@(x) fullfile(hypoxiaExampleDir, x), ...
            hypoxiaExampleName, hypoxiaExampleNoSwdName, hypoxiaExampleSwdName);

% Plot traces
% TODO: Let plot_traces accept 'AxesHandles' and use it
% Create figure
[figure1b, ax] = create_subplots(5, 1, 'AlwaysNew', true, 'FigExpansion', [1, 0.4]);

for iPlot = 1:numel(plotChannelNames)

    % Plot the appropriate trace or map
    if iPlot == 5           
        % Plot the log spectrum
        imagesc(ax(iPlot), timeInstantsMin, freqHz, abs(spectData));

        % Set a colormap
        colorMapSpectFile = matfile('/media/adamX/Settings_Matlab/spectrogram_colormap.mat');
        colorMapSpect = colorMapSpectFile.colorMap;
        colormap(ax(iPlot), colorMapSpect);

        % Flip the Y Axis so lower frequencies are at the bottom
        set(ax(iPlot), 'YDir', 'normal');

        % Show only 1-40 Hz
        set(ax(iPlot), 'YLim', [1, 50]);
    else
        plot(ax(iPlot), timeVecRelative, dataShifted(:, iPlot), ...
                'Color', 'k');
        if iPlot == 1
            set(ax(iPlot), 'YLim', [33, 40]);
        elseif iPlot == 2
            set(ax(iPlot), 'YLim', [-1, 1]);
        elseif iPlot == 3
            % Show only -0.5-0.5 mV
            set(ax(iPlot), 'YLim', [-0.5, 0.5]);
        elseif iPlot == 4
            set(ax(iPlot), 'YLim', [-0.5, 0.5]);
        end
    end

    % Plot a vertical line
    plot_vertical_line(0, 'AxesHandle', ax(iPlot), 'Color', 'k');
end

set(ax, 'XTick', [], 'YTick', []);
for iPlot = 1:numel(plotChannelNames); box(ax(iPlot), 'off'); end
%    set(ax(1:(numel(plotChannelNames)-1)), 'XTick', []);
%    set(ax, 'TickDir', 'out');
linkaxes(ax, 'x');

% Save the figure in all different zooms
set_figure_properties('FigHandle', figure1b, 'FigExpansion', [0.5, 1], ...
                        'ExpandFromDefault', false);
xlim(exampleNoSwdRelWindow);
save_all_figtypes(figure1b, figPathBase2, figTypes);
xlim(exampleSwdRelWindow);
save_all_figtypes(figure1b, figPathBase3, figTypes);

set_figure_properties('FigHandle', figure1b, 'FigExpansion', [2, 1], ...
                        'ExpandFromDefault', false);
for iPlot = 1:numel(plotChannelNames)
    % Plot vertical shades
    axes(ax(iPlot))
    plot_window_boundaries([exampleNoSwdRelWindow, exampleSwdRelWindow], ...
                            'BoundaryType', 'verticalShade');
end
xlim(exampleRelWindow);
save_all_figtypes(figure1b, figPathBase1, figTypes);    

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%