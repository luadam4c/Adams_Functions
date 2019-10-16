function handles = m3ha_plot_simulated_traces (varargin)
%% Plots simulated traces from NEURON output files
% Usage: handles = m3ha_plot_simulated_traces (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       handles     - TODO: Description of handles
%                   specified as a TODO
%
% Arguments:
%       varargin    - 'Directory': the directory to search in
%                   must be a string scalar or a character vector
%                   default == set in all_files.m
%                   - 'FileNames': paths to simulated data
%                   must be a string array or a cell array of character vectors
%                   default == detected from Directory
%                   - 'Extension': data file extension
%                   must be a string scalar or a character vector
%                   default == 'out'
%                   - 'XLimits': limits of x axis
%                               suppress by setting value to 'suppress'
%                   must be 'suppress' or a 2-element increasing numeric vector
%                   default == no restrictions
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == same as Directory
%                   - Any other parameter-value pair for plot_traces()
%
% Requires:
%       cd/all_files.m
%       cd/argfun.m
%       cd/construct_fullpath.m
%       cd/decide_on_colormap.m
%       cd/extract_columns.m
%       cd/load_neuron_outputs.m
%       cd/plot_traces.m
%       cd/set_figure_properties.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-10-14 Created by Adam Lu
% 

%% Hard-coded parameters
% Column numbers for simulated data
%   Note: Must be consistent with singleneuron4compgabab.hoc
TIME_COL_SIM = 1;
VOLT_COL_SIM = 2;
DEND1_COL_SIM = 3;
DEND2_COL_SIM = 4;
IDCLAMP_COL_SIM = 5;
GGABAB_COL_SIM = 6;
iCP_COL_SIM = 7;
IEXT_COL_SIM = 8;
ICA_COL_SIM = 9;
ITM_COL_SIM = 10;
ITMINF_COL_SIM = 11;
ITH_COL_SIM = 12;
ITHINF_COL_SIM = 13;
IH_COL_SIM = 14;
IHM_COL_SIM = 15;
IKA_COL_SIM = 16;
IAM1_COL_SIM = 17;
IAH1_COL_SIM = 18;
IAM2_COL_SIM = 19;
IAH2_COL_SIM = 20;
IKKIR_COL_SIM = 21;
IKIRM_COL_SIM = 22;
INAPNA_COL_SIM = 23;
INAPM_COL_SIM = 24;
INAPH_COL_SIM = 25;

% TODO: Use this
labelsAll = {'Time (ms)'; 'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
        'V_{dend2} (mV)'; 'I_{GABA_B} (nA)'; 'g_{GABA_B} (uS)'; ...
        'I_{cp} (nA)'; 'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
        'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
        'I_{h} (mA/cm^2)'; 'm_{h}'; 'I_{A} (mA/cm^2)'; ...
        'm_{1,A}'; 'h_{1,A}'; 'm_{2,A}'; 'h_{2,A}'; ...
        'I_{Kir} (mA/cm^2)'; 'm_{\infty,Kir}'; ...
        'I_{NaP} (mA/cm^2)'; 'm_{\infty,NaP}'; 'h_{NaP}'};

%% Default values for optional arguments
directoryDefault = '';          % set in all_files.m
fileNamesDefault = {};
extensionDefault = 'out';       % 
% xLimitsDefault = [];            % set later
xLimitsDefault = [2800, 4500];
outFolderDefault = '';          % set later
%TODO
expStr = '';
colorMap = [];
simMode = 'active';
residuals = [];
vVecsRec = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['fileNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Extension', extensionDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) isempty(x) || iscell(x) || ischar(x) && strcmpi(x, 'suppress') || ...
        isnumeric(x) && isvector(x) && length(x) == 2);
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
extension = iP.Results.Extension;
xLimits = iP.Results.XLimits;
outFolder = iP.Results.OutFolder;

% Keep unmatched arguments for the plot_traces() function
otherArguments = iP.Unmatched;

%% Preparation
% Use the present working directory for both inputs and output by default
if isempty(directory)
    directory = pwd;
end

% Set default output directory
if isempty(outFolder)
    outFolder = directory;
end

% Decide on input paths
if isempty(fileNames)
    [~, fileNames] = ...
        all_files('Directory', directory, 'Extension', extension);
else
    % Make sure they are fulle paths
    fileNames = construct_fullpath(fileNames, 'Directory', directory);
end

% Use the common prefix as the experiment string
if isempty(expStr)
    expStr = extract_fileparts(fileNames, 'commonprefix');
end

% Create an experiment identifier for title
expStrForTitle = strrep(expStr, '_', '\_');

% Count the number of files
nFiles = numel(fileNames);

% Decide on the colors for each row in the plots
colorMap = decide_on_colormap(colorMap, nFiles);

%% Do the job
% Load simulated data
simData = load_neuron_outputs('FileNames', fileNames);
% simDataOrig = load_neuron_outputs('FileNames', fileNames);

% If recorded data provided (tVecs not empty at this point), 
%   interpolate simulated data to match the time points of recorded data
% Note: This is necessary because CVODE (variable time step method) 
%       is applied in NEURON
% TODO
% if ~isempty(tVecs)
%     simData = cellfun(@(x, y) match_time_points(x, y), ...
%                         simDataOrig, tVecs, 'UniformOutput', false);
% else
%     simData = simDataOrig;
% end

% Extract vectors from simulated data
%   Note: these are arrays with 25 columns
if strcmpi(simMode, 'passive')
    [tVecs, vVecsSim, iVecsSim, vVecsDend1, vVecsDend2] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        IEXT_COL_SIM, DEND1_COL_SIM, DEND2_COL_SIM]);
elseif strcmpi(simMode, 'active')
    [tVecs, vVecsSim, gVecsSim, iVecsSim, icaVecsSim, ...
            itmVecsSim, itminfVecsSim, ithVecsSim, ithinfVecsSim, ...
            ihVecsSim, ihmVecsSim, ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
            iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
            inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
        extract_columns(simData, [TIME_COL_SIM, VOLT_COL_SIM, ...
                        GGABAB_COL_SIM, IEXT_COL_SIM, ...
                        ICA_COL_SIM, ITM_COL_SIM, ITMINF_COL_SIM, ...
                        ITH_COL_SIM, ITHINF_COL_SIM, ...
                        IH_COL_SIM, IHM_COL_SIM, ...
                        IKA_COL_SIM, IAM1_COL_SIM, IAH1_COL_SIM, ...
                        IAM2_COL_SIM, IAH2_COL_SIM, ...
                        IKKIR_COL_SIM, IKIRM_COL_SIM, ...
                        INAPNA_COL_SIM, INAPM_COL_SIM, INAPH_COL_SIM]);
end
% Find the indices of the x-axis limit endpoints
endPointsForPlots = find_window_endpoints(xLimits, tVecs);

% Prepare vectors for plotting
if strcmpi(simMode, 'passive')
    [tVecs, vVecsSim, vVecsRec, residuals, ...
        iVecsSim, vVecsDend1, vVecsDend2] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, vVecsSim, vVecsRec, residuals, ...
                iVecsSim, vVecsDend1, vVecsDend2);
elseif strcmpi(simMode, 'active')
    [tVecs, residuals, vVecsRec, vVecsSim, gVecsSim, iVecsSim, ...
        icaVecsSim, itmVecsSim, itminfVecsSim, ...
        ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
        ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
        iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
        inapnaVecsSim, inapmVecsSim, inaphVecsSim] = ...
        argfun(@(x) prepare_for_plotting(x, endPointsForPlots), ...
                tVecs, residuals, vVecsRec, vVecsSim, gVecsSim, iVecsSim, ...
                icaVecsSim, itmVecsSim, itminfVecsSim, ...
                ithVecsSim, ithinfVecsSim, ihVecsSim, ihmVecsSim, ...
                ikaVecsSim, iam1VecsSim, iah1VecsSim, ...
                iam2VecsSim, iah2VecsSim, ikkirVecsSim, ikirmVecsSim, ...
                inapnaVecsSim, inapmVecsSim, inaphVecsSim);
end

% Compute processed data
itm2hVecsSim = (itmVecsSim .^ 2) .* ithVecsSim;
itminf2hinfVecsSim = (itminfVecsSim .^ 2) .* ithinfVecsSim;

% Select data to plot
if strcmpi(simMode, 'passive')
    dataForOverlapped = {vVecsSim; vVecsDend1; vVecsDend2; iVecsSim};
elseif strcmpi(simMode, 'active')
    dataForOverlapped = {vVecsSim; gVecsSim; iVecsSim; ...
            icaVecsSim; itm2hVecsSim; itminf2hinfVecsSim; ...
            itmVecsSim; itminfVecsSim; ithVecsSim; ithinfVecsSim; ...
            ihVecsSim; ihmVecsSim; ...
            ikaVecsSim; iam1VecsSim; iah1VecsSim; ...
            iam2VecsSim; iah2VecsSim; ikkirVecsSim; ikirmVecsSim; ...
            inapnaVecsSim; inapmVecsSim; inaphVecsSim};
end

% Construct matching y labels
if strcmpi(simMode, 'passive')
    yLabelsOverlapped = {'V_{soma} (mV)'; 'V_{dend1} (mV)'; ...
                        'V_{dend2} (mV)'; 'I_{stim} (nA)'};
elseif strcmpi(simMode, 'active')
    yLabelsOverlapped = {'V_{soma} (mV)'; 'g_{GABA_B} (uS)'; ...
            'I_{stim} (nA)'; 'I_{Ca} (mA/cm^2)'; ...
            'm^2h_{T}'; 'm_{\infty}^2h_{\infty,T}'; ...
            'm_{T}'; 'm_{\infty,T}'; 'h_{T}'; 'h_{\infty,T}'; ...
            'I_{h} (mA/cm^2)'; 'm_{h}'; 'I_{A} (mA/cm^2)'; ...
            'm_{1,A}'; 'h_{1,A}'; 'm_{2,A}'; 'h_{2,A}'; ...
            'I_{Kir} (mA/cm^2)'; 'm_{\infty,Kir}'; ...
            'I_{NaP} (mA/cm^2)'; 'm_{\infty,NaP}'; 'h_{NaP}'};
end

% Construct matching time vectors
tVecsForOverlapped = repmat({tVecs}, size(dataForOverlapped));

% Decide on figure title and file name
figTitle = sprintf('Simulated traces for Experiment %s', expStrForTitle);
figName = fullfile(outFolder, [expStr, '_simulated.png']);

% Count the number of subplots and create figure
nSubPlots = numel(yLabelsOverlapped);
figHandle = set_figure_properties('AlwaysNew', true, ...
                'FigExpansion', [1, nSubPlots/4]);

% Plot overlapped traces
handles = ...
    plot_traces(tVecsForOverlapped, dataForOverlapped, ...
                'Verbose', false, 'PlotMode', 'parallel', ...
                'SubplotOrder', 'list', 'ColorMode', 'byTraceInPlot', ...
                'LegendLocation', 'suppress', ...
                'ColorMap', colorMap, 'XLimits', xLimits, ...
                'LinkAxesOption', 'x', 'XUnits', 'ms', ...
                'YLabel', yLabelsOverlapped, ...
                'FigTitle', figTitle, 'FigHandle', figHandle, ...
                'FigName', figName, otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecs = prepare_for_plotting(vecs, endPointsForPlots)
%% Prepare vectors for plotting

% Restrict vectors to xLimits to save time on plotting
vecs = extract_subvectors(vecs, 'Endpoints', endPointsForPlots);

% Combine vectors into matrices
vecs = force_matrix(vecs, 'AlignMethod', 'leftAdjustPad');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%