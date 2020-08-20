function [toBeAnalyzedAll, paramsAll, errormsg] = ...
                minEASE_read_params (xlInfo, varargin)
%% Convert info from Excel file into a structure for GUI
% Usage: [toBeAnalyzedAll, paramsAll, errormsg] = ...
%               minEASE_read_params (xlInfo, varargin)
% Outputs:
%       toBeAnalyzedAll - a cell array of whether to analyze this subdirectory 
%                           ('Y' or 'N')
%       paramsAll       - a cell array of input parameter structures
% Arguments:
%       xlInfo      - a cell array containing information to read
%                       the first row must be a header
%                   must be a cell array
%       varargin    - 'MaxNumWorkers': maximum number of workers for 
%                                       parfor for initial conditions
%                                   set to 0 if parfor not to be used
%                                   set to Inf to use maximum number of workers
%                   must be a nonnegative integer or Inf
%                   default == 0
%
% Requires:
%       ~/Adams_Functions/addpath_custom
%       ~/Adams_Functions/create_error_for_nargin
%       ~/Adams_Functions/find_in_strings.m
%       ~/Adams_Functions/isaninteger.m
%       ~/Adams_Functions/locate_functionsdir
%
% Used by:
%       cd/minEASE.m

% File History:
% ---------- Created by Mark P Beenhakker
% 2017-04-26 AL - Renamed guiTOexcelConverter.m -> minEASE_read_params.m
% 2017-04-26 AL - Added xlHeader as an argument 
% 2017-04-26 AL - Renamed analyze -> toBeAnalyzed
% 2017-04-26 AL - Now uses find_in_strings to find 
%                   the corresponding index in rowInfo
% 2017-04-26 AL - strcmp to strcmpi so that 'All' will work as good as 'all'
% 2017-05-20 AL - Renamed gui_selections -> params
% 2017-07-24 AL - The parameter 'BeforeBreakMs' is changed to 'BeforePeakMs'
% 2017-08-07 AL - Now sets default to 'all' when SweepsToAnalyze not given
% 2017-10-16 AL - Added baselineWindowMs, maxBelowBasePerc and minPscDecayMs
% 2018-01-28 AL - Added isdeployed
% 2018-02-08 AL - Changed params.minPscAmp to 
%                   params.minPscAmpAbs and params.minPscAmpRel
% 2018-02-16 AL - Changed skewnessCutoff and excessKurtosisCutoff to 
%                   zSkewnessThres and zExcessKurtosisThres
% 2018-02-19 MD - TODO: Remember to log all changes here
% 2018-02-20 AL - Added errormsg as an output
% 2018-02-20 AL - Now checks and reads the entire cell array returned by xlsread
% 2018-02-20 AL - Changed arguments from rowInfo and xlHeader to xlInfo
% 2018-04-02 MD - Added checks for all inputs for Excel files
% 2018-04-02 AL - Finished checks
% 2018-07-27 AL - Changed all instances of isdir() to isfolder()
% 2018-07-27 AL - Now Sweeps To Analyze does not have to exist


%% TODO: Put a list of parameters here
%% TODO: Change all to idxXXX followed by if ~isempty(idxXXX)


%% Hard-coded parameters
inputPattern1 = '\[[0-9]*:[0-9]*\]';    % input pattern #1: 
                                        %   for Files to Analyze & 
                                        %       Sweeps to Analyze


maxNumWorkersDefault = 0;       % don't use parfor by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of arguments (better error message than inputParser)
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'xlInfo', @iscell);     % cell array with info
addParameter(iP, 'MaxNumWorkers', maxNumWorkersDefault, ...
    @(x) assert(isinf(x) || isscalar(x) && isaninteger(x) && x >= 0, ...
                'MaxNumWorkers must be Inf or a nonnegative integer!'));

% Read parameter values from the input Parser
parse(iP, xlInfo, varargin{:});
maxNumWorkers = iP.Results.MaxNumWorkers;

% Extract the header
xlHeader = xlInfo(1, :);                % a cell array containing the header
%% TODO: Check if the header is a cell array of character arrays (iscellstr())

% Extract number of rows
nRows = size(xlInfo, 1);

% Initialize output
toBeAnalyzedAll = cell(nRows, 1);
paramsAll = cell(nRows, 1);

%% Find all relevant indices in xlHeader

% Find the index of "To Analyze (Y/N)" in the Excel header
idxToBeAnalyzed = ...           % whether this row will be analyzed
    find_in_strings({'Analyze', 'Y/N'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "To Analyze (Y/N)" doesn't exist
if isempty(idxToBeAnalyzed)
    errormsg = 'A ''To Analyze (Y/N)'' column must exist!!';
    return;
end

% Find the index of "Data Home Directory" in the Excel header
idxDataHomeDirectory = ...      % data home directory
    find_in_strings({'Home', 'Directory'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Data Home Directory" doesn't exist
if isempty(idxDataHomeDirectory)
    errormsg = 'A ''Data Home Directory'' column must exist!!';
    return;
end

% Find the index of "Data Subdirectory" in the Excel header
idxDataSubdirectory = ...       % data subdirectory
    find_in_strings({'Data', 'Subdirectory'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

                    
% Return with error message if "Data Subdirectory" doesn't exist
if isempty(idxDataSubdirectory)
    errormsg = 'A ''Data Subdirectory'' column must exist!!';
    return;
end

% Find the index of "Output directory" in the Excel header
idxOutputDirectory = ...        % output directory
    find_in_strings({'Output', 'Directory'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Output directory" doesn't exist
if isempty(idxOutputDirectory)
    errormsg = 'A ''Output Subdirectory'' column must exist!!';
    return;
end

% Find the index of "Files to Analyze" in the Excel header
idxFilesToAnalyze = ...         % files in the data subdirectory to analyze
    find_in_strings({'Files', 'Analyze'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Files to Analyze" doesn't exist
if isempty(idxFilesToAnalyze)
    errormsg = 'A ''Files to Analyze'' column must exist!!';
    return;
end

% Find the index of "Sweeps to Analyze" in the Excel header
idxSweepsToAnalyze = ...        % sweeps in each file to analyze
    find_in_strings({'Sweeps', 'Analyze'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Find the index of "Direction of PSC" in the Excel header
idxDirectionPsc = ...           % PSC direction to detect ('E' or 'I')
    find_in_strings({'Direction', 'PSC'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Direction of PSC" doesn't exist
if isempty(idxDirectionPsc)
    errormsg = 'A ''Direction of PSC'' column must exist!!';
    return;
end

% Find the index of "Lowpass Filter Cutoff Frequency" in the Excel header
idxLowpassCutoff = ...          % lowpass filter cutoff frequency (Hz)
    find_in_strings({'lowpass', 'cutoff'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Lowpass Filter Cutoff Frequency" doesn't exist
if isempty(idxLowpassCutoff)
    errormsg = 'A ''Lowpass Filter Cutoff Frequency'' column must exist!!';
    return;
end

% Find the index of "Lowpass Butterworth Filter Order" in the Excel header
idxLowpassNpoles = ...          % lowpass Butterworth filter order
    find_in_strings({'lowpass', 'order'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Lowpass Butterworth Filter Order" doesn't exist
if isempty(idxLowpassNpoles)
    errormsg = 'A ''Lowpass Butterworth Filter Order'' column must exist!!';
    return;
end

% Find the index of "Noise Window Size " in the Excel header
idxNoiseWindowSize = ...        % noise window size (samples)
    find_in_strings({'noise', 'window'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Noise Window Size" doesn't exist
if isempty(idxNoiseWindowSize)
    errormsg = 'A ''Noise Window Size'' column must exist!!';
    return;
end

% Find the index of "Noise Skewness Z-score Threshold" in the Excel header
idxZSkewnessThres = ...         % Gaussian noise skewness z-score threshold
    find_in_strings({'z', 'skew', 'thres'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Noise Skewness Z-score Threshold" doesn't exist
if isempty(idxZSkewnessThres)
    errormsg = 'A ''Noise Skewness Z-score Threshold'' column must exist!!';
    return;
end

% Find the index of "Noise Excess Kurtosis Z-score Threshold" in the Excel header
idxZExcessKurtosisThres = ...   % Gaussian noise excess kurtosis z-score threshold
    find_in_strings({'z', 'kurt', 'thres'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Noise Excess Kurtosis Z-score Threshold" doesn't exist
if isempty(idxZExcessKurtosisThres)
    errormsg = 'A ''Noise Excess Kurtosis Z-score Threshold'' column must exist!!';
    return;
end

% Find the index of "Signal to Noise ratio" in the Excel header
idxSignal2Noise = ...           % signal-to-noise ratio for an event
    find_in_strings({'Signal', 'Noise'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Signal to Noise ratio" doesn't exist
if isempty(idxSignal2Noise)
    errormsg = 'A ''Signal to Noise ratio'' column must exist!!';
    return;
end

% Find the index of "Minimum Amplitude Threshold" in the Excel header
idxMinEventThreshold = ...      % minimum event amplitude threshold (pA)
    find_in_strings({'min', 'amp', 'thres'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Minimum Amplitude Threshold" doesn't exist
if isempty(idxMinEventThreshold )
    errormsg = 'A ''Minimum Amplitude Threshold'' column must exist!!';
    return;
end

% Find the index of "Baseline Window Size" in the Excel header
idxBaselineWindowMs = ...       % window size for computing baseline (ms)
    find_in_strings({'baseline', 'window'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Baseline Window Size" doesn't exist
if isempty(idxBaselineWindowMs)
    errormsg = 'A ''Baseline Window Size'' column must exist!!';
    return;
end

% Find the index of "Maximum Below Baseline Percentage" in the Excel header
idxMaxBelowBasePerc = ...       % maximum below baseline percentage (%)
    find_in_strings({'max', 'baseline', 'perc'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Maximum Below Baseline Percentage" doesn't exist
if isempty(idxMaxBelowBasePerc)
    errormsg = 'A ''Maximum Below Baseline Percentage'' column must exist!!';
    return;
end

% Find the index of "Moving Average Filter Window" in the Excel header
idxSmoothWindowMs = ...         % moving average filter window (ms)
    find_in_strings({'moving', 'window'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Moving Average Filter Window" doesn't exist
if isempty(idxSmoothWindowMs)
    errormsg = 'A ''Moving Average Filter Window'' column must exist!!';
    return;
end

% Find the index of "Minimum PSC Amplitude (pA)" in the Excel header
idxMinPscAmpAbs = ...           % minimum allowable PSC amplitude (pA)
    find_in_strings({'min', 'psc', 'amp', 'pA'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);


% Find the index of "Minimum PSC Amplitude (%)" in the Excel header
idxMinPscAmpRel = ...           % minimum allowable PSC amplitude (% of maximum)
    find_in_strings({'min', 'psc', 'amp', '%'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if neither "Minimum PSC Amplitude (pA)"
%   nor "Minimum PSC Amplitude (%)" exists
if isempty(idxMinPscAmpRel) && isempty(idxMinPscAmpAbs)
    errormsg = ['A ''Minimum PSC Amplitude (pA)'' or ', ...
               'a ''Minimum PSC Amplitude (%)'' column must exist!!'];
    return;
end

% Find the index of "Maximum PSC 10-90% Rise Time" in the Excel header
idxMaxPscRiseMs = ...           % maximum allowable PSC rise time (ms)
    find_in_strings({'max', 'rise'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Maximum PSC 10-90% Rise Time" doesn't exist
if isempty(idxMaxPscRiseMs)
    errormsg = 'A ''Maximum PSC 10-90% Rise Time'' column must exist!!';
    return;
end

% Find the index of "Minimum PSC 50% Decay Time" in the Excel header
idxMinPscDecayMs = ...          % minimum allowable PSC decay time (ms)
    find_in_strings({'min', 'decay'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Minimum PSC 50% Decay Time" doesn't exist
if isempty(idxMinPscDecayMs)
    errormsg = 'A ''Minimum PSC 50% Decay Time'' column must exist!!';
    return;
end

% Find the index of "Maximum PSC 50% Decay Time" in the Excel header
idxMaxPscDecayMs = ...          % maximum allowable PSC decay time (ms)
    find_in_strings({'max', 'decay'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Maximum PSC 50% Decay Time" doesn't exist
if isempty(idxMaxPscDecayMs)
    errormsg = 'A ''Maximum PSC 50% Decay Time'' column must exist!!';
    return;
end

% Find the index of "Seal Test Window" in the Excel header
idxSealTestWindowMs = ...       % seal test window (ms)
    find_in_strings({'seal', 'window'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Seal Test Window" doesn't exist
if isempty(idxSealTestWindowMs)
    errormsg = 'A ''Seal Test Window'' column must exist!!';
    return;
end

% Find the index of "Total PSC Trace Length" in the Excel header
idxTraceLengthMs = ...          % total trace length for each PSC (ms)
    find_in_strings({'total', 'length'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Total PSC Trace Length" doesn't exist
if isempty(idxTraceLengthMs)
    errormsg = 'A ''Total PSC Trace Length'' column must exist!!';
    return;
end

% Find the index of "PSC Trace Length Before Peak" in the Excel header
idxBeforePeakMs = ...           % trace length before peak (ms)
    find_in_strings({'before', 'peak'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "PSC Trace Length Before Peak" doesn't exist
if isempty(idxBeforePeakMs)
    errormsg = 'A ''PSC Trace Length Before Peak'' column must exist!!';
    return;
end

% Find the index of "Plot Event Detection Flag" in the Excel header
idxPlotEventDetectionFlag = ... % whether to plot event detection
    find_in_strings({'Plot', 'Event', 'Detection'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);

% Return with error message if "Plot Event Detection Flag" doesn't exist
if isempty(idxPlotEventDetectionFlag)
    errormsg = 'A ''Plot Event Detection Flag'' column must exist!!';
    return;
end

% Find the index of "Plot Average Trace Flag" in the Excel header
idxPlotAverageTraceFlag = ...   % whether to plot average trace
    find_in_strings({'Plot', 'Average', 'Trace'}, xlHeader, ...
                        'SearchMode', 'substrings', 'IgnoreCase', true);
                    
% Return with error message if "Plot Average Trace Flag" doesn't exist
if isempty(idxPlotAverageTraceFlag)
    errormsg = 'A ''Plot Average Trace Flag'' column must exist!!';
    return;
end

% Extract the parameters for each row
if maxNumWorkers == 0
    for row = 2:nRows
        % Extract information for this row (this data subdirectory)
        rowInfo = xlInfo(row, :);               % info for this data subdirectory

        % Read rowInfo into toBeAnalyzedAll and paramsAll
        [errormsg, toBeAnalyzedAll{row-1}, paramsAll{row-1}] = ...
            read_to_params(row, rowInfo, inputPattern1, ...
                            idxToBeAnalyzed, ...
                            idxDataHomeDirectory, idxDataSubdirectory, ...
                            idxOutputDirectory, idxFilesToAnalyze, ...
                            idxSweepsToAnalyze, idxDirectionPsc, ...
                            idxLowpassCutoff, idxLowpassNpoles, ...
                            idxNoiseWindowSize, idxZSkewnessThres, ...
                            idxZExcessKurtosisThres, idxSignal2Noise, ...
                            idxMinEventThreshold, idxBaselineWindowMs, ...
                            idxMaxBelowBasePerc, idxSmoothWindowMs, ...
                            idxMinPscAmpAbs, idxMinPscAmpRel, ...
                            idxMaxPscRiseMs, idxMinPscDecayMs, ...
                            idxMaxPscDecayMs, idxSealTestWindowMs, ...
                            idxTraceLengthMs, idxBeforePeakMs, ...
                            idxPlotEventDetectionFlag, idxPlotAverageTraceFlag);
            
        % If errormsg is not empty, return
        if ~isempty(errormsg)
            return
        end
    end
else
    % Preallocate error messages
    errormsg = cell(1, nRows);

    % Run parfor loop
%    for row = 2:nRows
    parfor (row = 2:nRows, maxNumWorkers)
        % Extract information for this row (this data subdirectory)
        rowInfo = xlInfo(row, :);               % info for this data subdirectory

        % Read rowInfo into toBeAnalyzedAll and paramsAll
        [errormsg{row-1}, toBeAnalyzedAll{row-1}, paramsAll{row-1}] = ...
            read_to_params(row, rowInfo, inputPattern1, ...
                            idxToBeAnalyzed, ...
                            idxDataHomeDirectory, idxDataSubdirectory, ...
                            idxOutputDirectory, idxFilesToAnalyze, ...
                            idxSweepsToAnalyze, idxDirectionPsc, ...
                            idxLowpassCutoff, idxLowpassNpoles, ...
                            idxNoiseWindowSize, idxZSkewnessThres, ...
                            idxZExcessKurtosisThres, idxSignal2Noise, ...
                            idxMinEventThreshold, idxBaselineWindowMs, ...
                            idxMaxBelowBasePerc, idxSmoothWindowMs, ...
                            idxMinPscAmpAbs, idxMinPscAmpRel, ...
                            idxMaxPscRiseMs, idxMinPscDecayMs, ...
                            idxMaxPscDecayMs, idxSealTestWindowMs, ...
                            idxTraceLengthMs, idxBeforePeakMs, ...
                            idxPlotEventDetectionFlag, idxPlotAverageTraceFlag);
    end

    % If any of the errormsgs is not empty, return with that error
    for row = 2:nRows
        if ~isempty(errormsg{row})
            errormsg = errormsg{row};
            return
        end
    end
end

% Return an empty error message if successful
errormsg = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errormsg, toBeAnalyzed, params] = read_to_params(row, rowInfo, ...
                            inputPattern1, ...
                            idxToBeAnalyzed, ...
                            idxDataHomeDirectory, idxDataSubdirectory, ...
                            idxOutputDirectory, idxFilesToAnalyze, ...
                            idxSweepsToAnalyze, idxDirectionPsc, ...
                            idxLowpassCutoff, idxLowpassNpoles, ...
                            idxNoiseWindowSize, idxZSkewnessThres, ...
                            idxZExcessKurtosisThres, idxSignal2Noise, ...
                            idxMinEventThreshold, idxBaselineWindowMs, ...
                            idxMaxBelowBasePerc, idxSmoothWindowMs, ...
                            idxMinPscAmpAbs, idxMinPscAmpRel, ...
                            idxMaxPscRiseMs, idxMinPscDecayMs, ...
                            idxMaxPscDecayMs, idxSealTestWindowMs, ...
                            idxTraceLengthMs, idxBeforePeakMs, ...
                            idxPlotEventDetectionFlag, idxPlotAverageTraceFlag)

% Initialize error message as empty
errormsg = '';

% Initialize params as an empty struct
params = struct;

% Extract whether to analyze this subdirectory (must be 'Y' or 'N')
toBeAnalyzed = rowInfo{idxToBeAnalyzed};
toBeAnalyzed = validate_string(toBeAnalyzed, {'Yes', 'No'});

if isempty(toBeAnalyzed)
    errormsg = sprintf(['The toBeAnalyzed value for ', ...
                        'row %d is invalid!!'], row);
    return;
end

%% If not to be analyzed, keep params empty and return
if strcmp(toBeAnalyzed, 'No')
    errormsg = '';
    return;
end

%% Extract input/Output Directories
% Read in the data home directory to the params structure
params.dataHomeDirectory = rowInfo{idxDataHomeDirectory};

% Check that the data home directory exists
if ~isfolder(params.dataHomeDirectory)
    errormsg = sprintf(['The data home directory for ', ...
                        'row %d does not exist!!'], row);
    return;
end

% Read in the data subdirectory to the params structure
params.dataSubdirectory = rowInfo{idxDataSubdirectory};

% Check that the data subdirectory exists
if ~isfolder(fullfile(params.dataHomeDirectory, params.dataSubdirectory))
    errormsg = sprintf(['The data subdirectory for ', ...
                        'row %d does not exist!!'], row);
    return;
end

% Read in the output directory to the params structure
params.outputDirectory = rowInfo{idxOutputDirectory};
    %Checks that input is a character array for Output Directory (does not have to exist)
if ~ischar(params.outputDirectory)
    errormsg = sprintf(['The output directory for ', ...
                        'row %d does not exist!!'], row);
    return;
end

% Read in the 'To Analyze' (Y/N) option to the params structure
params.toBeAnalyzed = rowInfo{idxToBeAnalyzed};
    %Checks that input is either 'Yes' or 'No'

%% Files to analyze
tempStr = rowInfo{idxFilesToAnalyze};   % 'all' or '[1:2]', etc.
if ischar(tempStr) && strcmpi(tempStr, 'all')
                                        % if the input is 'all'
    params.filesToAnalyze = 'all';      % set to 'all'
elseif ischar(tempStr) && regexp(tempStr, inputPattern1)
                                        % if the input is of the form [%f:%f]
    % Scan for the first and last file numbers
    temp_cellarray = textscan(tempStr, '[%f:%f]');    
    startFiles = temp_cellarray{1};     % file # of first file
    endFiles = temp_cellarray{2};       % file # of last file
    
    % Set to a vector of file numbers to analyze
    params.filesToAnalyze = startFiles:endFiles;
else
    errormsg = sprintf(['Files To Analyze should be of the form', ...
                        '[idxStart:idxEnd] (e.g. [1:10])!!']);
    return;
end

%% Sweeps to analyze (default: 'all')
if ~isempty(idxSweepsToAnalyze)
    tempStr = rowInfo{idxSweepsToAnalyze};  % 'all' or '[1:2]', etc.
    if ischar(tempStr) && strcmpi(tempStr, 'all')
                                            % if the input is 'all'
        params.sweepsToAnalyze = 'all';     % set to 'all'
    elseif ischar(tempStr) && regexp(tempStr, inputPattern1)
                                            % if the input is of the form [%f:%f]
        % Scan for the first and last file numbers
        temp_cellarray = textscan(tempStr, '[%f:%f]');    
        startFiles = temp_cellarray{1};     % file # of first file
        endFiles = temp_cellarray{2};       % file # of last file
        
        % Set to a vector of file numbers to analyze
        params.sweepsToAnalyze = startFiles:endFiles;
    else
        errormsg = sprintf(['Sweeps To Analyze should be of the form', ...
                            '[idxStart:idxEnd] (e.g. [1:10])!!']);
        return;
    end
else
    params.sweepsToAnalyze = 'all';         % set to 'all'
end

%% Parameters for synaptic event detection
% Event direction
params.directionPsc = rowInfo{idxDirectionPsc};
if validate_string(params.directionPsc, {'Excitatory', 'Inhibitory', 'E', 'I'})
else 
    errormsg = sprintf(['Direction of PSC is not specified as ', ...
                        'either ''Excitatory'' or ''Inhibitory''']);
    return;
end

% Parameters for lowpass filtering
params.lowpassCutoff = rowInfo{idxLowpassCutoff};
if isempty(params.lowpassCutoff) || ...
    ~(isnumeric(params.lowpassCutoff) && all(params.lowpassCutoff > 0))
    errormsg = sprintf(['Lowpass filter cutoff frequency not specified ', ...
                        'or is not a positive number!']);
    return;
end

params.lowpassNpoles = rowInfo{idxLowpassNpoles};
if isempty (params.lowpassNpoles) || ...
    ~(isnumeric(params.lowpassNpoles) && all(params.lowpassNpoles > 0) ...
    && floor(params.lowpassNpoles) == params.lowpassNpoles)
    errormsg = sprintf(['Lowpass Butterworth filter order not specified ', ...
                        'or is not a positive integer!']);
    return;
end

% Parameters for computing Gaussian noise level
params.noiseWindowSize = rowInfo{idxNoiseWindowSize};
if isempty(params.noiseWindowSize) || ~(isnumeric(params.noiseWindowSize)...
     && all(params.noiseWindowSize > 0) ...
     && floor(params.noiseWindowSize) == params.noiseWindowSize)
     errormsg = sprintf(['Noise window size not specified ', ...
                        'or is not a positive integer!']);
    return;
end

params.zSkewnessThres = rowInfo{idxZSkewnessThres};
if isempty(params.zSkewnessThres) || ~(isnumeric(params.zSkewnessThres) ...
    && all(params.zSkewnessThres >= 0))
     errormsg = sprintf(['Noise skewness Z-score threshold not specified ', ...
                        'or is not a nonnegative number!']);
    return;
end

params.zExcessKurtosisThres = rowInfo{idxZExcessKurtosisThres};
if isempty(params.zExcessKurtosisThres) || ...
    ~(isnumeric(params.zExcessKurtosisThres) ...
    && all(params.zExcessKurtosisThres >= 0))
     errormsg = sprintf(['Noise excess Kurtosis Z-score threshold ', ...
                        'not specified or is not a nonnegative number!']);
    return;
end

% Parameters for detecting directional events
params.signal2Noise = rowInfo{idxSignal2Noise};    
if isempty(params.signal2Noise) || ~(isnumeric(params.signal2Noise) ...
    && all(params.signal2Noise > 0))
     errormsg = sprintf(['Signal to noise ratio not specified ', ...
                        'or is not a positive number!']);
    return;
end

params.minEventThreshold = rowInfo{idxMinEventThreshold};    
if isempty(params.minEventThreshold) || ...
    ~(isnumeric(params.minEventThreshold) && all(params.minEventThreshold >= 0))
     errormsg = sprintf(['Minimum threshold amplitude for event not ', ...
                        'specified or is not a nonnegative number!']);
    return;
end

params.baselineWindowMs = rowInfo{idxBaselineWindowMs};       
if isempty(params.baselineWindowMs) || ~(isnumeric(params.baselineWindowMs) ...
    && all(params.baselineWindowMs > 0))
     errormsg = sprintf(['Baseline window size not specified ', ...
                        ' or is not a positive number!']);
    return;
end

params.maxBelowBasePerc = rowInfo{idxMaxBelowBasePerc};    
if isempty (params.maxBelowBasePerc) || ~(isnumeric(params.maxBelowBasePerc) ...
    && params.maxBelowBasePerc >= 0 && params.maxBelowBasePerc <= 100)
     errormsg = sprintf(['Maximum below baseline percentage not specified ', ...
                        'or is not a value between 0 and 100!']);
    return;
end

params.smoothWindowMs = rowInfo{idxSmoothWindowMs};
if isempty(params.smoothWindowMs) || ~(isnumeric(params.smoothWindowMs) ...
    && params.smoothWindowMs > 0)
    errormsg = sprintf(['Moving average filter window not specified ', ...
                        'or is not a positive number!']);
    return;
end

% Parameters for classifying PSCs
if ~isempty(idxMinPscAmpAbs)
    params.minPscAmpAbs = rowInfo{idxMinPscAmpAbs};
elseif ~isempty(idxMinPscAmpRel)
    params.minPscAmpRel = rowInfo{idxMinPscAmpRel};
else
    error('Minimum PSC amplitude missing!');
end

params.maxPscRiseMs = rowInfo{idxMaxPscRiseMs};
if isempty(params.maxPscRiseMs) || ~(isnumeric(params.maxPscRiseMs) ...
    && params.maxPscRiseMs >= 0)
    errormsg = sprintf(['Maximum PSC amplitude not specified ', ...
                        'or is not a nonnegative number']);
    return;
end

params.minPscDecayMs = rowInfo{idxMinPscDecayMs};
if isempty(params.minPscDecayMs) || ~(isnumeric(params.minPscDecayMs) ...
    && params.minPscDecayMs >= 0)
    errormsg = sprintf(['Minimum PSC 50-percent decay time not specified ', ...
                        'or is not a nonnegative number!']);
    return;
end

params.maxPscDecayMs = rowInfo{idxMaxPscDecayMs};
if isempty(params.maxPscDecayMs) || ~(isnumeric(params.maxPscDecayMs) ...
    && params.maxPscDecayMs >= 0)
    errormsg = sprintf(['Maximum PSC 50-percent decay time not specified ', ...
                        'or is not a nonnegative number!']);
    return;
end

tempStr2 = ...                          % seal test window (ms) '[1000, 1050]'
    rowInfo{idxSealTestWindowMs};
if isempty(tempStr2)
    errormsg = sprintf('Seal test window not specified!');
    return;
elseif ~ischar(tempStr2)
    errormsg = sprintf('Seal test window is not a character array!');
    return;
end
temp_cellarray = textscan(tempStr2, '[%f, %f]');    % scan for the first and 
                                                    %   last numbers
params.sealTestWindowMs = [temp_cellarray{1}, temp_cellarray{2}];
if isempty(params.sealTestWindowMs) ...
   || ~(isnumeric(params.sealTestWindowMs) ...
   && all(params.sealTestWindowMs >= 0) ...
   && numel(params.sealTestWindowMs) == 2)
    errormsg = sprintf(['Seal test window not specified ', ...
                        'or is not a nonnegative number with two inputs!']);
    return;
end

% Parameters for computing average PSC trace
params.traceLengthMs = rowInfo{idxTraceLengthMs};
if isempty(params.traceLengthMs) || ~(isnumeric(params.traceLengthMs) ...
    && params.traceLengthMs >= 0)
    errormsg = sprintf(['Total PSC trace length is not specified or is ', ...
                        'not a nonnegative number!']);
    return;
end

params.beforePeakMs = rowInfo{idxBeforePeakMs};
if isempty(params.beforePeakMs) || ~(isnumeric(params.beforePeakMs) ...
    && params.beforePeakMs >= 0)
    errormsg = sprintf(['PSC trace length before peak not specified or is ', ...
                        'not a nonnegative number!']);
    return;
end

% Flags
params.plotEventDetectionFlag =  rowInfo{idxPlotEventDetectionFlag};
 if isempty(params.plotEventDetectionFlag) || ...
    ~islogical(params.plotEventDetectionFlag)
    errormsg = sprintf(['Plot event detection flag not defined ', ...
                        'or is not a logical input!']);
    return;
end

params.plotAverageTraceFlag = rowInfo{idxPlotAverageTraceFlag};
if isempty(params.plotAverageTraceFlag) || ...
    ~islogical(params.plotAverageTraceFlag)
    errormsg = sprintf(['Plot average trace flag not defined ', ...
                        'or is not a logical input!']);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% save('/home/mark/matlab_temp_variables/guiTOexcelTEMP')
% lo%Checks that input is a character array for Data Home Directory
assert(ischar(params.dataHomeDirectory...
    ),'Invalid input detected in ''Data Home Directory'' of Excel file.');
ad('/home/mark/matlab_temp_variables/guiTOexcelTEMP')
% gui_selectionsOLD = load('/media/shareX/share/Mark/mini_test/guis') ;

colonFiles = findstr(':', rowInfo{6});
    %% If not to be analyzed, keep paramsAll{row-1} empty and continue
startFiles = str2num(rowInfo{6}(2:colonFiles-1));
endFiles = str2num(rowInfo{6}(colonFiles+1:end-1));

colonCap = findstr(':', rowInfo{19}) ;
startCap = str2num(rowInfo{19}(2:colonCap-1)) ;
endCap = str2num(rowInfo{19}(colonCap+1:end-1)) ;

elseif isempty(rowInfo) || isempty(xlHeader)
    error('First two inputs cannot be empty!');
elseif ~iscell(rowInfo)
    error('First input must be a cell array!');
elseif ~iscell(xlHeader)
    error('Second input must be a cell array!');

% gui_selections.findSeals = rowInfo{19};

params.minBaselineDiff = ...            % minimum baseline difference (pA)
    rowInfo{find_in_strings({'baseline', 'diff'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.crudeRegionSize = ...            % burst region size (events)
    rowInfo{find_in_strings({'region', 'size'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.minSpikesPerBurst = ...          % minimum number of spikes in a burst
    rowInfo{find_in_strings({'spikes', 'burst'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.maxIsiMs = ...                   % maximum inter-spike interval (ms) 
    rowInfo{find_in_strings({'spike', 'interval'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.startDet = ...                   % start time for detection (sec)        %%%TODO: Examine
    rowInfo{find_in_strings({'start', 'detect'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.endDet = ...                     % end time for detection (sec or "end") %%%TODO: Examine
    rowInfo{find_in_strings({'end', 'detect'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();
params.beforeBreakMs = ...              % trace length before breakpoint (ms)
    rowInfo{find_in_strings({'before', 'break'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
    %%% TODO: check input with assert() or validateattributes() or validatestring();

params.minPscAmp = ...                  % minimum allowable PSC amplitude (pA)
    rowInfo{find_in_strings({'min', 'psc', 'amp'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
function [toBeAnalyzed, params, errormsg] = minEASE_read_params (rowInfo, xlHeader)

params.skewnessCutoff = ...             % Gaussian noise skewness cutoff
    rowInfo{find_in_strings({'skew', 'cut'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};
params.excessKurtosisCutoff = ...       % Gaussian noise excess kurtosis cutoff
    rowInfo{find_in_strings({'kurt', 'cut'}, xlHeader, ...
            'SearchMode', 'substrings', 'IgnoreCase', true)};

function [toBeAnalyzed, params, errormsg] = minEASE_read_params (rowInfo, xlHeader)

%% If not to be analyzed, return empty structure for params
if toBeAnalyzed == 'N'
    params = struct;
    return;
end

addRequired(iP, 'xlHeader', @iscell);   % cell array with header

    %% If not to be analyzed, return empty structure for params

if nargin < 2

        'Invalid input detected in ''Data Subdirectory'' of Excel file.');

    %% If not to be analyzed, keep paramsAll{row-1} empty and continue
    if toBeAnalyzedAll{row-1} == 'N'
        continue
    end

% Extract whether to analyze this subdirectory (must be 'Y' or 'N')
toBeAnalyzedAll{row-1} = rowInfo{idxToBeAnalyzed};
toBeAnalyzedAll{row-1} = validatestring(toBeAnalyzedAll{row-1}, {'Y', 'N'});

% Only update params structure if to be analyzed
if toBeAnalyzedAll{row-1} == 'Y'
end

%Checks that input is a character array for Data Home Directory
assert(ischar(params.dataHomeDirectory...
    ),'Invalid input detected in ''Data Home Directory'' of Excel file.');

elseif ~isinteger(params.lowpassNpoles)
    floor(params.lowpassNpoles);
else ~isinteger(params.noiseWindowSize)
      floor(params.noiseWindowSize);

% Return with error message if "Sweeps to Analyze" doesn't exist
% 2018-07-27 Don't do this
if isempty(idxSweepsToAnalyze)
    errormsg = 'A ''Sweeps to Analyze'' column must exist!!';
    return;
end

% Template starts: All columns of the Excel file should be processed like this:
% Template ends


%}

