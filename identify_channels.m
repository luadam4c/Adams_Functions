function [channelTypes, channelUnits, channelLabels] = ...
                identify_channels (abfdata, varargin)
%% Assigns voltage, current or conductance to each channel (2nd dim) in abfdata
% Usage: [channelTypes, channelUnits, channelLabels] = ...
%               identify_channels (abfdata, varargin)
% Explanation:
%       TODO: Explain strategy
%       For example, to find the index for the current channel, use:
%           idxCurrent = strcmp('Current', channelTypes);
%
% Outputs:
%       channelTypes    - type assigned to each channel, possibly:
%                           'Voltage', 'Current', 'Conductance' or 'Other'
%                       specified as a row cell array with the
%                           number of elements same as the length of the 
%                           2nd dimension of abfdata
%       channelUnits    - units assigned to each channel
%                       specified as a row cell array with the
%                           number of elements same as the length of the 
%                           2nd dimension of abfdata
%       channelLables   - labels assigned to each channel
%                       specified as a row cell array with the
%                           number of elements same as the length of the 
%                           2nd dimension of abfdata
% Arguments:
%       abfdata     - raw data matrix from abf2load
%                       Note: assumes the second dimension is the channel
%                   must be a numeric 2d or 3d array
%       varargin    - 'ExpMode': experiment mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'EEG'   - EEG data; x axis in seconds; y-axis in uV
%                       'patch' - patch data; x axis in ms; y-axis in mV
%                   default == 'EEG' for 2d data 'patch' for 3d data
% Requires:
%       cd/apply_iteratively.m
%
% Used by:
%       cd/parse_abf.m
%       TODO: use parse_abf instead for the following:
%       cd/plot_FI.m
%       cd/combine_sweeps.m
%       /home/Matlab/minEASE/minEASE.m
%
% File History:
% 2017-03-15 Created by BT, adapted from correct_patchannelLabels.m
% 2017-06-15 AL - Changed tabs -> spaces & modified margins
% 2017-06-15 AL - Changed vcc -> vig (for voltage, current, conductance)
% 2017-06-15 AL - Changed vig -> channelTypes, now a cell array
% 2018-09-17 AL - Renamed variables
% 2018-09-18 AL - Added an input parser
% 2018-09-18 AL - Now combines channelTypes and channelUnits into channelLabels
% 2018-09-18 AL - Now considers the case for EEG data
% 2018-09-19 AL - Added the case for 4 channels
% 2018-09-20 AL - Added the case for 3 channels when the 3rd channel
%                   is StimCopy
% 2018-09-25 AL - Now detects voltage channel by minVoltageMv and maxVoltageMv
%

%% Hard-coded parameters
validExpModes = {'EEG', 'patch', ''};
minConductancePs = 2000;    % minimum possible conductance in pS
stimCopyRange = 5;          % range of stim copy channel is 5
tolerance = 0.1;            % relative tolerance when comparing values

% Set the possible channel types, units, labels
channelTypeEEG = 'Voltage';
channelUnitsEEG = 'uV';
channelLabelEEG = 'EEG Amp (uV)';
channelTypePatchVoltage = 'Voltage';
channelUnitsPatchVoltage = 'mV';
channelLabelPatchVoltage = 'Voltage (mV)';
channelTypePatchCurrent = 'Current';
channelUnitsPatchCurrent = 'pA';
channelLabelPatchCurrent = 'Current (pA)';
channelTypePatchConductance = 'Conductance';
channelUnitsPatchConductance = 'nS';
channelLabelPatchConductance = 'Conductance (nS)';
channelTypeStimCopy = 'Voltage';
channelUnitsStimCopy = 'V';
channelLabelStimCopy = 'Stim Copy (V)';
channelTypeStimulator = 'Current';
channelUnitsStimulator = 'uA';
channelLabelStimulator = 'Stim (uA)';

%% Default values for optional arguments
expModeDefault = '';        % set later

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

% Add required inputs to the Input Parser
addRequired(iP, 'abfdata', ...
    @(x) validateattributes(x, {'numeric'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ExpMode', expModeDefault, ...
    @(x) isempty(x) || any(validatestring(x, validExpModes)));

% Read from the Input Parser
parse(iP, abfdata, varargin{:});
expMode = validatestring(iP.Results.ExpMode, validExpModes);

%% Extract from arguments
nChannels = size(abfdata, 2);       % number of channels
nDimensions = ndims(abfdata);       % number of dimensions in data

% Find data dimensions and make sure it is <= 3
if nDimensions > 3
    error('Cannot parse data with more than 3 dimensions!');
end

% Set experiment mode (if not provided) based on data dimensions
if isempty(expMode)
    switch nDimensions
    case 2
        expMode = 'EEG';
    case 3
        expMode = 'patch';
    otherwise
        error('nDimensions unrecognize!');
    end
end

%% Compute stuff
% Maximum absolute ranges for each channel
ranges = zeros(1, nChannels);

% Maximum values for each channel
maxValues = zeros(1, nChannels);
minValues = zeros(1, nChannels);

% Loop through each channel
parfor j = 1:nChannels
    % Extract all traces in this channel
    vecAll = squeeze(abfdata(:, j, :));
    
    % Compute the maximum value over all traces for this channel
    maxValues(j) = apply_iteratively(@max, vecAll);

    % Compute the minimum value over all traces for this channel
    minValues(j) = apply_iteratively(@min, vecAll);

    % Compute the maximum range over all traces for this channel
    ranges(j) = abs(maxValues(j) - minValues(j));
end

% Compute mean values for each channel
avgs = mean(mean(abfdata, 1), 3);

% Compute minimum and maximum values for each channel
minimums = min(min(abfdata, [], 1), [], 3);
maximums = max(max(abfdata, [], 1), [], 3);

%% Classify each channel
channelTypes = cell(1, nChannels);      % channel types
channelUnits = cell(1, nChannels);      % channel units
channelLabels = cell(1, nChannels);     % channel labels
if strcmpi(expMode, 'EEG')
    % EEG channels are all the same
    parfor iChannel = 1:nChannels
        channelTypes{iChannel} = channelTypeEEG;
        channelUnits{iChannel} = channelUnitsEEG;
        channelLabels{iChannel} = channelLabelEEG;
    end
elseif strcmpi(expMode, 'patch') && nChannels == 1
    % Assume that only PatchCurrent or PatchVoltage is present
    idxVoltage = identify_voltage_or_current(minimums, maximums, ranges);

    % Set channel types and units accordingly
    if idxVoltage == 1
        channelTypes{1} = channelTypePatchVoltage;
        channelUnits{1} = channelUnitsPatchVoltage;
        channelLabels{1} = channelLabelPatchVoltage;
    else
        channelTypes{1} = channelTypePatchCurrent;
        channelUnits{1} = channelUnitsPatchCurrent;
        channelLabels{1} = channelLabelPatchCurrent;
    end
elseif strcmpi(expMode, 'patch') && nChannels == 2
    % Assume that only PatchCurrent and PatchVoltage are present
    [idxVoltage, idxCurrent] = ...
        identify_voltage_or_current(minimums, maximums, ranges);

    % Set channel types and units accordingly
    channelTypes{idxVoltage} = channelTypePatchVoltage;
    channelUnits{idxVoltage} = channelUnitsPatchVoltage;
    channelLabels{idxVoltage} = channelLabelPatchVoltage;
    channelTypes{idxCurrent} = channelTypePatchCurrent;
    channelUnits{idxCurrent} = channelUnitsPatchCurrent;
    channelLabels{idxCurrent} = channelLabelPatchCurrent;
elseif strcmpi(expMode, 'patch') && nChannels == 3
    % Ranges of first two channels are similar enough, assume both are voltage
    if abs(ranges(1) - ranges(2)) < 5
        % TODO: AL - what about check the autocorrelation vector? 
        %   then you can be sure it's recording the same channel
        channelTypes = {channelTypePatchVoltage, ...
                        channelTypePatchVoltage, ...
                        channelTypePatchCurrent};
        channelUnits = {channelUnitsPatchVoltage, ...
                        channelUnitsPatchVoltage, ...
                        channelUnitsPatchCurrent};
        channelLabels = {channelLabelPatchVoltage, ...
                        channelLabelPatchVoltage, ...
                        channelLabelPatchCurrent};
    else    
        % The third channel is either StimCopy or conductance
        % Look for the range value that is closest to 5 or 0. 
        %   This will be assigned as StimCopy channel if it is within tolerance
        %   Otherwise, the third channel is assumed to be conductance
        [~, idxTemp] = min(abs(ranges - stimCopyRange));
        if abs(ranges(idxTemp) - stimCopyRange) < tolerance || ...
            abs(avgs(idxTemp)) < tolerance
            idxStimCopy = idxTemp;
            channelTypes{idxStimCopy} = channelTypeStimCopy;
            channelUnits{idxStimCopy} = channelUnitsStimCopy;
            channelLabels{idxStimCopy} = channelLabelStimCopy;
            idxThird = idxStimCopy;
        else
            % Assume greatest average is "Conductance" (generally positive)       
            %   Then since all ranges are positive, assign to -Inf to 
            %       remove conductance from subsequent detection
            [~, idxConductance] = max(avgs);
            channelTypes{idxConductance} = channelTypePatchConductance;

            % Conductance usually exceeds 2000 when picoSiemens
            if maxValues(idxConductance) > minConductancePs    
                channelUnits{idxConductance} = 'pS';
                channelLabels{idxConductance} = 'Conductance (pS)';
            else
                % Use default units for conductance
                channelUnits{idxConductance} = channelUnitsPatchConductance;
                channelLabels{idxConductance} = channelLabelPatchConductance;
            end

            idxThird = idxConductance;
        end

        % Get the remaining indices
        remaining = ones(1, 3);
        remaining(idxThird) = 0;
        indRemaining = find(remaining);

        % Assume the remaining two channels are PatchCurrent and PatchVoltage
        [idxVoltageTemp, idxCurrentTemp] = ...
            identify_voltage_or_current(minimums(indRemaining), ...
                            maximums(indRemaining), ranges(indRemaining));
        idxVoltage = indRemaining(idxVoltageTemp);
        idxCurrent = indRemaining(idxCurrentTemp);

        % Set channel types and units accordingly
        channelTypes{idxVoltage} = channelTypePatchVoltage;
        channelUnits{idxVoltage} = channelUnitsPatchVoltage;
        channelLabels{idxVoltage} = channelLabelPatchVoltage;
        channelTypes{idxCurrent} = channelTypePatchCurrent;
        channelUnits{idxCurrent} = channelUnitsPatchCurrent;
        channelLabels{idxCurrent} = channelLabelPatchCurrent;
    end
elseif strcmpi(expMode, 'patch') && nChannels == 4
    % Assume the first two channels are PatchCurrent and PatchVoltage
    [idxVoltage, idxCurrent] = ...
        identify_voltage_or_current(minimums(1:2), maximums(1:2), ranges(1:2));

    % Set channel types and units accordingly
    channelTypes{idxVoltage} = channelTypePatchVoltage;
    channelUnits{idxVoltage} = channelUnitsPatchVoltage;
    channelLabels{idxVoltage} = channelLabelPatchVoltage;
    channelTypes{idxCurrent} = channelTypePatchCurrent;
    channelUnits{idxCurrent} = channelUnitsPatchCurrent;
    channelLabels{idxCurrent} = channelLabelPatchCurrent;

    % Remove from subsequent detection
    ranges(1:2) = -Inf;

    % Look for the range value that is closest to 5 or 0. 
    %   This will be assigned as StimCopy channel if it is within tolerance
    [~, idxTemp] = min(abs(ranges - stimCopyRange));
    if abs(ranges(idxTemp) - stimCopyRange) < tolerance || ...
        abs(avgs(idxTemp)) < tolerance
        idxStimCopy = idxTemp;
    else
        idxStimCopy = 3;
    end

    % The other index is Stimulator
    idxStimulator = 2 + (3 - (idxStimCopy - 2));

    % Set channel types and units accordingly
    channelTypes{idxStimCopy} = channelTypeStimCopy;
    channelUnits{idxStimCopy} = channelUnitsStimCopy;
    channelLabels{idxStimCopy} = channelLabelStimCopy;
    channelTypes{idxStimulator} = channelTypeStimulator;
    channelUnits{idxStimulator} = channelUnitsStimulator;
    channelLabels{idxStimulator} = channelLabelStimulator;
elseif strcmpi(expMode, 'patch')
    error('Too many channels!');
else
    error('Channels cannot be identified!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idxVoltage, idxCurrent] = ...
            identify_voltage_or_current (minimums, maximums, ranges)
%% Identify which of two channels is voltage/current

% Hard-coded parameters
minVoltageMv = -120;        % minimum possible voltage in mV
maxVoltageMv = 50;          % maximum possible voltage in mV

% Make sure there are only two indices
if numel(minimums) > 2 || numel(maximums) > 2 || numel(ranges) > 2
    idxVoltage = NaN;
    idxCurrent = NaN;
    return;
end

% Determine whether each channel is within the 'Voltage range'
inVoltageRange = minimums >= minVoltageMv & maximums <= maxVoltageMv;

% If both or none are within range, use the range of values to differentiate
%   Otherwise, 
if all(inVoltageRange) || ~any(inVoltageRange)
    % The index with the maximum range is PatchCurrent
    [~, idxCurrent] = max(ranges);

    % The other index is PatchVoltage
    idxVoltage = 3 - idxCurrent;
else
    % The index within the 'Voltage range' is PatchVoltage
    idxVoltage = find(inVoltageRange);

    % The other index is PatchCurrent
    idxCurrent = 3 - idxVoltage;        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% AL - This is a bit confusing
%       vig     - numeric array of size 2 or 3, 
%                   with the following possible elements:
%                   1 - Voltage
%                   2 - Current
%                   3 - Conductance
%               e.g. [2 1 3] -> 1st channel: Current, 2nd channel: Voltage, 3rd channel: Conductance
%
vig = zeros(1, nChannels);      % array marker for voltage, current, conductance

% recording mode is [Voltage Voltage Current]

allLabels = {'Voltage (mV)', 'Current (pA)', 'Conductance (nS)'};
channelLabels{idxVoltage} = allLabels{1};
channelLabels{idxCurrent} = allLabels{2};
channelLabels = {allLabels{1}, allLabels{1}, allLabels{2}};
channelLabels{idxVoltage} = allLabels{1};
channelLabels{idxCurrent} = allLabels{2};
channelLabels{idxConductance} = 'Conductance (pS)';
channelLabels{idxConductance} = allLabels{3};

allTypes = {'Voltage', 'Current', 'Conductance'};
allUnits = {'mV', 'pA', 'nS'};

% Remove current from detection, voltage is the other index
ranges(idxCurrent) = -Inf;
[~, idxVoltage] = max(ranges);

% Construct the channel labels from channelTypes and channelUnits
channelLabels = cellfun(@(x, y) [x , ' (', y, ')'], ...
                        channelTypes, channelUnits, ...
                        'UniformOutput', false);

% Assume next highest range is "Current"
[~, idxCurrent] = max(ranges);
ranges(idxCurrent) = -Inf;

% Remove current from ranges, leaving remaining channel as "Voltage"
[~, idxVoltage] = max(ranges);

% The index with the maximum range is PatchCurrent
[~, idxCurrent] = max(ranges(1:2));

% The other index is PatchVoltage
idxVoltage = 3 - idxCurrent;

ranges(idxConductance) = -Inf;
ranges(idxStimCopy) = -Inf;

maxValues(j) = max(max(vecAll));
ranges(j) = abs(max(max(vecAll)) - min(min(vecAll)));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%