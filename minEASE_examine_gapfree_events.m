function [toQuit, toRedo, eventInfoOut, eventClassOut, isCheckedOut] = ...
                minEASE_examine_gapfree_events (eventInfoLast, eventClassLast, isCheckedLast, dataRaw, dataLowpass, tVec, outputDirectory, outputLabel, directionPsc, varargin)
%% Examine detected events and allow additions and deletions
% Usage: [toQuit, toRedo, eventInfoOut, eventClassOut, isCheckedOut] = ...
%               minEASE_examine_gapfree_events (eventInfoLast, eventClassLast, isCheckedLast, dataRaw, dataLowpass, tVec, outputDirectory, outputLabel, directionPsc, varargin)
% Explanation: 
%       TODO
% Outputs:
%       TODO
% Side Effects:
%       TODO
% Arguments:    
%       TODO
%       varargin    - 'NoiseLevel': Noise level for the signal to compare with
%                   must be a positive scalar
%                   default == root-mean-square level of Gaussian part of data
%                   - 'ZoomWindowMs': window for examining events (ms)
%                   must be a numeric positive scalar
%                   default == 50 ms
%                   - 'ToPrompt': whether to prompt for decisions
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   
%
% Requires:
%       cd/minEASE_gui_examine_events.m
%       cd/compute_rms_Gaussian.m
%
% Used by:
%       cd/minEASE.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-06-06 AL - Renamed goldfinger_minis.m -> minEASE_examine_gapfree_events.m
% 2017-06-12 AL - Passed directionPsc to GUI for adjust_peaks
% 2017-07-24 AL - Moved code to dlmwrite_with_header.m
% 2017-07-24 AL - Instead of using Auto versions, now uses Last versions
%                   of eventInfo, eventClass, isChecked
% 2017-07-25 AL - Now saves sweep labels too
% 2017-10-15 AL - Moved saving events to minEASE.m
% 2018-01-28 AL - Added isdeployed
% 2018-02-13 AL - Added NoPrompts
% 2018-02-16 AL - Now passes noiseLevel to minEASE_gui_examine_events.m
% 2018-02-23 AL - Now checks whether GUI is closed every 0.1 seconds
% 2018-02-27 AL - Changed NoPrompts to toPrompt
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
%

% Shared with minEASE_gui_examine_events.m
global eventInfo eventClass isChecked                   % used
global finalCommand                                     % used

%% Hard-coded parameters
PAUSE_TIME = 0.1;

%% Default parameters for TODO
zoomWindowMsDefault = 50;           % default window for examining events (ms)
toPromptDefault = true;             % whether to prompt for decisions by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 7
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NoiseLevel', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ZoomWindowMs', zoomWindowMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ToPrompt', toPromptDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
noiseLevel   = iP.Results.NoiseLevel;
zoomWindowMs = iP.Results.ZoomWindowMs;
toPrompt     = iP.Results.ToPrompt;

% If no noiseLevel provided, use the root mean square of the 
%   Gaussian part of the original data trace as the default
%   Note: this should not happen in minEASE
if isempty(noiseLevel)
    noiseLevel = compute_rms_Gaussian(data);            % default Gaussian noise level
end

%% Initialize outputs
toQuit = false;
toRedo = true;
eventInfo = eventInfoLast;
eventClass = eventClassLast;
isChecked = isCheckedLast;

%% Create and show GUI
%% TODO: Pass in other optional arguments
[hfig.gui] = minEASE_gui_examine_events(dataRaw, dataLowpass, tVec, ...
                                outputLabel, directionPsc, ...
                                'NoCloseReq', ~toPrompt, ...
                                'NoiseLevel', noiseLevel);

% Wait for user to close GUI before continuing
while ishandle(hfig.gui)    % while GUI is still there
    % Pause for PAUSE_TIME seconds
    pause(PAUSE_TIME);          
end

% Upon clicking "DONE!" (closes GUI), radiobutton tag is returned
switch finalCommand
case 'NextTrace'
    toRedo = false;
case 'StartOver'
    %toRedo = true;
case 'Quit'
    toQuit = true;
end

%% TODO: Don't close GUI but replot instead
% Make sure all other figures are closed
close all force

% Return new eventInfo, eventClass & isChecked
eventInfoOut = eventInfo;
eventClassOut = eventClass;
isCheckedOut = isChecked;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% ccc
% load('/home/mark/matlab_temp_variables/goldfinger_minisVars')
end

skip = 0;
if skip ==0
end

% plot(tVec, dataRaw, 'color', rgb('silver'))
% hold on
% 
% hpanel=uipanel('Position',[0 .05 2 .95]);
% hscrollbar=uicontrol('style', 'slider', 'Units', 'normalized', 'Position',[0 0 1 .05], 'callback',@hscroll_CallbackMini);
% axes('parent',hpanel, 'outerPosition',[.25 0 .5 1])

%     timeIDX = cursorData(i, 1);
%     correctCursorData(i,3) = maxV;
%     correctCursorData(i,4) = maxVidx;  
%     correctCursorData(i,3) = tVec(maxVidx);
%     correctCursorData(i,4) = maxV;   

%     plot(cursorData(:, 1), cursorData(:, 2), 'ro');    
 % set(gcf, 'Units', 'normalized', 'Position', [0 0 0.9 0.7])
% set(gcf, 'Units', 'normalized', 'Position', [0 0.1 0.23 0.3])
% make_my_figure_fit_horizontal(6, 16)

%% print fig
% figPN = strcat('/home/barrettlab/public_html/detectorama/', outputDirectory, '/');
% mkdir(figPN)
% figFN = sprintf('Cell_%03d',cellNum);
% figFNfull = strcat(figPN, figFN);
% print(figFNfull, '-dpdf', '-r200')

%     load('/home/mark/matlab_temp_variables/detectedMinis.mat');    
%     baseVals = [];

uiRangeSamples = find(tVec >= uiRangeMs + tVec(1), 1);

function [runThru, toRedo, uiRangeMs, quitButt] = goldfinger_minis(dataRaw, dataFilt, tVec, cellNum, outputDirectory, runThru, toRedo, uiRangeMs, directionPsc, baseVals)

if rad1 == 1 && rad2 == 0 && rad3 == 0 && rad4 == 0
if rad1 == 0 && rad2 == 1 && rad3 == 0 && rad4 == 0
if rad1 == 0 && rad2 == 0 && rad3 == 1 && rad4 == 0
if rad1 == 0 && rad2 == 0 && rad3 == 0 && rad4 == 1

    rad1 = get(r1, 'Value');
    rad2 = get(r2, 'Value');
    rad3 = get(r3, 'Value');
    rad4 = get(r4, 'Value');

if rad1
if rad2
if rad3
if rad4

posFig1 = [0, 0.2, 0.5, 0.4];       % position of Figure 1 in normalized units

moddedEvents = [];
save('/home/mark/matlab_temp_variables/goldfinger_minisVars')
load('/home/mark/matlab_temp_variables/detectedMinis.mat');
load('/home/mark/matlab_temp_variables/baselineMinis.mat')
save('/home/mark/matlab_temp_variables/moddedMinis', 'moddedEvents');

% Correct the cursor dataRaw
for i = 1:size(cursorData, 1)
    if i > 25
        x = 1;              % TODO: What's this?
    end
    timeIDX = find(tVec >= cursorData(i, 1), 1);
    chunk = timeIDX - uiRangeSamples:timeIDX + uiRangeSamples;
    switch directionPsc
    case {'E', 'Excitatory'}
        maxV = min(dataRaw(chunk));
    case {'I', 'Inhibitory'}
        maxV = max(dataRaw(chunk));
    end
    maxVidx = find(dataRaw(chunk) == maxV, 1) + timeIDX-uiRangeSamples-1;
    correctCursorData(i, 1) = tVec(maxVidx);
    correctCursorData(i, 2) = maxV; 
    if i == size(cursorData, 1)     % TODO: What for?
        clear cursorData
        cursorData = correctCursorData;
    end
end

if ~isempty(correctCursorData)
    save('/home/mark/matlab_temp_variables/detectedMinisCorrect.mat', 'correctCursorData');
    pause(1)
end

if ~isempty(cursorData)
    plot(baseVals(:, 1), baseVals(:, 2), 'x', ...
        'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
    plot(cursorData(:, 1), cursorData(:, 2), 'o', ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 4);
end

runThru = 0;
if ~runThru
end

uiRange = str2num(get(winSize, 'String'))/1000;

    cursorData = [];
    cursorData = [];
    cursorData = [];
    cursorData = [];

siMs = tVec(2) - tVec(1);
uiRangeSamples = round(uiRangeMs/siMs);

quitButt = get(btnQuit, 'Value');
uiRangeMs = str2num(get(handles.winSize, 'String'));

case 'Add Events'
    load('/home/mark/matlab_temp_variables/moddedMinis');  
    [cursorDataMod, baseValsMods] = ...
        addEvents(dataRaw, tVec, correctCursorData, ...
                    moddedEvents, baseVals, directionPsc, uiRangeMs);
    cursorData = cat(1,correctCursorData, cursorDataMod);
    baseVals = cat(1, baseVals, baseValsMods);
    save('/home/mark/matlab_temp_variables/detectedMinis', 'cursorData')
    save('/home/mark/matlab_temp_variables/baselineMinis', 'baseVals')
case 'Remove Events'
    load('/home/mark/matlab_temp_variables/detectedMinis.mat');
    load('/home/mark/matlab_temp_variables/moddedMinis');    
    [cursorDataMod, baseValsMods] = deleteEvents(cursorData, moddedEvents, baseVals);
    cursorData = cursorDataMod;
    baseVals = [];
    baseVals = baseValsMods;
    save('/home/mark/matlab_temp_variables/detectedMinis', 'cursorData')
    save('/home/mark/matlab_temp_variables/baselineMinis', 'baseVals')

    load('/home/mark/matlab_temp_variables/autoMinis');
    save('/home/mark/matlab_temp_variables/detectedMinis', 'cursorData')
    pause(1)

    if ~exist('xPOINTS', 'var') && isempty(cursorData)
        correctCursorData = [];
        pause(1)
    end
    
    
switch handles.ButtonGroup.SelectedObject.Tag

% Wait for user response before continuing
uiwait(gcf);

% Upon clicking "DONE!" (calling uiresume(gcbf)), get button values

    fid = fopen(csvFileNameWHeader, 'w');
    for j = 1:numel(outputCellHeader)
        fprintf(fid, '%s, ', outputCellHeader{j});
    end
    fprintf(fid, '\n');
    if ~isempty(outputMatrix)
        if numel(outputCellHeader) ~= size(outputMatrix, 2)
            error('Header is no the correct length!');
        end
        for i = 1:size(outputMatrix, 1)
            for j = 1:size(outputMatrix, 2)
                fprintf(fid, '%g, ', outputMatrix(i, j));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);

                    % Save siMs & nSamples
                    siMs = tVec(2) - tVec(1);       % sampling interval (ms)
                    nSamples = length(tVec);        % number of samples

%                   - 'NoPrompts': whether to suppress prompts
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
addParameter(iP, 'NoPrompts', noPromptsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
noPrompts    = iP.Results.NoPrompts;

%}   
