function hGUI = minEASE_gui_examine_events (dataRaw, dataLowpass, tVec, ...
                                    outputLabel, directionPsc, varargin)
%% A graphic user interface for examining detected events
% Usage: hGUI = minEASE_gui_examine_events (dataRaw, dataLowpass, tVec, ...
%                                   outputLabel, directionPsc, varargin)
% Explanation: 
%       TODO
% Outputs:
%       TODO
% Side Effects:
%       TODO
% Arguments:    
%       TODO
%       varargin    - 'ZoomWindowInit': window for examining events (ms)
%                   must be a numeric positive scalar
%                   default == 50 ms
%                   - 'CorrRangeInit': range for cursor correction (ms)
%                   must be a numeric positive scalar
%                   default == 1 ms
%                   - 'NoCloseReq': whether to suppress close request function
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
% Requires:
%       cd/minEASE_plot_gapfree_events.m
%       /home/Matlab/Downloaded_Functions/rgb.m
%       cd/compute_rms_Gaussian.m
%       cd/validate_string.m
%       cd/increment_editbox.m
%       cd/find_custom.m
%       cd/adjust_peaks.m
%       cd/my_closereq.m
%       cd/get_idxEnd.m
%
% Used by:
%       cd/minEASE_examine_gapfree_events.m

% File History:
% ---------- Created by Koji Takahashi & Mark P Beenhakker
% 2017-06-06 AL - Renamed myuiDetectMinis.m -> minEASE_gui_examine_events.m
% 2017-06-06 AL - Created full GUI with the figure not overlapping buttons
% 2017-06-07 AL - Added class number and event number selection
% 2017-06-07 AL - Removed eventInfo & eventClass from arguments and make 
%                   them global variables
% 2017-06-08 AL - Create a checkbox for the zoomed event's isChecked status;
%                   gray out when no events zoomed
% 2017-06-09 AL - Added hPlot to manage groups of plots
% 2017-06-09 AL - Display filled circles for checked events and 
%                   allow to be removed if unchecked
% 2017-06-09 AL - Removed lastClassInfo, nLastClass from global variables
% 2017-06-09 AL - Updated extract_class arguments and outputs
% 2017-06-09 AL - Made classNoStr, eventNoStr, zoomWindowMs global variables
% 2017-06-09 AL - Replaced classNoStrLast & eventNoLast by eventNoInAllLast
% 2017-06-09 AL - Now deals with out of range issues in callback functions
% 2017-06-09 AL - Made the buttonUps and buttonDowns go circular
% 2017-06-09 AL - Made commonly shared gui handles global variables for speed
% 2017-06-09 AL - Moved code to increment_editbox.m
% 2017-06-10 AL - Fixed empty string: '' is 0x0 whereas an empty String field
%                   is 1x0, so strcmpi() will give 0. Use isempty() instead.
% 2017-06-10 AL - Added the ability to change class
% 2017-06-11 AL - Added the ability to remove events
% 2017-06-12 AL - Added MainGUI_WindowScrollWheelFcn to increment zoom window
%                   upon mouse scroll
% 2017-06-12 AL - Added the ability to add events
% 2017-06-12 AL - Simplified initialization
% 2017-06-12 AL - Now checks if examined event is changed before calling 
%                   zoom_event(), update_checkbox_status() & annotate_event()
% 2017-06-13 AL - Fixed nuances for adding events
% 2017-06-14 AL - Add event is now interruptible (abortable)
% 2017-06-14 AL - Added 'all' pushbuttons
% 2017-06-14 AL - Added pushbuttons that will increment/decrement events and 
%                   mark checked at the same time (<< and >>)
% 2017-06-14 AL - Added pageup/pagedown to increment/decrement events and 
%                   mark checked at the same time
% 2017-06-14 AL - Added pushbuttons that will increment/decrement zoom windows 
%                   in increments of 1 ms & 100 ms
% 2017-06-14 AL - Added pageup/pagedown to increment/decrement zoom windows 
%                   in increments of 5 ms & 50 ms
% 2017-06-14 AL - Added window scroll function to increment/decrement 
%                   zoom windows in 10 ms
% 2017-06-14 AL - Added reset button for zoom window
% 2017-07-31 AL - Now skips to the next unchecked event when using keyboard
%                   to increment/decrement event number
% 2017-08-03 AL - Made font sizes a little bit smaller for compatibility
% 2017-08-07 AL - Fixed calculation of rise10Vals & rise90Vals
% 2017-08-07 AL - Now recomputes IEIs, ISIs and decay times after 
%                   adding/changing/removing an event
% 2017-08-07 AL - Removed events now have NaNs for IEIs, ISIs and decay times 
%                   (for consistency)
% 2017-08-07 AL - Now skips events of class 8 when checking previous
%                   or next event
% 2017-08-30 AL - Fixed ButtonRemove_Callback by making nSamples global
% 2017-10-16 AL - Change Slow Decay -> Wrong Decay
% 2017-10-17 AL - Changed Button Group to Button Group 2 and 
%                   added Button Group 1 to rank events by time, 
%                   by amplitude or by decay time
% 2017-10-17 AL - Added eventRank & rankThisClass and 
%                   now increments events according to rankThisClass
% 2017-10-18 AL - Added 'Rise' to Button Group 1
% 2017-10-19 AL - Moved verify_classNoNew() into change_class()
%                   so that it is executed after computing 
%                   next-event-dependent statistics in case an event is added
% 2018-01-28 AL - Added isdeployed
% 2018-01-29 AL - Removed verify_classNoNew() and added classify_PSC()
% 2018-01-29 AL - Now reclassifies previous and next PSCs 
%                   when events are removed
% 2018-01-30 AL - Added checkboxes for "All of class" and "Don't prompt"
% 2018-01-30 AL - Added doForAll and ToPrompt
% 2018-01-30 AL - Now suppresses prompt whenever changing class in bulk
% 2018-01-30 AL - Implemented bulk remove and bulk change
% 2018-01-30 AL - Created update_peak_markers() and allowed update
%                   to happen after bulk change
% 2018-01-30 AL - Replaced all waitfor with uiwait
% 2018-02-02 AL - Changed minLine and maxLine to be the minimum 
%                   and maximum of the y-axis at the 'all' plot
% 2018-02-02 AL - Autoscaled the y-axis manually because of the gray line
% 2018-02-03 AL - y-axis limits are now 130% of full data range
% 2018-02-03 AL - Changed ButtonUp2Zoom and ButtonDown2Zoom to use x2 and x0.5
%                   increments, respectively
% 2018-02-13 AL - Added NoCloseReq
% 2018-02-23 AL - Added Clean Mode
% 2018-02-23 AL - Now allows user to change the cursor correction range 
%                   while adding events or under clean mode
% 2018-02-23 AL - Added 'HitTest', 'off' to all plots
% 2018-02-27 AL - Now gives focus to editEventNo every time the class number 
%                   is changed. The GUI is also initialized with this focus
% 2018-02-27 AL - Now gives focus to editEventNo after zooming with mouse scroll
% 2018-02-27 AL - Now does not jump to a new event when removing 
%                   an event not zoomed
% 2018-02-27 AL - Fixed bug on incrementing event numbers 
%                   in a class with no events
% 2018-02-28 AL - Now does not change displayed class number for all cases
% 2018-02-28 AL - The 'insert' key is now the same as the Change button
%                   when cursor is in editEventNo
% 2018-02-28 AL - All mouse clicks and key presses will now give focus 
%                   to editEventNo at the very end, except when incrementing
%                   or decrementing an edit field
% 2018-02-28 AL - Fixed bulk changing so that it will still work 
%                   when Class is at 'All'
% 2018-03-01 AL - Fixed bug in Add Events so that new event will be in focus
%                   after being added
% 2018-03-01 AL - Now allows the user to exit the program during event addition
% 2018-03-01 AL - Now automatically turns buttonClean off 
%                   when buttonAdd is clicked
% 2018-03-01 AL - Now does not allow user to add an event 
%                   that overlaps with a previous event
% 2018-03-01 AL - Fixed bug after event addition is aborted
% 2018-03-15 AL - Changed corrRangeInitDefault from 1 -> 0.5 ms 
%                   per Paula's request
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
% 2018-08-03 AL - Updated legend to turn 'AutoUpdate' off for R2017a and beyond
% TODO: 2018-08-03 AL - Added button to remove events completely
% TODO: Add scroll bar?
% TODO: Update previous and next events that are not PSCs when 
%       events are removed
% TODO: Save eventInfo eventClass isChecked in a matfile whenever they are
%           changed?
%

global posTextBox fontSizeTextBox fontWeightTextBox     % updated
global textColorTextBox bgColorTextBox                  % updated
global classNumbers classColorMap classLabels           % updated
global markerSize markerEdgeWidth                       % updated
global colorButton colorButtonOn                        % updated
global nEvents nClass nSamples siMs                     % updated
global hPlot                                            % updated
global zoomWindowMs corrRangeMs classNoNewUser          % updated
global zoomWindowInit                                   % updated
global noiseLevel                                       % updated
global eventNoInAllLast                                 % updated
global doForAll toPrompt byForce                        % updated
global gui                                              % updated
global editClassNo editEventNo editZoomWindow           % updated 
global dispNThisClass dispNEvents                       % updated
global checkboxChecked checkboxNoPrompt                 % updated
global editCorrRange editAddClass                       % updated
global buttonClean                                      % updated
global buttonGroup1                                     % updated
global radioButtons2 buttonGroup2 buttonDone            % updated
global panelLeft allUiControl                           % updated
global panelRight sweepFigure                           % updated
global swpT swpD eOrIFactor                             % updated
global allowAddClick allowRemoveClick                   % updated
global minD maxD                                        % updated
global minY maxY                                        % updated
global abort                                            % updated
global isIncrDecr                                       % updated
global eventInfo eventClass isChecked                   % used

% Plot parameters
markerSize = 10;                        % marker size
markerEdgeWidth = 2;                    % marker edge width
posTextBox = [0.14, 0.76, 0.2, 0.2];    % position of annotation text box
fontSizeTextBox = 12;                   % font size of annotation text box
fontWeightTextBox = 'bold';             % font weight of annotation text box
textColorTextBox = 'Indigo';            % text color of annotation text box
bgColorTextBox = 'Bisque';              % background color of text box

% Class numbers, colors and labels
classNumbers = 1:8;                     % class numbers go from 1 through 8
classColors = {'Lime', 'Red', 'Gold', ...
              'Cyan', 'DarkOrange', 'Purple', ...
              'Violet', 'Gray'};
classLabels = {'Class 1: Type I PSC', 'Class 2: Type II PSC', ...
               'Class 3: Type III PSC', 'Class 4: Wrong Decay', ...
               'Class 5: Slow Rise', 'Class 6: Too Small', ...
               'Class 7: In Seal Test', 'Class 8: Deleted'};

% Colors
colorPanelLeft = 'Bisque';      % background color of Left Panel
colorPanelRight = 'Cornsilk';   % background color of Right Panel
colorSlider = 'Khaki';          % background color of sliders
colorButton = 'Lime';           % color of pushbuttons & togglebuttons
colorButtonOn = 'Red';          % color of togglebuttons when on
colorText = 'Indigo';           % foreground color for texts
colorTextButton = 'DodgerBlue'; % foreground color for big buttons

% Font names
textFontNameDefault = 'Arial';
textFontNameButton = 'FixedWidth';

% Font sizes units
textFontUnitsDefault = 'normalized';
textFontUnitsButton = 'normalized';

% Font sizes in normalized units
textFontSizeDefault = 0.5;
textFontSizeButton = 0.4;
radiobuttonFontSize = 0.7;
checkboxFontSize = 0.6;

% Font weights
textFontWeightDefault = 'bold';
textFontWeightButton = 'bold';

% Outer position of MainGUI figure in normalized units, 
%       [left bottom width height]
outPosGui = [0.1, 0.1, 0.8, 0.8];

% Position of panels within MainGUI figure in normalized units
posPanelLeft  = [0.0, 0, 0.2, 1];
posPanelRight = [0.2, 0, 0.8, 1];

% Position of objects within PanelLeft in normalized units
posTextNEvents     = [0.00, 0.950, 0.80, 0.04];
posDispNEvents     = [0.80, 0.950, 0.20, 0.04];
posTextClassNo     = [0.00, 0.910, 1.00, 0.04];
posButtonAllClass  = [0.40, 0.900, 0.20, 0.02];
posEditClassNo     = [0.40, 0.860, 0.20, 0.04];
posButtonDownClass = [0.30, 0.860, 0.10, 0.04];
posButtonUpClass   = [0.60, 0.860, 0.10, 0.04];
posTextNThisClass  = [0.00, 0.800, 0.80, 0.04];
posDispNThisClass  = [0.80, 0.800, 0.20, 0.04];
posTextEventNo     = [0.00, 0.760, 1.00, 0.04];
posButtonAllNo     = [0.40, 0.750, 0.20, 0.02];
posEditEventNo     = [0.40, 0.710, 0.20, 0.04];
posButtonDown2No   = [0.20, 0.710, 0.10, 0.04];
posButtonDownNo    = [0.30, 0.710, 0.10, 0.04];
posButtonUpNo      = [0.60, 0.710, 0.10, 0.04];
posButtonUp2No     = [0.70, 0.710, 0.10, 0.04];
posButtonGroup1    = [0.00, 0.660, 1.00, 0.033];
posTextZoomWindow  = [0.00, 0.610, 1.00, 0.04];
posButtonResetZoom = [0.40, 0.600, 0.20, 0.02];
posEditZoomWindow  = [0.40, 0.560, 0.20, 0.04];
posButtonDown2Zoom = [0.15, 0.560, 0.15, 0.04];
posButtonDownZoom  = [0.30, 0.560, 0.10, 0.04];
posButtonUpZoom    = [0.60, 0.560, 0.10, 0.04];
posButtonUp2Zoom   = [0.70, 0.560, 0.15, 0.04];
posCheckboxChecked = [0.15, 0.510, 0.70, 0.04];
posButtonClean     = [0.05, 0.435, 0.40, 0.07];
posButtonAdd       = [0.55, 0.435, 0.40, 0.07];
posTextCorrRange   = [0.00, 0.390, 0.80, 0.04];
posEditCorrRange   = [0.80, 0.395, 0.15, 0.04];
posTextAddClass    = [0.00, 0.350, 0.80, 0.04];
posEditAddClass    = [0.80, 0.355, 0.15, 0.04];
posButtonRemove    = [0.05, 0.280, 0.40, 0.07];
posButtonChange    = [0.55, 0.280, 0.40, 0.07];
posCheckboxForAll  = [0.00, 0.240, 0.50, 0.04];
posCheckboxNoPrompt= [0.50, 0.240, 0.50, 0.04];
posCheckboxByForce = [0.00, 0.200, 0.50, 0.04];
posCheckboxReserved= [0.50, 0.200, 0.50, 0.04];
posButtonGroup2    = [0.00, 0.100, 1.00, 0.10];
posButtonDone      = [0.05, 0.015, 0.90, 0.07];

%{
% TODO: Position of objects within PanelRight in normalized units
posButtonRightPanel      = [0.01, 0.95, 0.15, 0.04];
%}

% Position of buttons within ButtonGroup1 in normalized units
posRadioButtons1 = {[0.00, 0.0, 0.25, 1], ...
                   [0.25, 0.0, 0.25, 1], ...
                   [0.50, 0.0, 0.25, 1], ...
                   [0.75, 0.0, 0.25, 1]};

% Position of buttons within ButtonGroup2 in normalized units
posRadioButtons2 = {[0.1, 0.66, 0.8, 0.33], ...
                   [0.1, 0.33, 0.8, 0.33], ...
                   [0.1, 0.00, 0.8, 0.33]};

%% Default values for optional parameters
eventNoInitDefault = 'all';         % default event number to examine
classNoInitDefault = 'all';         % default event class number to examine
zoomWindowInitDefault = 50;         % default window for examining events (ms)
corrRangeInitDefault = 0.5;         % default range for cursor correction (ms)
classNoNewInitDefault = 1;          % default event class number to add/change
noCloseReqDefault = false;          % whether to suppres close request function
                                    %   by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NoiseLevel', [], ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'EventNoInit', eventNoInitDefault, ...
    @(x) assert((ischar(x) || isstring(x)) && strcmpi(x, 'all') || ...
                isnumeric(x) && isscalar(x) && mod(x, 1) == 0 && x > 0, ...
                ['EventNoInit must be either ''all'' ', ...
                 'or a positive integer scalar!']));
addParameter(iP, 'ClassNoInit', classNoInitDefault, ...
    @(x) assert((ischar(x) || isstring(x)) && strcmpi(x, 'all') || ...
                isnumeric(x) && isscalar(x) && mod(x, 1) == 0 && x > 0, ...
                ['ClassNoInit must be either ''all'' ', ...
                 'or a positive integer scalar!']));
addParameter(iP, 'ZoomWindowInit', zoomWindowInitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'CorrRangeInit', corrRangeInitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ClassNoNewInit', classNoNewInitDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'NoCloseReq', noCloseReqDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
noiseLevel = iP.Results.NoiseLevel;
eventNoInit = iP.Results.EventNoInit;
classNoInit = iP.Results.ClassNoInit;
zoomWindowInit = iP.Results.ZoomWindowInit;
corrRangeInit = iP.Results.CorrRangeInit;
classNoNewUser = iP.Results.ClassNoNewInit;
noCloseReq = iP.Results.NoCloseReq;

% Extract from arguments
swpD = dataLowpass;                         % current sweep data (pA)
swpT = tVec;                                % current sweep time
minD = min(swpD);                           % maximum sweep data value (pA)
maxD = max(swpD);                           % minimum sweep data value (pA)

% Compute from arguments
nEvents = size(eventInfo, 1);               % total number of events detected
nClass = numel(classNumbers);               % total number of classes
nSamples = length(swpT);                    % total number of samples
siMs = swpT(2) - swpT(1);                   % sampling interval in ms

% If no noiseLevel provided, use the root mean square of the 
%   Gaussian part of the original data trace as the default
if isempty(noiseLevel)
    noiseLevel = compute_rms_Gaussian(swpD);        % default Gaussian noise level
end

% Construct color map
classColorMap = zeros(nClass, 3);
for iClass = 1:nClass
    classColorMap(iClass, :) = rgb(classColors{iClass});
end

% Determine string for class number
if isnumeric(classNoInit)
    classNoStr = num2str(classNoInit);
else
    classNoStr = classNoInit;
end

% Set the excitatory or inhibitory factor
switch directionPsc
case {'Excitatory', 'E'}
    eOrIFactor = -1;
case {'Inhibitory', 'I'}
    eOrIFactor = 1;
otherwise
    error('This PSC direction is not supported!\n');
end

% Determine string for event number
if isnumeric(eventNoInit)
    eventNoStr = num2str(eventNoInit);
else
    eventNoStr = eventNoInit;
end

% Initialize windows
zoomWindowMs = zoomWindowInit;
corrRangeMs = corrRangeInit;

%% Change default uicontrol properties
%   Note: restore in minEASE_examine_gapfree_events.m after closing GUI
set(0, 'DefaultUicontrolFontUnits', textFontUnitsDefault, ...
       'DefaultUicontrolFontSize', textFontSizeDefault, ...
       'DefaultUicontrolFontName', textFontNameDefault, ...
       'DefaultUicontrolFontWeight', textFontWeightDefault, ...
       'DefaultUicontrolUnits', 'normalized', ...
       'DefaultFigureUnits', 'normalized');

%% Create main GUI figure
gui = figure('Visible', 'off', 'Tag', 'MainGUI', ...
             'Name', sprintf('minEASE (c) 2017 for %s', outputLabel), ...
             'OuterPosition', outPosGui, ...
             'WindowScrollWheelFcn', @MainGUI_WindowScrollWheelFcn);

% Set a close request function if not noPrompts
if ~noCloseReq
    set(gui, 'CloseRequestFcn', {@my_closereq, 'Yes', 'minEASE'});
end

%% Create left panel
panelLeft = ...
    uipanel(gui, 'Tag', 'PanelLeft', 'Position', posPanelLeft, ...
                 'Title', 'Control Panel', ...
                 'BackgroundColor', rgb(colorPanelLeft));

% Display the number of events detected
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextNEvents', ...
                     'String', 'Total number of events:', ...
                     'Position', posTextNEvents);

dispNEvents = ...
    uicontrol(panelLeft, 'Style', 'text', 'Tag', 'DispNEvents', ...
                         'String', nEvents, ...
                         'Position', posDispNEvents);

% Create an editable text box for choosing event class to examine
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextClassNo', ...
                     'String', 'Zoom in on event class:', ...
                     'Position', posTextClassNo);

editClassNo = ...
    uicontrol(panelLeft, 'Style', 'edit', 'Tag', 'EditClassNo', ...
                         'String', classNoStr, ...
                         'Position', posEditClassNo, ...
                         'KeyPressFcn', @EditClassNo_KeyPressFcn, ...
                         'Callback', @EditClassNo_Callback);

% Create pushbuttons for incrementing/decrementing event class or setting 'all'
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonAllClass', ...
                     'String', '---------------', ...
                     'Position', posButtonAllClass, ...
                     'Callback', @ButtonAllClass_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDownClass', ...
                     'String', '<', ...
                     'Position', posButtonDownClass, ...
                     'Callback', @ButtonDownClass_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUpClass', ...
                     'String', '>', ...
                     'Position', posButtonUpClass, ...
                     'Callback', @ButtonUpClass_Callback);

% Display the number of events of this class detected
%   Note: This is updated after calling EditClassNo_Callback
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextNThisClass', ...
                     'String', 'Number in this class:', ...
                     'Position', posTextNThisClass);

dispNThisClass = ...
    uicontrol(panelLeft, 'Style', 'text', 'Tag', 'DispNThisClass', ...
                         'String', '', ...
                         'Position', posDispNThisClass);

% Create an editable text box for choosing event number to examine
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextEventNo', ...
              'String', 'Event number in this class:', ...
              'Position', posTextEventNo);

editEventNo = ...
    uicontrol(panelLeft, 'Style', 'edit', 'Tag', 'EditEventNo', ...
                         'String', eventNoStr, ...
                         'Position', posEditEventNo, ...
                         'KeyPressFcn', @EditEventNo_KeyPressFcn, ...
                         'Callback', @EditEventNo_Callback);

% Create pushbuttons for incrementing/decrementing event number or setting 'all'
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonAllNo', ...
                     'String', '---------------', ...
                     'Position', posButtonAllNo, ...
                     'Callback', @ButtonAllNo_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDownNo', ...
                     'String', '<', ...
                     'Position', posButtonDownNo, ...
                     'Callback', @ButtonDownNo_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUpNo', ...
                     'String', '>', ...
                     'Position', posButtonUpNo, ...
                     'Callback', @ButtonUpNo_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDown2No', ...
                     'String', '<<', ...
                     'Position', posButtonDown2No, ...
                     'Callback', @ButtonDown2No_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUp2No', ...
                     'String', '>>', ...
                     'Position', posButtonUp2No, ...
                     'Callback', @ButtonUp2No_Callback);

% Create button group 1
buttonGroup1 = ...
    uibuttongroup(panelLeft, 'Tag', 'ButtonGroup1', ...
                             'Position', posButtonGroup1, ...
                             'BackgroundColor', rgb(colorPanelLeft), ...
                    'SelectionChangedFcn', @ButtonGroup1_SelectionChangedFcn);

% Create radio buttons for each way of ordering events
tagRadioButtons1 = {'Time', 'Amplitude', 'Rise', 'Decay'};
textRadioButtons1 = {'Time', 'Amp', 'Rise', 'Dec'};
for iButton = 1:numel(tagRadioButtons1)
    uicontrol(buttonGroup1, 'Style', 'radiobutton', ...
                           'Tag', tagRadioButtons1{iButton}, ...
                           'String', textRadioButtons1{iButton}, ...
                           'FontSize', radiobuttonFontSize, ...
                           'Position', posRadioButtons1{iButton}, ...
                           'HandleVisibility', 'off');
end

% Create an editable text box for choosing the window for examining events (ms)
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextZoomWindow', ...
                     'String', 'Zoom Window Size (ms):', ...
                     'Position', posTextZoomWindow);
editZoomWindow = ...
    uicontrol(panelLeft, 'Style', 'edit', 'Tag', 'EditZoomWindow', ...
                         'String', num2str(zoomWindowMs), ...
                         'Position', posEditZoomWindow, ...
                         'KeyPressFcn', @EditZoomWindow_KeyPressFcn, ...
                         'Callback', @EditZoomWindow_Callback);

% Create pushbuttons for incrementing/decrementing zoom window size (ms)
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonResetZoom', ...
                     'String', '---------------', ...
                     'Position', posButtonResetZoom, ...
                     'Callback', @ButtonResetZoom_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDownZoom', ...
                     'String', '-1', ...
                     'Position', posButtonDownZoom, ...
                     'Callback', @ButtonDownZoom_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUpZoom', ...
                     'String', '+1', ...
                     'Position', posButtonUpZoom, ...
                     'Callback', @ButtonUpZoom_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDown2Zoom', ...
                     'String', 'x.5', ...
                     'Position', posButtonDown2Zoom, ...
                     'Callback', @ButtonDown2Zoom_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUp2Zoom', ...
                     'String', 'x2', ...
                     'Position', posButtonUp2Zoom, ...
                     'Callback', @ButtonUp2Zoom_Callback);

% Create check box for isChecked status of zoomed event
%   Note: 'Enable' status will be updated after calling EditEventNo_Callback
checkboxChecked = ...
    uicontrol(panelLeft, 'Style', 'checkbox', 'Tag', 'CheckboxChecked', ...
                         'String', 'Event checked', ...
                         'FontSize', checkboxFontSize, ...
                         'Position', posCheckboxChecked, ...
                         'Callback', @CheckboxChecked_Callback, ...
                         'Enable', 'off');

% Create togglebutton for Clean Mode
buttonClean = ...
    uicontrol(panelLeft, 'Style', 'togglebutton', 'Tag', 'ButtonClean', ...
                         'String', 'Clean', ...
                         'Position', posButtonClean, ...
                         'Callback', @ButtonClean_Callback);    

% Create togglebutton for Add Event
uicontrol(panelLeft, 'Style', 'togglebutton', 'Tag', 'ButtonAdd', ...
                     'String', 'Add Event', ...
                     'Position', posButtonAdd, ...
                     'Callback', @ButtonAdd_Callback);

% Create an editable text box for choosing the range for cursor correction
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextCorrRange', ...
                     'String', 'Cursor Correction (ms):', ...
                     'Position', posTextCorrRange);

editCorrRange = ...
    uicontrol(panelLeft, 'Style', 'edit', 'Tag', 'EditCorrRange', ...
                         'String', num2str(corrRangeMs), ...
                         'Position', posEditCorrRange, ...
                         'KeyPressFcn', @EditCorrRange_KeyPressFcn, ...
                         'Callback', @EditCorrRange_Callback);

% Create an editable text box for choosing the event class to add
uicontrol(panelLeft, 'Style', 'text', 'Tag', 'TextAddClass', ...
                     'String', 'Class to Add/Change:', ...
                     'Position', posTextAddClass);

editAddClass = ...
    uicontrol(panelLeft, 'Style', 'edit', 'Tag', 'EditAddClass', ...
                         'String', num2str(classNoNewUser), ...
                         'Position', posEditAddClass, ...
                         'KeyPressFcn', @EditAddClass_KeyPressFcn, ...
                         'Callback', @EditAddClass_Callback);

% Create pushbutton for Remove Event
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonRemove', ...
                     'String', 'Remove', ...
                     'Position', posButtonRemove, ...
                     'Callback', @ButtonRemove_Callback);    

% Create pushbutton for Change Class
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonChange', ...
                     'String', 'Change', ...
                     'Position', posButtonChange, ...
                     'Callback', @ButtonChange_Callback);    

% Create check box for doForAll status
uicontrol(panelLeft, 'Style', 'checkbox', 'Tag', 'CheckboxForAll', ...
                     'String', 'All of class', ...
                     'FontSize', checkboxFontSize, ...
                     'Position', posCheckboxForAll, ...
                     'Callback', @CheckboxForAll_Callback, ...
                     'Enable', 'on');
                         
% Create check box for toPrompt status
checkboxNoPrompt = ...
    uicontrol(panelLeft, 'Style', 'checkbox', 'Tag', 'CheckboxNoPrompt', ...
                         'String', 'Don''t prompt', ...
                         'FontSize', checkboxFontSize, ...
                         'Position', posCheckboxNoPrompt, ...
                         'Callback', @CheckboxNoPrompt_Callback, ...
                         'Enable', 'on');

% Create check box for byForce status
uicontrol(panelLeft, 'Style', 'checkbox', 'Tag', 'CheckboxByForce', ...
                     'String', 'By Force', ...
                     'FontSize', checkboxFontSize, ...
                     'Position', posCheckboxByForce, ...
                     'Callback', @CheckboxByForce_Callback, ...
                     'Enable', 'off');
%                     'Enable', 'on');

% Create button group
buttonGroup2 = ...
    uibuttongroup(panelLeft, 'Tag', 'ButtonGroup2', ...
                             'Position', posButtonGroup2, ...
                             'BackgroundColor', rgb(colorPanelLeft));
%                  'SelectionChangedFcn', @buttonGroupSelection);

% Create radio buttons for each way of exiting GUI
tagRadioButtons2 = {'NextTrace', 'StartOver', 'Quit'};
textRadioButtons2 = {'Good! Next Trace!', 'Start Over ...', 'I Quit ><'};
for iButton = 1:numel(tagRadioButtons2)
    radioButtons2(iButton) = ...
            uicontrol(buttonGroup2, 'Style', 'radiobutton', ...
                                    'Tag', tagRadioButtons2{iButton}, ...
                                    'String', textRadioButtons2{iButton}, ...
                                    'FontSize', radiobuttonFontSize, ...
                                    'Position', posRadioButtons2{iButton}, ...
                                    'HandleVisibility', 'off');
end

% Create pushbutton for DONE
buttonDone = ...
    uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDone', ...
                         'String', 'DONE!', ...
                         'Position', posButtonDone, ...
                         'Callback', {@ButtonDone_Callback, buttonGroup2});

%% Create right panel
panelRight = ...
    uipanel(gui, 'Tag', 'PanelRight', 'Position', posPanelRight, ...
                 'Title', 'Current Sweep', ...
                 'BackgroundColor', rgb(colorPanelRight));

%{
% Create pushbutton for TODO
uicontrol(panelRight, 'Style', 'pushbutton', 'Tag', 'ButtonLast', ...
                     'String', 'Go to last event!', ...
                     'Position', posButtonLast, ...
                     'Callback', {@ButtonLast_Callback, buttonGroup2});
%}

% Create sweep figure
sweepFigure = ...
    axes('Parent', panelRight, 'Tag', 'SweepFigure', ...
         'ButtonDownFcn', @SweepFigure_ButtonDownFcn);
hPlot = minEASE_plot_gapfree_events(swpT, dataRaw, swpD, ...
                            eventInfo, eventClass, outputLabel, ...
                            'ClassNumbers', classNumbers, ...
                            'ClassColors', classColors, ...
                            'ClassLabels', classLabels, ...
                            'LegendLocation', 'northeast', ...
                            'MarkerSize', markerSize, ...
                            'LineWidth', markerEdgeWidth);

% Change xlimits to full time length
xlimits = [swpT(1) - siMs, swpT(end)];
xlim(xlimits);

% Change y-axis limits to 130s% of full data range
medianData = (minD + maxD) / 2;
halfRange = (maxD - minD) / 2;
ylim(medianData + 1.3 * halfRange * [-1, 1]);

% Do not allow the sweep figure to respond to left clicks for adding or removing 
%   by default
allowAddClick = false;
allowRemoveClick = false;

% Construct a context menu under the sweepFigure
%   that lets the user switch to pan mode or zoom mode by right-clicking
sweepFigureContextMenu = uicontextmenu;
set(sweepFigure, 'UIContextMenu', sweepFigureContextMenu);
uimenu('Parent', sweepFigureContextMenu, 'Label', 'Turn on pan mode',...
                'Callback', 'pan(gcbf, ''on'')');
uimenu('Parent', sweepFigureContextMenu, 'Label', 'Turn on zoom mode',...
                'Callback', 'zoom(gcbf, ''on'')');

% Construct a context menu (shows up by right-clicking) under panelRight
%   that lets the user switch to and from pan mode or zoom mode
panelRightContextMenu = uicontextmenu;
set(panelRight, 'UIContextMenu', panelRightContextMenu);
uimenu('Parent', panelRightContextMenu, 'Label', 'Turn on pan mode',...
                'Callback', 'pan(gcbf, ''on'')');
uimenu('Parent', panelRightContextMenu, 'Label', 'Turn on zoom mode',...
                'Callback', 'zoom(gcbf, ''on'')');
%{
%% TODO: The following is still not working under zoom or pan modes
%        Might need to turn interactive modes off from the MATLAB tool bar
uimenu('Parent', panelRightContextMenu, 'Label', 'Turn off pan mode', ...
                'Callback', 'pan(gcbf, ''off'')');
uimenu('Parent', panelRightContextMenu, 'Label', 'Turn off zoom mode', ...
                'Callback', 'zoom(gcbf, ''off'')');
%}

% TODO: Add the following to the default context menu 
%       under pan modes and zoom modes       
%{
hPan = pan;
panContextMenu = uicontextmenu;
uimenu('Parent', panContextMenu, 'Label', 'Turn off pan mode', ...
        'Callback', 'pan(gcbf, ''off'')');
uimenu('Parent', panContextMenu, 'Label', 'Switch to zoom mode', ...
        'Callback', 'zoom(gcbf, ''on'')');
set(hPan, 'UIContextMenu', panContextMenu);
%}
%{
hZoom = zoom;
zoomContextMenu = uicontextmenu;
uimenu('Parent', zoomContextMenu, 'Label', 'Switch to pan mode', ...
        'Callback', 'pan(gcbf, ''on'')');
uimenu('Parent', zoomContextMenu, 'Label', 'Turn off zoom mode', ...
        'Callback', 'zoom(gcbf, ''off'')');
set(hZoom, 'UIContextMenu', zoomContextMenu);
%}

% Set background and foreground colors for all objects in left panel
%   except pushbuttons
pLChildren = setdiff(setdiff(allchild(panelLeft), buttonGroup1), buttonGroup2);
bgChildren = union(allchild(buttonGroup1), allchild(buttonGroup2));
pLAllChildren = union(bgChildren, pLChildren);
pushbuttons = findall(gui, 'Style', 'pushbutton');
togglebuttons = findall(gui, 'Style', 'togglebutton');
bigbuttons = union(pushbuttons, togglebuttons);
allUiControl = union(pLAllChildren, pushbuttons);

set(setdiff(pLAllChildren, bigbuttons), ...
              'BackgroundColor', rgb(colorPanelLeft), ...
              'ForegroundColor', rgb(colorText));

% Set font sizes and units, etc. for pushbuttons
set(bigbuttons, 'FontUnits', textFontUnitsButton, ...
                'FontSize', textFontSizeButton, ...
                'FontName', textFontNameButton, ...
                'FontWeight', textFontWeightButton, ...
                'BackgroundColor', rgb(colorButton), ...
                'ForegroundColor', rgb(colorTextButton));

%% Initialize GUI
% Initialize the last examined overall event number to be the first one
eventNoInAllLast = 1;

% Set a global variable to 0 to allow for interruption
abort = false;

% Initialize other states
isIncrDecr = false;                     % whether incrementing or decrementing
                                        %   an edit field

% Extract maximum and minimum Y axis range from sweepFigure
ylimits = get(sweepFigure, 'Ylim'); % use the y axis limits in the 'all' plot
minY = ylimits(1); maxY = ylimits(2);

% Initialize arrays
%   Note: Initializing to gobjects(nEvents) takes a lot of time (4.3 seconds)
%           so try cell array instead. Not sure if this will be slower in the 
%           long term though
hPlot.peaksFilled = cell(nEvents);      % peaks that are checked and filled
hPlot.cursorMark = cell(4);             % cursor marks for adding events

% Initialize peaksFilled based on isChecked status
for iEvent = 1:nEvents
    if isChecked(iEvent)
        mark_event_checked(iEvent, checkboxChecked);
    end
end

% Initialize eventNoInAll
update_eventNoInAll(eventNoStr);

% Initialize thisClassInfo, nThisClass, eventRank, rankThisClass
ButtonGroup1_SelectionChangedFcn(buttonGroup1);
EditClassNo_Callback(editClassNo);

% Initialize doForAll, toPrompt and byForce stati
doForAll = false;
toPrompt = true;
byForce = false;

% Make GUI visible and print message in standard output
set(gui, 'Visible', 'on');
fprintf('minEASE (c) 2017 for %s is ready!\n\n', outputLabel);
drawnow;                                            % Update GUI

%% Restore default uicontrol properties
%   Note: these are set in minEASE_gui_examine_events.m
set(0, 'DefaultUicontrolFontUnits', 'remove', ...
       'DefaultUicontrolFontSize', 'remove', ...
       'DefaultUicontrolFontName', 'remove', ...
       'DefaultUicontrolFontWeight', 'remove', ...
       'DefaultUicontrolUnits', 'remove', ...
       'DefaultFigureUnits', 'remove');

% Return handle to GUI
hGUI = gui;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_eventNoInAll(eventNoStr, eventNo)
%% Update eventNoInAll and eventNoInAllLast
%   Note: eventNoInAll is the event number of the currently examined event
%         eventNoInAllLast is the event number of the last examined event

global eventNoInAll eventNoInAllLast                    % updated
global eventInfo                                        % used
global thisClassInfo                                    % used

if nargin < 2
    eventNo = str2double(eventNoStr);
end

if strcmpi(eventNoStr, 'all')               % to show all events
    % Set eventNoInAll to NaN
    eventNoInAll = NaN;
elseif isempty(eventNoStr)                  % to not update plot
    % Set eventNoInAll to 0
    eventNoInAll = 0;
else                                        % to show a particular event
    % Update the number of this event in all events
    eventNoInAll = find(eventInfo(:, 2) == thisClassInfo(eventNo, 2), 1);

    % Update the last examined overall event number
    eventNoInAllLast = eventNoInAll;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mark_event_checked(thisEventNoInAll, checkboxChecked)
%% Update isChecked status & fill circle of event peak in sweep figure
%   Set checkbox value only if checkboxChecked provided

global isChecked                                        % updated
global hPlot                                            % updated
global sweepFigure                                      % used
global swpT swpD                                        % used
global eventInfo eventClass                             % used
global classColorMap                                    % used
global markerSize markerEdgeWidth                       % used

if nargin >= 2
    set(checkboxChecked, 'Value', get(checkboxChecked, 'Max'));
end

% Update isChecked status
isChecked(thisEventNoInAll) = true;

% Find the class number and color of this event
classNo = eventClass(thisEventNoInAll);
thisColor = classColorMap(classNo, :);

% Fill circle of event peak in sweep figure
if isempty(hPlot.peaksFilled{thisEventNoInAll})
    % Find the peak index of this event
    idxPeak = eventInfo(thisEventNoInAll, 2);

    % Plot filled circle of event peak
    hPlot.peaksFilled{thisEventNoInAll} = ...
        plot(sweepFigure, swpT(idxPeak), swpD(idxPeak), ...
             'o', 'Color', thisColor, 'MarkerFaceColor', thisColor, ...
             'MarkerSize', markerSize, 'LineWidth', markerEdgeWidth, ...
             'DisplayName', ['Event #', num2str(thisEventNoInAll)], ...
             'HitTest', 'off');
end

% Update GUI
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mark_event_unchecked(thisEventNoInAll, checkboxChecked)
%% Update isChecked status & unfill circle of event peak in sweep figure
%   Set checkbox value only if checkboxChecked provided

global isChecked                                        % updated
global hPlot                                            % updated

if nargin >= 2
    set(checkboxChecked, 'Value', get(checkboxChecked, 'Min'));
end

% Update isChecked status
isChecked(thisEventNoInAll) = false;

% Unfill circle of event peak in sweep figure
delete(hPlot.peaksFilled{thisEventNoInAll});

% Reinitialize as empty array
hPlot.peaksFilled{thisEventNoInAll} = [];

% Update GUI
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zoom_event
%% Zoom to event by changing x-axis limits and replotting the event in red
%   Note: Assumes one of the following is true:
%       (1) eventNoStr is 'all' or '' (eventNoInAll is NaN or 0)
%       (2) eventNoInAll is a number within range
% TODO: Manually autoscale y axis here

global swpT swpD                                        % used
global minD maxD                                        % used
global nSamples siMs                                    % used
global zoomWindowMs                                     % used
global eventInfo                                        % used
global eventNoInAll                                     % used
global sweepFigure                                      % used

% Make sweep figure current figure
axes(sweepFigure);                                      % make current figure

% Find xlimits and corresponding ylimits
if isnan(eventNoInAll)                      % to show all events
    % Change xlimits to full data length
    xlimits = [swpT(1) - siMs, swpT(end)];

    % Find the data value extrema
    minDataNow = minD;
    maxDataNow = maxD;
elseif eventNoInAll == 0                    % to not update plot
    % Exit function
    return;
else                                        % to show a particular event
    % Convert window size for examining events to samples
    zoomWindowSamples = round(zoomWindowMs/siMs);

    % Find the peak index of this event
    idxPeak = eventInfo(eventNoInAll, 2);
    
    % The left limit is half of zoomWindowSamples from the peak index,
    %   or 0 if too small or nSamples - zoomWindowSamples if too large
    left = max(0, idxPeak - floor(zoomWindowSamples/2));
    if left + zoomWindowSamples > nSamples
        left = nSamples - zoomWindowSamples;
    end

    % Change xlimits to show exactly zoomWindowSamples samples
    xlimits = [swpT(left + 1), swpT(left + zoomWindowSamples)];

    % Find the data value extrema and median
    minDataNow = min(swpD((left + 1):(left + zoomWindowSamples)));
    maxDataNow = max(swpD((left + 1):(left + zoomWindowSamples)));
end

% Change the x-axis
xlim(xlimits);

% Change y-axis limits to 130% of full range
medDataNow = (minDataNow + maxDataNow) / 2;
halfRangeNow = (maxDataNow - minDataNow) / 2;
ylim(medDataNow + 1.3 * halfRangeNow * [-1, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_checkbox_status
%% Update to the correct checkbox status if a new event is selected
%   Note: Assumes one of the following is true:
%       (1) eventNoStr is 'all' or '' (eventNoInAll is NaN or 0)
%       (2) eventNoInAll is a number within range

global isChecked                                        % used
global eventNoInAll                                     % used
global checkboxChecked                                  % used

% Find checkbox status for examined event
if isnan(eventNoInAll) || eventNoInAll == 0         % no single event examined
    % Disable CheckboxChecked 
    set(checkboxChecked, 'Enable', 'off');
else                                                % single event selected
    % Update CheckboxChecked status for this event
    if isChecked(eventNoInAll)
        % Check checkbox
        set(checkboxChecked, 'Value', get(checkboxChecked, 'Max'));
    else
        % Uncheck checkbox
        set(checkboxChecked, 'Value', get(checkboxChecked, 'Min'));
    end

    % Enable CheckboxChecked 
    set(checkboxChecked, 'Enable', 'on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_gray_line(idxPeak)
% Draw a gray dashed line at a given peak index
%   Note: this is not always at the center of the screen, 
%         so annotation() is not useful here

global hPlot                                            % updated
global swpT                                             % used
global minY maxY                                        % used

% Use the new y axis limits for the gray dashed line
minLine = minY; maxLine = maxY;             % use the overall y axis range

% Put the line in hPlot
hPlot.grayLine = line([swpT(idxPeak), swpT(idxPeak)], ...
                      [minLine, maxLine], ...
                      'Color', rgb('Gray'), 'LineStyle', '--', ...
                      'DisplayName', 'Current peak time', ...
                      'HitTest', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function annotate_event
%% Annotate event information
%   Note: Assumes one of the following is true:
%       (1) eventNoStr is 'all' or '' (eventNoInAll is NaN or 0)
%       (2) eventNoInAll is a number within range

global hPlot                                            % used, then updated
global swpT swpD                                        % used
global siMs                                             % used
global eventInfo                                        % used
global eventNoInAll                                     % used
global posTextBox fontSizeTextBox fontWeightTextBox     % used
global textColorTextBox bgColorTextBox                  % used
global panelRight sweepFigure                           % used

% Make sweep figure current figure
axes(sweepFigure);                                      % make current figure

% Remove previous gray dashed line if exists
if isfield(hPlot, 'grayLine')
    delete(hPlot.grayLine);
end

% Remove previous replot if exists
if isfield(hPlot, 'thisData')
    delete(hPlot.thisData);
end

% Remove previous annotation text box if exists
if isfield(hPlot, 'textBox')
    delete(hPlot.textBox);
end

% If a new event is selected, 
%   (1) Plot a gray dashed line
%   (2) Replot event in red
%   (3) Display peak amplitude, 10-90% rise time and 50% decay time
if ~isnan(eventNoInAll) && eventNoInAll ~= 0        % an event is selected

    % Find the peak index of this event
    idxPeak = eventInfo(eventNoInAll, 2);

    % Draw a gray dashed line
    draw_gray_line(idxPeak);

    % Replot the event from breakpoint to "full decay point" in red
    %   If full decay point doesn't exist, 
    %   replot to just before next event breakpoint
    %   If event is removed, plot from breakpoint to peak only
    idxBreak = eventInfo(eventNoInAll, 1);       % index of event breakpoint
    fullDecayTime = eventInfo(eventNoInAll, 11); % full decay time (samples)
    interStimulusInterval = eventInfo(eventNoInAll, 9);
                % the number of samples from this peak to next event breakpoint
    idxEnd = get_idxEnd(idxPeak, fullDecayTime, interStimulusInterval);

    hPlot.thisData = plot(swpT(idxBreak:idxEnd), ...
                          swpD(idxBreak:idxEnd), ...
                          'r-', 'LineWidth', 2, ...
                          'DisplayName', 'Current event', ...
                          'HitTest', 'off');

    % Display annotation textbox with peak info
    hPlot.textBox = ...
        annotation(panelRight, 'textbox', posTextBox, ...
                   'String', {['Peak amplitude: ', ...
                            num2str(eventInfo(eventNoInAll, 5), 3), ...
                            ' pA         '], ...
                            ['10-90% rise time: ', ...
                            num2str(eventInfo(eventNoInAll, 7) * siMs, 3), ...
                            ' ms         '], ...
                            ['50% decay time: ', ...
                            num2str(eventInfo(eventNoInAll, 10) * siMs, 3), ...
                            ' ms         ']}, ...
                   'Color', rgb(textColorTextBox), ...
                   'FontSize', fontSizeTextBox, ...
                   'FontWeight', fontWeightTextBox, ...
                   'BackgroundColor', rgb(bgColorTextBox), ...
                   'FitBoxToText', 'on', ...
                   'HorizontalAlignment', 'left', ...
                   'VerticalAlignment', 'middle', ...
                   'HitTest', 'off');

    % Update GUI
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_peak_markers(classesChanged)

global hPlot                                            % updated
global swpT swpD                                        % used
global eventClass eventInfo                             % used
global classNumbers classColorMap classLabels           % used
global markerSize markerEdgeWidth                       % used

% Update event peak markers for any changed event class
for iClass = unique(classesChanged)
    % Remove previous plot of event peaks for this class
    if ~isa(hPlot.eventPeaks(iClass), 'matlab.graphics.GraphicsPlaceholder')
        % Delete plot
        delete(hPlot.eventPeaks(iClass));
    
        % Reinitialize as GraphicsPlaceholder
        hPlot.eventPeaks(iClass) = gobjects;
    end

    % Find peak indices of this class of events
    indPeaksIClass = eventInfo(eventClass == classNumbers(iClass), 2);

    % Plot peaks of this class of events
    h = plot(swpT(indPeaksIClass), swpD(indPeaksIClass), ...
             'o', 'Color', classColorMap(iClass, :), ...
             'MarkerSize', markerSize, 'LineWidth', markerEdgeWidth, ...
             'DisplayName', classLabels{iClass}, 'HitTest', 'off');
    if ~isempty(h)
        hPlot.eventPeaks(iClass) = h;
    end
end

% Refresh legend for valid plots
%   Note: hPlot.eventPeaks has dimensions 1, nClass
hAll = [hPlot.dataRaw, hPlot.dataLowpass, ...
        hPlot.eventBreaks, hPlot.eventPeaks];       % all plot handles
validPlots = arrayfun(@(x) ~isa(x, 'matlab.graphics.GraphicsPlaceholder'), ...
                      hAll);
hLegend = legend(hAll(validPlots));
set(hLegend, 'AutoUpdate', 'off');                  % For R2017a and beyond

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [classNoAuto, qString] = classify_PSC (eventNoThis, classNoUser, strMod)
% Automatically classify PSC and check if the same as user input

global eventClass eventInfo                             % used
global nEvents                                          % used

% Find the number of the next event, skipping all removed events
eventNoNext = eventNoThis + 1;
while eventNoNext ~= nEvents + 1 && eventClass(eventNoNext) == 8
    eventNoNext = eventNoNext + 1;
end

% Find the number of the previous event, skipping all removed events
eventNoPrev = eventNoThis - 1; 
while eventNoPrev ~= 0 && eventClass(eventNoPrev) == 8
    eventNoPrev = eventNoPrev - 1;
end

% Check if PSC is "incomplete":
%   interstimulus interval to next breakpoint < "full decay time"
if eventNoNext == nEvents + 1
    % Completeness of the last event cannot be determined
    isIncomplete = false;
else
    % If incomplete, "full decay time" cannot be calculated so is NaN
    isIncomplete = isnan(eventInfo(eventNoThis, 11));
end

% Check if PSC is "premature":
%   preceding interstimulus interval < preceding "full decay time"
if eventNoPrev == 0
    % First event cannot be premature
    isPremature = false;
else
    % In premature, the preceding "full decay time" is NaN
    isPremature = isnan(eventInfo(eventNoPrev, 11));
end

% Check if PSC is "Type II":
%   Incomplete PSCs that are not also premature
isTypeTwo = isIncomplete & ~isPremature;

% Check if PSC is "Type III":
%   Premature PSCs that follow Type II PSCs or other Type III PSCs
if eventNoPrev == 0
    % First event cannot be Type III
    isTypeThree = false;
else
    isTypeThree = isPremature & ...
                    (eventClass(eventNoPrev) == 2 | ...
                     eventClass(eventNoPrev) == 3);
end

% Decide on auto-classified class number
if isTypeTwo
    classNoAuto = 2;
    typeStr = 'Type II';
    reasonStr = ['it is incomplete (ISI < full decay) ', ...
                 'and NOT premature (previous event''s ISI', ...
                 ' < full decay).'];
elseif isTypeThree
    classNoAuto = 3;
    typeStr = 'Type III';
    reasonStr = ['it is premature (previous event''s ISI', ...
                 ' < full decay) and is preceded by ', ...
                 'a Type II PSC or a Type III PSC.'];
else
    classNoAuto = 1;
    typeStr = 'Type I';
    reasonStr = ['it is neither ''incomplete and NOT premature'' ', ...
                 'nor ''premature and preceded by a Type II/III PSC'''];
end        

% Check if user input is the same as automatically classified class number
if classNoAuto ~= classNoUser
    % Need to reclassify
    qString = {['Seems like the ', strMod, ...
                ' PSC is actually ', typeStr, ...
                ' because ', reasonStr], ...
               ['Would you like to change to ', typeStr, ...
                ' instead?']};
else
    % Don't need to reclassify
    qString = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [classNoNew, prevClassNoOld, prevClassNoAuto, nextClassNoOld, nextClassNoAuto] = change_class(eventNoThis, classNoNewUser, classNoOld)
% Change class number for given event to classNoNew
%   Note: assumes that eventNoThis is a number in range
%           and that classNoNewUser is different from classNoOld

global eventClass                                       % updated
global eventInfo                                        % used, then updated
global swpD                                             % used
global nEvents                                          % used
global editClassNo                                      % used
global checkboxChecked                                  % used
global eOrIFactor noiseLevel                            % used
global nSamples                                         % used
global nextClassNoOld                                   % used
global eventNoInAll                                     % used
global toPrompt                                         % used

% Check if the class numbers are different
if classNoNewUser == classNoOld
    error(['Something wrong with the code.\n', ...
            ' classNoNewUser should be different from classNoOld!']);
end

% Mark event unchecked (removes any previously filled marker)
mark_event_unchecked(eventNoThis, checkboxChecked);

% Find the number of the next event, skipping all removed events
eventNoNext = eventNoThis + 1;
while eventNoNext ~= nEvents + 1 && eventClass(eventNoNext) == 8
    eventNoNext = eventNoNext + 1;
end

% Find the number of the previous event, skipping all removed events
eventNoPrev = eventNoThis - 1; 
while eventNoPrev ~= 0 && eventClass(eventNoPrev) == 8
    eventNoPrev = eventNoPrev - 1;
end

% If changing FROM or TO class 8, update the next-event-dependent statistics:
%   interEventInterval, interStimulusInterval, halfDecayTime, fullDecayTime
%   for the current event and for the previous event
if (classNoOld ~= 8 && classNoNewUser == 8) || ...
    (classNoOld == 8 && classNoNewUser ~= 8)   % event was added or removed

    % Find breakpoint and peak indices for next event, or nSamples
    if eventNoNext == nEvents + 1        % if next event doesn't exist
        % Use nSamples
        idxPeakNext = nSamples;
        idxBreakNext = nSamples;
    else                                    % if next event does exist
        idxPeakNext = eventInfo(eventNoNext, 2);    % next breakpoint index
        idxBreakNext = eventInfo(eventNoNext, 1);   % next peak index
    end

    % If changing TO class 8, make the next-event-dependent statistics:
    %   interEventInterval, interStimulusInterval, halfDecayTime, fullDecayTime
    %   all NaN
    % If changing FROM class 8, compute the next-event-dependent statistics:
    %   interEventInterval, interStimulusInterval, halfDecayTime, fullDecayTime
    %   for this event
    if classNoOld ~= 8 && classNoNewUser == 8       % event is removed
        % Make interEventInterval, interStimulusInterval, 
        %   halfDecayTime, fullDecayTime all NaN
        eventInfo(eventNoThis, 8:11) = NaN * ones(1, 4);

        % Update displayed event information for this event
        %   if this event is the currently zoomed event
        if eventNoThis == eventNoInAll
            annotate_event;
        end
    elseif classNoOld == 8 && classNoNewUser ~= 8   % event is added
        % Compute next-event-dependent statistics for this event
        thisInfo = eventInfo(eventNoThis, :);  % event info for this event
        idxBreak = thisInfo(1);                 % breakpoint index
        idxPeak = thisInfo(2);                  % peak index
        thisInfo(8) = idxPeakNext - idxPeak;    % inter-event interval (samples)
        thisInfo(9) = idxBreakNext - idxPeak;% inter-stimulus interval (samples)
        halfDecayVal = thisInfo(4) - eOrIFactor * 0.5 * thisInfo(5);
                                                        % 50% decay value
        thisInfo(10) = find_custom(swpD(idxPeak:idxBreakNext) * eOrIFactor < ...
                                   halfDecayVal * eOrIFactor, 1, 'first', ...
                                   'ReturnNan', true) - 1;
                                                        % 50% decay time (samples)
        fullDecayVal = thisInfo(3) + eOrIFactor * noiseLevel;
                                                        % full decay value
        thisInfo(11) = find_custom(swpD(idxPeak:idxBreakNext) * eOrIFactor < ...
                                   fullDecayVal * eOrIFactor, 1, 'first', ...
                                   'ReturnNan', true) - 1;
                                                        % full decay time (samples)
        eventInfo(eventNoThis, :) = thisInfo;

        % Update displayed event information for this event
        %   if this event is the currently zoomed event
        if eventNoThis == eventNoInAll
            annotate_event;
        end
    end

    % Change the next-event-dependent statistics for the previous event
    if eventNoPrev ~= 0                     % if previous event exists
        % Extract info for the previous event
        idxPeakPrev = eventInfo(eventNoPrev, 2);    % previous peak index
        valBreakPrev = eventInfo(eventNoPrev, 3);   % prev breakpoint value
        valPeakPrev = eventInfo(eventNoPrev, 4);    % previous peak value
        ampPeakPrev = eventInfo(eventNoPrev, 5);    % previous peak amplitude

        % Find the next event of the previous event
        if classNoOld ~= 8 && classNoNewUser == 8       % event is removed
            % Use the next event
            idxPeakNextOfPrev = idxPeakNext;
            idxBreakNextOfPrev = idxBreakNext;
        elseif classNoOld == 8 && classNoNewUser ~= 8   % event is added
            % Use this event
            idxPeakNextOfPrev = idxPeak;
            idxBreakNextOfPrev = idxBreak;
        end

        % Compute the new interEventInterval, interStimulusInterval, 
        %   halfDecayTime, fullDecayTime for the previous event
        eventInfo(eventNoPrev, 8) = idxPeakNextOfPrev - idxPeakPrev;
        eventInfo(eventNoPrev, 9) = idxBreakNextOfPrev - idxPeakPrev;
        halfDecayValPrev = valPeakPrev - eOrIFactor * 0.5 * ampPeakPrev;
                                        % previous peak 50% decay value
        eventInfo(eventNoPrev, 10) = ...
            find_custom(swpD(idxPeakPrev:idxBreakNextOfPrev) * eOrIFactor < ...
                        halfDecayValPrev * eOrIFactor, 1, 'first', ...
                        'ReturnNan', true) - 1;     
                                        % previous peak 50% decay time (samples)
        fullDecayValPrev = valBreakPrev + eOrIFactor * noiseLevel;
                                        % previous peak full decay value
        eventInfo(eventNoPrev, 11) = ...
            find_custom(swpD(idxPeakPrev:idxBreakNextOfPrev) * eOrIFactor < ...
                        fullDecayValPrev * eOrIFactor, 1, 'first', ...
                        'ReturnNan', true) - 1;
                                        % prev peak full decay time (samples)
    end
end

% Decide on actual class number to add:
%   If changing to a PSC, check whether the type makes sense
%   and prompt user to change if not
if any(classNoNewUser == 1:3)   % if changing to a PSC (Classes 1~3)
    % Try automatically classifying the PSC and check if its the same
    [classNoAuto, qString] = ...
        classify_PSC(eventNoThis, classNoNewUser, 'current');

    % If classNoNewUser does not agree with auto-classification result,
    %   prompt user to decide whether to use auto-classified class number
    if ~isempty(qString)
        qTitle = 'Prompt for Class Decision';
        choice1 = 'YES, ok boss, whatever ...';
        choice2 = 'NO, listen to what I''m telling ya!';
        if toPrompt
            answer = questdlg(qString, qTitle, choice1, choice2, choice1);
        else
            answer = choice1;
        end
        switch answer
        case choice1
            % New class number is the auto-classified class number
            classNoNew = classNoAuto;
        case choice2
            % New class number is the user input
            classNoNew = classNoNewUser;
        otherwise
            error('Problem with code!');
        end
    else
        % New class number is the same as user input
        classNoNew = classNoNewUser;
    end
else                            % if not changing to a PSC
    % New class number is the same as user input
    classNoNew = classNoNewUser;
end

if classNoNewUser == classNoOld
    error(['Something wrong with the code.\n', ...
            ' classNoNewUser should be different from classNoOld!']);
end

% Only continue if class number is indeed changed
if classNoNew == classNoOld
    prevClassNoOld = [];
    prevClassNoAuto = [];
    nextClassNoOld = [];
    nextClassNoAuto = [];
    return
end

% Change eventClass for this event to new value
eventClass(eventNoThis) = classNoNew;

% Determine the need to reclassify the previous event according to the 
%   following truths derived from the definition of Type II and III PSCs:
%   (1) Type II PSCs never precede a Type I PSC or a Type II PSC,
%           but it could precede non-PSC events
%   (2) Type III PSCs always follow a Type II PSC or a Type III PSC
reclassifyPrevFlag = false;             % do not reclassify by default
if eventNoPrev ~= 0                     % a previous event exists
    % Retrieve the previous event's class number
    prevClassNoOld = eventClass(eventNoPrev);

    % If the previous event is a Type II PSC, truth #1 is violated 
    %   if one of the following is true:
    %   (1) this event is forcibly changed to a Type I PSC or a Type II PSC.
    %   (2) this event is removed or added and the next event occurs 
    %           before the previous event decays completely
    if prevClassNoOld == 2
        qString = {};           % initialize prompt as an empty string
        if any(classNoNew == [1, 2])
                                % user decides this PSC is Type I or Type II
            %   Since a Type II PSC is a PSC (not Class 4~8) and not premature 
            %       (not Class 3), the previous PSC should really be Type I
            prevClassNoAuto = 1;

            % Set up prompt
            qString = {'Type II PSCs never precede a Type I or II PSC.', ...
                       'Recommend changing previous event to a Type I PSC.', ...
                       'Proceed?'};
        elseif classNoNew == 8 || classNoOld == 8  
                                % user adds or removes this event
            % Try automatically classifying the previous PSC
            [prevClassNoAuto, qString] = ...
                classify_PSC(eventNoPrev, prevClassNoOld, 'previous');
        else
            % Don't need to reclassify previous event
        end

        % Only prompt if needed
        if ~isempty(qString)
            % Prompt user to decide whether to change class of previous event
            qTitle = 'Prompt for Previous Event Class Decision';
            choice1 = 'YES, change to this class';
            choice2 = 'NO, I want my old class!';
            if toPrompt
                answer = questdlg(qString, qTitle, choice1, choice2, choice1);
            else
                answer = choice1;
            end
            switch answer
            case choice1
                % Reclassify previous event
                reclassifyPrevFlag = true;
            case choice2
                % Don't reclassify previous event
            otherwise
                error('Problem with code!');
            end
        else
            % Don't need to reclassify previous peak
        end
    else
        % Don't need to reclassify previous peak
    end
else
    % No previous event to reclassify
end

% Reclassify the previous event
if reclassifyPrevFlag
    % New class number is the auto-classified class number
    eventClass(eventNoPrev) = prevClassNoAuto;
    
    % Mark event unchecked
    mark_event_unchecked(eventNoPrev);

    % Update GUI
    drawnow;
else
    % Don't need to reclassify previous peak, set prevClassNoOld 
    %    and prevClassNoAuto to be empty to prevent plot update
    prevClassNoOld = [];
    prevClassNoAuto = [];
end

% Determine the need to reclassify the next event according to the truths above
reclassifyNextFlag = false;             % do not reclassify by default
if eventNoNext ~= nEvents + 1           % a next event exists
    % Retrieve the next event's class number
    nextClassNoOld = eventClass(eventNoNext);

    % If the next event is a Type I PSC or a Type II PSC, truth #1 is violated 
    %   if one of the following is true:
    %   (1) this event is forcibly changed to a Type II PSC.
    %   (2) this event is removed or added and the next event occurs 
    %           before the previous event decays completely
    % If the next event is a Type III PSC, truth #2 is violated 
    %   if one of the following is true: 
    %   (1) this event is forcibly changed to NOT a Type II PSC
    %       nor a Type III PSC
    %   (2) this event is removed or added and the next event occurs 
    %           after the previous event decays completely
    qString = {};               % initialize prompt as an empty string
    if any(nextClassNoOld == [1, 2])
        if classNoNew == 2      % user decides this PSC is Type II
            % Since a Type I or II PSC is a PSC (not Class 4~8), 
            %   the next PSC should really be a Type III PSC
            nextClassNoAuto = 3;

            % Set up prompt
            qString = {'Type II PSCs never precede a Type I or II PSC.', ...
                       'Recommend changing next event to a Type III PSC.', ...
                       'Proceed?'};
        elseif classNoNew == 8 || classNoOld == 8  
                                % user adds or removes this event
            % Try automatically classifying the next PSC
            [nextClassNoAuto, qString] = ...
                classify_PSC(eventNoNext, nextClassNoOld, 'next');
        else
            % Don't need to reclassify next event
        end
    elseif nextClassNoOld == 3
        if ~any(classNoNew == [2, 3])
            %   Since a Type III PSC is a PSC (not Class 4~8), 
            %       the next PSC should really be a Type I PSC or a Type II PSC
            if eventNoThis ~= nEvents - 1 && ...   % if not second to last event
                eventClass(eventNoThis + 2) == 3   % if next next event is 
                                                    %   a Type III PSC
                % The next PSC should really be a Type II PSC
                nextClassNoAuto = 2;
                typeStr = 'Type II';
            else
                % The next PSC should really be a Type I PSC
                nextClassNoAuto = 1;
                typeStr = 'Type I';
            end
            qString = {['Type III PSCs always follow a ', ...
                        'Type II PSC or a Type III PSC.'], ...
                       ['Recommend changing next event to a ', ...
                        typeStr, ' PSC.'], 'Proceed?'};
            % Don't need to reclassify next peak
        elseif classNoNew == 8 || classNoOld == 8  
                                % user adds or removes this event
            % Try automatically classifying the next PSC
            [nextClassNoAuto, qString] = ...
                classify_PSC(eventNoNext, nextClassNoOld, 'next');
        else        
            % Don't need to reclassify next event
        end
    else
        % Don't need to reclassify next peak
    end

    % Only prompt if needed
    if ~isempty(qString)
        % Prompt user to decide whether to change class of next event
        qTitle = 'Prompt for Next Event Class Decision';
        choice1 = 'YES, change to this class';
        choice2 = 'NO, I want my old class!';
        if toPrompt
            answer = questdlg(qString, qTitle, choice1, choice2, choice1);
        else
            answer = choice1;
        end
        switch answer
        case choice1
            % Reclassify next event
            reclassifyNextFlag = true;
        case choice2
            % Don't reclassify next event
        otherwise
            error('Problem with code!');
        end
    else
        % Don't need to reclassify next peak
    end
else
    % No next event to reclassify
end

% Reclassify the next event
if reclassifyNextFlag
    % New class number is the auto-classified class number
    eventClass(eventNoNext) = nextClassNoAuto;
    
    % Mark event unchecked
    mark_event_unchecked(eventNoNext);

    % Update GUI
    drawnow;
else
    % Don't need to reclassify next peak, set nextClassNoOld 
    %    and nextClassNoAuto to be empty to prevent plot update
    nextClassNoOld = [];
    nextClassNoAuto = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function remove_by_force(eventNoToRemove)
%% Remove an event by force
% TODO

global eventInfo                                        % updated
global isChecked                                        % updated
global eventClass                                       % used, then updated

% Retrieve old class number
classNoOld = eventClass(eventNoToRemove);

% Remove the event from eventInfo, eventClass and isChecked
% TODO

% Update nEvents and nEventsThisClass and its display
% TODO

% Update eventNoInAll
% TODO

% Update zoomed event
% TODO

% Update properties of previous event
% TODO

% Update the markers for those classes that are changed
classesChanged = classNoOld;
update_peak_markers(classesChanged);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function remove_event(eventNoToRemove)
%% Remove an event

global eventClass                                       % used, then updated
global editClassNo                                      % used
global doForAll                                         % used
global toPrompt                                         % used
global byForce                                          % used
global eventInfo                                        % used
global checkboxChecked checkboxNoPrompt                 % used
% global eventNoInAll                                     % used

% Prompt user if about to remove by force
if toPrompt && byForce
    % Display question dialogue box: 
    answer = questdlg({'Are you sure about removing this event permanently?', ...
                        'You can uncheck the ''By Force'' checkbox if not.'}, ...
                      'Prompt for Remove Event by Force', 'Yes');
    switch answer
    case 'Yes'
        % Proceed
    case {'No', 'Cancel'}
        % Do nothing
        return;
    otherwise
        error('Problem with code!');
    end
end

% Retrieve old class number
classNoOld = eventClass(eventNoToRemove);

% Do nothing if the old class is already 8 and user does not to remove
%   the event by force
if ~byForce && classNoOld == 8
    % Display error dialog box and return
    errormsg = 'Event already removed!';
    dlgname = 'Remove Event Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
    return;
end

% Change class number for the event or all events of same class to 8
if ~byForce && doForAll
    % Do not prompt: %TODO: remove this when prompt is implemented with rezooming
    set(checkboxNoPrompt, 'Value', get(checkboxNoPrompt, 'Max'));
    CheckboxNoPrompt_Callback(checkboxNoPrompt);
    
    % Get the event numbers of all events of this class
    eventNosToChange = find(eventClass == classNoOld);

    % Get the number of events to change
    nEventsToChange = length(eventNosToChange);

    % Change class number for all events of this class to 8
    classesChanged = [];
    for iEventsToChange = 1:nEventsToChange
        % Retrieve the number of this event in all events
        eventNoThis = eventNosToChange(iEventsToChange);

        % Change the class of this event
        %   Note: Must follow change_class() by update_peak_markers()
        [classNoNew, prevClassNoOld, prevClassNoAuto, ...
            nextClassNoOld, nextClassNoAuto] = ...
            change_class(eventNoThis, 8, classNoOld);

        % Add to classesChanged
        classesChanged = [classesChanged, classNoOld, classNoNew, ...
                             prevClassNoOld, prevClassNoAuto, ...
                             nextClassNoOld, nextClassNoAuto];

        % TODO: Mark event checked?
        % mark_event_checked(eventNoThis, checkboxChecked);
    end

    % Update the markers for those classes that are changed
    update_peak_markers(classesChanged);

    % If event is currently zoomed, 
    %   set displayed class number to 8 
    % TODO: When to do this?
    % if eventNoToRemove == eventNoInAll
        % set(editClassNo, 'String', num2str(8));
    % end
elseif ~byForce
    % Change class number for this event to 8
    %   Note: Must follow change_class() by update_peak_markers()
    [classNoNew, prevClassNoOld, prevClassNoAuto, ...
        nextClassNoOld, nextClassNoAuto] = ...
        change_class(eventNoToRemove, 8, classNoOld);

    % Put all the changed classes together 
    %   and update the markers for those that are changed
    classesChanged = [classNoOld, classNoNew, ...
                         prevClassNoOld, prevClassNoAuto, ...
                         nextClassNoOld, nextClassNoAuto];
    update_peak_markers(classesChanged);

    % Mark event checked
    mark_event_checked(eventNoToRemove, checkboxChecked);

    % If event is currently zoomed, 
    %   set displayed class number to old class number
    %   so that user can move to the next event in the class of interest
    % TODO: When to do this?
    % if eventNoToRemove == eventNoInAll
        % set(editClassNo, 'String', num2str(classNoOld));
    % end
elseif doForAll
    % Remove all events of this class by force
    % TODO
else
    % Remove this event by force
    remove_by_force(eventNoToRemove)
end

% Execute callback function for EditClassNo 
%   (this updates thisClassInfo, rankThisClass & nThisClass and 
%       update xlimits, checkbox, annotation if event is changed)
EditClassNo_Callback(editClassNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Callbacks for minEASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MainGUI_WindowScrollWheelFcn(~, callbackData)
%% Executes on mouse scroll in main GUI

global nSamples siMs                                    % used
global editZoomWindow                                   % used
global editEventNo                                      % used

% Increment zoom window size and execute callback function for EditZoomWindow
increment_editbox(editZoomWindow, 1, nSamples * siMs, ...
                  10 * callbackData.VerticalScrollCount, {});
EditZoomWindow_Callback(editZoomWindow);

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonAllClass_Callback(~, ~)
%% Executes on press of ButtonAllClass

global editClassNo                                      % used

% Set class number string to 'all'
set(editClassNo, 'String', 'all');

% Execute callback function for EditClassNo
EditClassNo_Callback(editClassNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDownClass_Callback(~, ~)
%% Executes on press of ButtonDownClass

global isIncrDecr                                       % updated
global nClass                                           % used
global editClassNo                                      % used

% Decrement class number
increment_editbox(editClassNo, 1, nClass, -1, {'all'});

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditClassNo
EditClassNo_Callback(editClassNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonUpClass_Callback(~, ~)
%% Executes on press of ButtonUpClass

global isIncrDecr                                       % updated
global nClass                                           % used
global editClassNo                                      % used

% Increment class number
increment_editbox(editClassNo, 1, nClass, 1, {'all'});

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditClassNo
EditClassNo_Callback(editClassNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditClassNo_KeyPressFcn(hObject, callbackData)
%% Executes on key press in EditClassNo

global nClass                                           % used
global isIncrDecr                                       % updated

% If up arrow/pageup or down arrow/pagedown is pressed, 
%   increment or decrement, respectively. 
switch callbackData.Key
case {'downarrow', 'pagedown'}
    increment_editbox(hObject, 1, nClass, -1, {'all'});
case {'uparrow', 'pageup'}
    increment_editbox(hObject, 1, nClass, 1, {'all'});
otherwise
    % Otherwise, do nothing
    return;
end

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditClassNo
EditClassNo_Callback(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditClassNo_Callback(hObject, ~)
%% Executes on text change in EditClassNo
%   Note: will also execute callback function for EditEventNo, 
%           which also updates xlimits, checkbox, annotation

global thisClassInfo nThisClass rankThisClass           % updated, then used
global isIncrDecr                                       % used, then updated
global eventInfo eventClass eventRank                   % used 
global nClass                                           % used
global eventNoInAllLast                                 % used 
global editEventNo dispNThisClass                       % used
global toPrompt                                         % used

% Update classNoStr
classNoStr = get(hObject, 'String');

% Convert classNoStr to a number
%   Note: str2double() will return NaN if not a number
classNo = str2double(classNoStr);

% If not within range, display error and force change to 'all'
if ~strcmpi(classNoStr, 'all') && ...
    (isnan(classNo) || mod(classNo, 1) ~= 0 || ...
    classNo < 1 || classNo > nClass)
    % Display error dialog box
    errormsg = sprintf(['Class number must either be ''all'' ', ...
                          'or a positive integer between 1 and %d!'], ...
                          nClass);
    dlgname = 'Class Number Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));

    % Force change to 'all'
    classNoStr = 'all';
    classNo = NaN;
    set(hObject, 'String', classNoStr);
end

% Update thisClassInfo & rankThisClass
%   NOTE: thisClassInfo is the eventInfo of the class that is displayed,
%           not necessarily the class of the event that is selected
%   NOTE: rankThisClass is the eventRank of the class that is displayed
%           not necessarily the class of the event that is selected
if strcmpi(classNoStr, 'all')
    thisClassInfo = eventInfo;
    rankThisClass = eventRank;
else
    thisClassInfo = eventInfo(eventClass == classNo, :);
    if ~isempty(eventRank)
        rankThisClass = eventRank(eventClass == classNo);
    else
        rankThisClass = [];
    end
end

% Update number of events in this class
%   NOTE: nThisClass is the number of events of the class that is displayed
nThisClass = size(thisClassInfo, 1);        % number of events of this class

% Update dispNThisClass
set(dispNThisClass, 'String', num2str(nThisClass));

% Give EditEventNo focus when not incrementing/decrementing
if isIncrDecr
    % Give focus to hObject (editClassNo)
    uicontrol(hObject);

    % Don't restore isIncrDecr state yet because EditEventNo_Callback
    %   will be called
else
    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(editEventNo);
end

% If event number string is not 'all' (trying to examine events), 
%   find the event of this new class closest to the last examined event
if ~strcmpi(get(editEventNo, 'String'), 'all')
    if nThisClass == 0                            
        % If there are no events in this class, prompt user for what to do
        qString = 'No events are in this class! What to do?';
        qTitle = 'Prompt for Class Number';
        choice1 = 'Ignore';
        choice2 = 'Return to last examined event';
        choice3 = 'Show all events';
        if toPrompt
            answer = questdlg(qString, qTitle, ...
                              choice1, choice2, choice3, choice1);
        else
            answer = choice1;
        end
        switch answer
        case choice1
            % Update the event number
            set(editEventNo, 'String', '');

            % Execute callback function for EditEventNo
            EditEventNo_Callback(editEventNo);
        case choice2
            % Update the class number to last examined class
            set(hObject, 'String', eventClass(eventNoInAllLast));

            % Execute callback function for EditClassNo again
            %   Warning: this recursive usage is ok 
            %            only if no more code follow the loop
            EditClassNo_Callback(hObject);
        case choice3
            % Update the event number
            set(editEventNo, 'String', 'all');

            % Execute callback function for EditEventNo
            EditEventNo_Callback(editEventNo);
        end
    else
        % Find peak index of last examined event
        idxLastPeak = eventInfo(eventNoInAllLast, 2);
        
        % Find all peak indices of new class of events
        idxThisPeaks = thisClassInfo(:, 2);

        % Find the event of the new class with peak index closest to original
        [~, eventNoNew] = min(abs(idxThisPeaks - idxLastPeak));

        % Update the event number
        set(editEventNo, 'String', num2str(eventNoNew));

        % Execute callback function for EditEventNo (includes updating xlimits)
        EditEventNo_Callback(editEventNo);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonAllNo_Callback(~, ~)
%% Executes on press of ButtonAllNo

global editEventNo                                      % used

% Set event number string to 'all'
set(editEventNo, 'String', 'all');

% Execute callback function for EditEventNo
EditEventNo_Callback(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDownNo_Callback(~, ~)
%% Executes on press of ButtonDownNo

global nThisClass rankThisClass                         % used
global editEventNo                                      % used

% Decrement event number
increment_editbox(editEventNo, 1, nThisClass, -1, {'all', ''}, ...
                    'ValueNo', rankThisClass);

% Execute callback function for EditEventNo
EditEventNo_Callback(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonUpNo_Callback(~, ~)
%% Executes on press of ButtonUpNo

global nThisClass rankThisClass                         % used
global editEventNo                                      % used

% Increment event number
increment_editbox(editEventNo, 1, nThisClass, 1, {'all', ''}, ...
                    'ValueNo', rankThisClass);

% Execute callback function for EditEventNo
EditEventNo_Callback(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDown2No_Callback(~, ~)
%% Executes on press of ButtonDown2No

global nThisClass rankThisClass                         % used
global editEventNo checkboxChecked                      % used
global eventNoInAll                                     % used

% If coming from an examined event, mark it checked
if ~isnan(eventNoInAll) && eventNoInAll ~= 0
    mark_event_checked(eventNoInAll, checkboxChecked);
end

% Decrement event number
increment_editbox(editEventNo, 1, nThisClass, -1, {'all', ''}, ...
                    'ValueNo', rankThisClass);

% Execute callback function for EditEventNo
EditEventNo_Callback(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonUp2No_Callback(~, ~)
%% Executes on press of ButtonUp2No

global nThisClass rankThisClass                         % used
global editEventNo checkboxChecked                      % used
global eventNoInAll                                     % used

% If coming from an examined event, mark it checked
if ~isnan(eventNoInAll) && eventNoInAll ~= 0
    mark_event_checked(eventNoInAll, checkboxChecked);
end

% Increment event number
increment_editbox(editEventNo, 1, nThisClass, 1, {'all', ''}, ...
                    'ValueNo', rankThisClass);

% Execute callback function for EditEventNo
EditEventNo_Callback(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditEventNo_KeyPressFcn(hObject, callbackData)
%% Executes on key press in EditEventNo
%   Differs from button presses in that checked events are skipped

global nThisClass rankThisClass                         % used
global checkboxChecked                                  % used
global eventNoInAll                                     % used
global eventInfo isChecked                              % used
global thisClassInfo                                    % used

% If page up/down is pressed and 
%   if coming from an examined event, mark it checked
if any(strcmpi(callbackData.Key, {'pagedown', 'pageup'})) && ...
    ~isnan(eventNoInAll) && eventNoInAll ~= 0
    mark_event_checked(eventNoInAll, checkboxChecked);
end

% If up arrow/pageup or down arrow/pagedown is pressed, 
%   increment or decrement, respectively, until an unchecked event is reached. 
switch callbackData.Key
case {'downarrow', 'pagedown'}
    % Decrement once
    newStr = increment_editbox(hObject, 1, nThisClass, -1, {'all', ''}, ...
                                'ValueNo', rankThisClass);

    % Decrement more if the new event exists and is already checked
    newNo = str2double(newStr);
    while ~isnan(newNo) && ~isempty(thisClassInfo) && ...
            isChecked(find(eventInfo(:, 2) == thisClassInfo(newNo, 2), 1))
        % Decrement once more
        newStr = increment_editbox(hObject, 1, nThisClass, -1, {'all', ''}, ...
                                    'ValueNo', rankThisClass);
        newNo = str2double(newStr);        
    end
case {'uparrow', 'pageup'}
    % Increment once
    newStr = increment_editbox(hObject, 1, nThisClass, 1, {'all', ''}, ...
                                'ValueNo', rankThisClass);

    % Increment more if the new event exists and is already checked
    newNo = str2double(newStr);
    while ~isnan(newNo) && ~isempty(thisClassInfo) && ...
            isChecked(find(eventInfo(:, 2) == thisClassInfo(newNo, 2), 1))
        % Increment once more
        newStr = increment_editbox(hObject, 1, nThisClass, 1, {'all', ''}, ...
                                    'ValueNo', rankThisClass);
        newNo = str2double(newStr);        
    end
case {'delete'}
    % Remove event
    ButtonRemove_Callback;

    % Increment to next unchecked event
    newStr = get(hObject, 'String');
    newNo = str2double(newStr);
    while ~isnan(newNo) && ...
            isChecked(find(eventInfo(:, 2) == thisClassInfo(newNo, 2), 1))
        % Increment once more
        newStr = increment_editbox(hObject, 1, nThisClass, 1, {'all', ''}, ...
                                    'ValueNo', rankThisClass);
        newNo = str2double(newStr);        
    end
case {'insert'}
    % Change event
    ButtonChange_Callback;

    % Increment to next unchecked event
    newStr = get(hObject, 'String');
    newNo = str2double(newStr);
    while ~isnan(newNo) && ...
            isChecked(find(eventInfo(:, 2) == thisClassInfo(newNo, 2), 1))
        % Increment once more
        newStr = increment_editbox(hObject, 1, nThisClass, 1, {'all', ''}, ...
                                    'ValueNo', rankThisClass);
        newNo = str2double(newStr);        
    end
otherwise
    % Otherwise, do nothing
    return;
end

% Execute callback function for EditEventNo
EditEventNo_Callback(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditEventNo_Callback(hObject, ~)
%% Executes on text change in EditEventNo

global isIncrDecr                                       % used, then updated
global nThisClass                                       % used
global eventNoInAll eventNoInAllLast                    % used

% Update text of EditEventNo
eventNoStr = get(hObject, 'String');

% Convert to a number
%   Note: str2double() will return NaN if not a number
eventNo = str2double(eventNoStr);

% If not a number, validate string
if isnan(eventNo) && ~isempty(eventNoStr)
    % Validate eventNoStr and return empty string if not a match
    eventNoStr = validate_string(eventNoStr, {'all'});

    % If nonsense, change eventNoStr to 'nonsense'
    if isempty(eventNoStr)
        eventNoStr = 'nonsense';
    end
    
    % Update GUI
    set(hObject, 'String', eventNoStr);
end

% If eventNoStr is out of range, display error and 
%   set eventNoStr to either a valid value or ''
if ~strcmpi(eventNoStr, 'all') && ~isempty(eventNoStr) && ...
    (strcmpi(eventNoStr, 'nonsense') || ... % gibberish
     nThisClass == 0 || ...                 % none in this class
     ~isnan(eventNo) && (mod(eventNo, 1) ~= 0 || ...
                         eventNo < 1 || eventNo > nThisClass))
                                            % out of range
    % Display error dialog box
    if nThisClass == 0
        errormsg = 'There are no events in this class!';
    else
        errormsg = sprintf(['Event number must either be ''all'' ', ...
                              'or a positive integer between 1 and %d!'], ...
                              nThisClass);
    end
    dlgname = 'Event Number Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
                 
    % Set eventNoStr according to user input
    if strcmpi(eventNoStr, 'nonsense') || nThisClass == 0    
                                        % if not a number or none in this class
        % Set eventNoStr to ''
        eventNoStr = '';
    else
        if mod(eventNo, 1) ~= 0         % if not an integer
            % Find closest integer
            closestInteger = round(eventNo);
            
            % Force eventNo within range
            eventNo = min(nThisClass, max(1, closestInteger));            
        elseif eventNo < 1                      % if too small
            % Set eventNo to minimum value
            eventNo = 1;
        elseif eventNo > nThisClass             % if too large
            % Set eventNo to maximum value
            eventNo = nThisClass;
        end
        % Update eventNoStr
        eventNoStr = num2str(eventNo);
    end
    
    % Update GUI
    set(hObject, 'String', eventNoStr);
end

% Save previous instance of eventNoInAll & eventNoInAllLast
eventNoInAllPrev = eventNoInAll;
eventNoInAllLastPrev = eventNoInAllLast;

% Update eventNoInAll & eventNoInAllLast
update_eventNoInAll(eventNoStr, eventNo);

% Update event display if the examined event is changing from 'all' or ''
%   or if the new event is not the same as the last examined event
if isnan(eventNoInAllPrev) || eventNoInAllPrev == 0 || ...
    eventNoInAll ~= eventNoInAllLastPrev
    % Update xlimits of sweep figure
    zoom_event;

    % Update to the correct checkbox status
    update_checkbox_status;

    % Annotate event information
    annotate_event;
end

% Give hObject (editEventNo) focus if not incrementing/decrementing
%   from elsewhere
if isIncrDecr
    % Restore isIncrDecr state
    isIncrDecr = false;
else
    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(hObject);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonGroup1_SelectionChangedFcn(hObject, callbackData)
%% Executes on selection change of radiobuttons in button group 1

global eventRank rankThisClass                          % updated
global nEvents eventInfo eventClass                     % used
global editClassNo                                      % used
global editEventNo                                      % used

% Get radiobutton selection (what to rank events by)
rankBy = hObject.SelectedObject.Tag;

% Update eventRank
switch rankBy
case 'Time'
%    % Events in eventInfo are already ranked by time
%    eventRank = (1:nEvents)';
    % Make eventRank empty
    eventRank = [];
case 'Amplitude'
    % Sort the events in eventInfo by amplitude in ascending order
    [~, origIndices] = sort(eventInfo(:, 5));

    % Give a rank to each event
    eventRank(origIndices) = (1:nEvents)';
case 'Rise'
    % Sort the events in eventInfo by 10-90% rise time in ascending order
    [~, origIndices] = sort(eventInfo(:, 7));

    % Give a rank to each event
    eventRank(origIndices) = (1:nEvents)';
case 'Decay'
    % Sort the events in eventInfo by 50% decay time in ascending order
    [~, origIndices] = sort(eventInfo(:, 10));

    % Give a rank to each event
    eventRank(origIndices) = (1:nEvents)';
end

% Update rankThisClass
if ~isempty(eventRank)
    classNoStr = get(editClassNo, 'String');
    classNo = str2double(classNoStr);
    if strcmpi(classNoStr, 'all')
        rankThisClass = eventRank;
    else
        rankThisClass = eventRank(eventClass == classNo);
    end
else
    rankThisClass = [];
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonResetZoom_Callback(~, ~)
%% Executes on press of ButtonResetZoom

global zoomWindowInit                                   % used
global editZoomWindow                                   % used

% Reset zoom window size
set(editZoomWindow, 'String', num2str(zoomWindowInit));

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(editZoomWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDownZoom_Callback(~, ~)
%% Executes on press of ButtonDownZoom

global isIncrDecr                                       % updated
global nSamples siMs                                    % used
global editZoomWindow                                   % used

% Decrement zoom window size
increment_editbox(editZoomWindow, 1, nSamples * siMs, -1, {});

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(editZoomWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonUpZoom_Callback(~, ~)
%% Executes on press of ButtonUpZoom

global isIncrDecr                                       % updated
global nSamples siMs                                    % used
global editZoomWindow                                   % used

% Increment zoom window size
increment_editbox(editZoomWindow, 1, nSamples * siMs, 1, {});

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(editZoomWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDown2Zoom_Callback(~, ~)
%% Executes on press of ButtonDown2Zoom

global isIncrDecr                                       % updated
global nSamples siMs                                    % used
global editZoomWindow                                   % used

% Decrement zoom window size
increment_editbox(editZoomWindow, 0, nSamples * siMs, -2, {}, ...
                  'Operation', 'multiply');

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(editZoomWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonUp2Zoom_Callback(~, ~)
%% Executes on press of ButtonUp2Zoom

global isIncrDecr                                       % updated
global nSamples siMs                                    % used
global editZoomWindow                                   % used

% Increment zoom window size
increment_editbox(editZoomWindow, 0, nSamples * siMs, 2, {}, ...
                  'Operation', 'multiply');

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(editZoomWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditZoomWindow_KeyPressFcn(hObject, callbackData)
%% Executes on key press in EditZoomWindow

global isIncrDecr                                       % updated
global nSamples siMs                                    % used

% If up or down arrow is pressed, increment or decrement, respectively. 
switch callbackData.Key
case 'downarrow'
    increment_editbox(hObject, 1, nSamples * siMs, -5, {});
case 'uparrow'
    increment_editbox(hObject, 1, nSamples * siMs, 5, {});
case 'pagedown'
    increment_editbox(hObject, 1, nSamples * siMs, -50, {});
case 'pageup'
    increment_editbox(hObject, 1, nSamples * siMs, 50, {});
otherwise
    % Otherwise, do nothing
    return;
end

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditZoomWindow
EditZoomWindow_Callback(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditZoomWindow_Callback(hObject, ~)
%% Executes on text change in EditZoomWindow

global zoomWindowMs                                     % used, then updated
global isIncrDecr                                       % used, then updated
global nSamples siMs                                    % used
global editEventNo                                      % used

% Save previous zoomWindowMs
zoomWindowMsPrev = zoomWindowMs;

% Update window size for examining events
zoomWindowMs = str2double(get(hObject, 'String'));

% If zoom window size not a number or not within range, display error
%   and choose a valid zoom window size
if isnan(zoomWindowMs) || ...
    zoomWindowMs <= 0 || zoomWindowMs > siMs * nSamples
    % Display error dialog box
    errormsg = sprintf(['Window size for examining events ', ...
                          'must be a positive number ', ...
                          'less than %d!'], siMs * nSamples);
    dlgname = 'Zoom Window Size Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
                 
    % Choose a valid zoom window size
    if isnan(zoomWindowMs)                  % if zoom window size not a number
        % Use previous zoomWindowMs
        zoomWindowMs = zoomWindowMsPrev;
    elseif zoomWindowMs <= 0                % if zoom window size too small
        % Use 1 ms
        zoomWindowMs = 1;
    elseif zoomWindowMs > siMs * nSamples   % if zoom window size too large
        % Use siMs * nSamples
        zoomWindowMs = siMs * nSamples;
    end

    % Update GUI
    set(hObject, 'String', num2str(zoomWindowMs));   
end

% Update xlimits of sweep figure
zoom_event;

% Give EditEventNo focus when not incrementing/decrementing
if isIncrDecr
    % Give focus to hObject (editZoomWindow)
    uicontrol(hObject);

    % Restore isIncrDecr state
    isIncrDecr = false;
else
    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(editEventNo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckboxChecked_Callback(hObject, ~)
%% Executes on checking or unchecking CheckboxChecked
%   Note: Assumes that checkbox is disabled when eventNoInAll is not in range

global eventNoInAll                                     % used
global editEventNo                                      % used

% Get checkbox value
valCheckboxChecked = get(hObject, 'Value');             % checkbox value

% Update isChecked status for this event
switch valCheckboxChecked
case get(hObject, 'Max')
    % Update isChecked status & fill circle of event peak in sweep figure
    mark_event_checked(eventNoInAll);
case get(hObject, 'Min')
    % Update isChecked status & unfill circle of event peak in sweep figure
    mark_event_unchecked(eventNoInAll);
otherwise
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonClean_Callback(hObject, ~)
%% Executes on toggle of ButtonClean

global abort                                            % updated
global allowRemoveClick                                 % updated
global allowAddClick                                     % used
global sweepFigure                                      % used
global editEventNo                                      % used
global colorButton colorButtonOn                        % used

% Get toggle state of add event button and act accordingly
togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
switch togstate
case togmax                 % if user toggles ON
    % Turn on mouse click removal for sweepFigure
    allowRemoveClick = true;

    % Mouse click add should not be on for sweepFigure
    %   since ButtonClean is disabled when ButtonAdd is clicked
    if allowAddClick
        error('error with code!');
    end

    % Change text on the button to 'Cleaning' and change the color to red
    set(hObject, 'String', 'Cleaning');                 % change button text
    set(hObject, 'BackgroundColor', rgb(colorButtonOn));% change button color
case togmin                 % if user toggles OFF
    % Turn off mouse click removal for sweepFigure
    allowRemoveClick = false;

    % Reset state, text and color of button
    set(hObject, 'String', 'Clean');                    % reset button text
    set(hObject, 'BackgroundColor', rgb(colorButton));  % reset button color
otherwise
    error('togstate unrecognised!\n\n');
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonAdd_Callback(hObject, ~)
%% Executes on toggle of ButtonAdd

global abort                                            % updated and used
global nEventsOld eventNoInAllOld eventNoInAllLastOld   % updated and used
global eventInfoOld eventClassOld isCheckedOld          % updated and used
global bOrPFactor                                       % updated
global allowAddClick                                    % updated
global eventInfo eventClass isChecked                   % used, then updated
global eventNoInAll eventNoInAllLast                    % used, then updated
global nEvents                                          % used, then updated
global allowRemoveClick                                 % used
global swpT swpD eOrIFactor                             % used
global colorButton colorButtonOn                        % used
global sweepFigure                                      % used
global allUiControl                                     % used
global editClassNo editEventNo                          % used
global editCorrRange                                    % used
global radioButtons2 buttonGroup2 buttonDone            % used
global dispNEvents                                      % used
global corrRangeMs classNoNewUser                       % used
global buttonGroup1                                     % used
global buttonClean                                      % used
global idxNew                                           % used
global hPlot                                            % used
global markerSize markerEdgeWidth                       % used
global checkboxChecked                                  % used

% Initialize a variable that is only true if an event is added successfully
success = false;

% Don't allow user to remove click at the same time
if allowRemoveClick
    % Turn off mouse click removal for sweepFigure
    allowRemoveClick = false;

    % Reset state, text and color of button
    set(buttonClean, 'String', 'Clean');                    % reset button text
    set(buttonClean, 'BackgroundColor', rgb(colorButton));  % reset button color
end

% Do nothing if not prepared
if isnan(corrRangeMs) || isnan(classNoNewUser)
    % Create error message
    if isnan(corrRangeMs)       % if no cursor correction range (ms) specified
        errormsg = 'Cursor correction range not specified!';
    elseif isnan(classNoNewUser)    % if no class to change specified
        errormsg = 'Class number to change to not specified!';
    end

    % Display error dialog box and return
    dlgname = 'Add Event Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
    return;
end

% Find all other buttons, sliders, checkboxes & editable text fields
%   except this button, the buttons for exiting the program, and
%   the cursor correction range
tempUiControl = setdiff(allUiControl, hObject);
tempUiControl = setdiff(tempUiControl, buttonDone);
tempUiControl = setdiff(tempUiControl, buttonGroup2);
tempUiControl = setdiff(tempUiControl, radioButtons2);
pLotherchildren = setdiff(tempUiControl, editCorrRange);

% Get toggle state of add event button and act accordingly
togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
switch togstate
case togmax                 % if user toggles ON
    % Store enable/disable UIControls in UserData
    % TODO: the following is only necessary if there are other enable/disables
    %       in the program
    % for k = 1:length(pLotherchildren)
    %     set(pLotherchildren(k), 'UserData', get(pLotherchildren(k), 'Enable'));
    % end

    % Store original values of variables to be updated
    nEventsOld = nEvents;
    eventNoInAllOld = eventNoInAll;
    eventNoInAllLastOld = eventNoInAllLast;
    eventInfoOld = eventInfo;
    eventClassOld = eventClass;
    isCheckedOld = isChecked;

    % Change text on the button to 'Adding' and change the color to red
    set(hObject, 'String', 'Adding');                   % change button text
    set(hObject, 'BackgroundColor', rgb(colorButtonOn));% change button color

    % Disable all other UIControls   
    set(pLotherchildren, 'Enable', 'off');    

    % Give SweepFigure focus
    axes(sweepFigure);

    % Attempt to add an event
    while ~success
        % Open a help dialogue box
        helpmsg = ['Please mark the time of the breakpoint', ...
                     ' as accurately as possible'];
        dlgname = 'Breakpoint Selection';
        uiwait(helpdlg(helpmsg, dlgname));

        % Turn off zoom or pan modes
        zoom(gcbf, 'off');
        pan(gcbf, 'off');

        % Set the breakpoint or peak factor to -1 (breakpoint)
        bOrPFactor = -1;

        % Allow sweep figure to respond to mouse left clicks for adding
        allowAddClick = true;

        % Wait for user input on breakpoint time
        uiwait(gcf);

        % If the next input occurs after the user aborts, exit the function
        if abort
            % Reset abort
            abort = false;
            return
        end

        % Obtain result of SweepFigure_ButtonDownFcn()
        idxBreak = idxNew;

        % Open a help dialogue box
        helpmsg = ['Please mark the time of the peak', ...
                     ' as accurately as possible'];
        dlgname = 'Peak Selection';
        uiwait(helpdlg(helpmsg, dlgname));

        % Turn off zoom or pan modes
        zoom(gcbf, 'off');
        pan(gcbf, 'off');

        % Set the breakpoint or peak factor to 1 (peak)
        bOrPFactor = 1;

        % Allow sweep figure to respond to mouse left clicks for adding
        allowAddClick = true;

        % Wait for user input on breakpoint time
        uiwait(gcf);

        % If the next input occurs after the user aborts, exit the function
        if abort
            % Reset abort
            abort = false;
            return
        end
   
        % Obtain result of SweepFigure_ButtonDownFcn()
        idxPeak = idxNew;

        % Locate the previous event
        eventNoPrev = find(eventInfoOld(:, 2) <= idxPeak, 1, 'last');
        if isempty(eventNoPrev)
            eventNoPrev = 0;
        end

        % Get the peak index of the previous event
        idxPeakPrev = eventInfoOld(eventNoPrev, 2);

        % Locate the next event
        eventNoNext = find(eventInfoOld(:, 2) >= idxPeak, 1, 'first');
        if isempty(eventNoNext)
            eventNoPrev = nEventsOld + 1;
        end

        % Get the breakpoint index of the next event
        idxBreakNext = eventInfoOld(eventNoNext, 1);

        % Check breakpoint and peak index relative to
        %   the previous and the next event
        if idxBreak >= idxPeak || ...
            idxBreak <= idxPeakPrev || idxPeak >= idxBreakNext
            % Event addition was not successful
            success = false;
            
            % Decide on error message
            if idxBreak >= idxPeak
                errormsg = 'Peak time cannot precede breakpoint time!';
            elseif idxBreak <= idxPeakPrev
                errormsg = ['The new breakpoint cannot be ', ...
                            'before the previous peak!'];
            elseif idxPeak >= idxBreakNext
                errormsg = ['The new peak cannot be ', ...
                            'after the next breakpoint!'];
            end

            % Display error dialogue box
            dlgname = 'Add Event Error';
            uiwait(errordlg(errormsg, dlgname, 'modal'));
        else
            % Event addition was successful
            success = true;
        end

        % Remove cursor marks if not successful
        if ~success
            if isfield(hPlot, 'cursorMark')
                for k = 1:numel(hPlot.cursorMark)
                    delete(hPlot.cursorMark{k});
                end
            end
        end
    end

    % This is the point of no return! Set Interruptible to be 'off'
    set(hObject, 'Interruptible', 'off');

    % First change class number to 'all' and call EditClassNo_Callback
    %   Note: this should only update thisClassInfo, rankThisClass & nThisClass
    set(editClassNo, 'String', 'all');
    EditClassNo_Callback(editClassNo)

    % Initialize this event as a removed event (Class 8)
    thisInfo = zeros(1, 11);
    thisInfo(1) = idxBreak;                         % breakpoint index
    thisInfo(2) = idxPeak;                          % peak index
    thisInfo(3) = swpD(idxBreak);                   % breakpoint value
    thisInfo(4) = swpD(idxPeak);                    % peak value
    thisInfo(5) = abs(thisInfo(4) - thisInfo(3));   % peak amplitude
    thisInfo(6) = idxPeak - idxBreak;               % 0-100% rise time (samples)
    rise10Val = thisInfo(3) + eOrIFactor * 0.1 * thisInfo(5);
                                                    % 10% of event peak value
    rise90Val = thisInfo(3) + eOrIFactor * 0.9 * thisInfo(5);   
                                                    % 90% of event peak value
    rise10Idx = find_custom(swpD(idxBreak:idxPeak) * eOrIFactor < ...
                                 rise10Val * eOrIFactor, 1, 'last', ...
                                 'ReturnNan', true);
    rise90Idx = find_custom(swpD(idxBreak:idxPeak) * eOrIFactor > ...
                                 rise90Val * eOrIFactor, 1, 'first', ...
                                 'ReturnNan', true);
    thisInfo(7) = rise90Idx - rise10Idx;            % 10-90% rise time (samples)
    thisInfo(8) = NaN;                          % inter-event interval (samples)
    thisInfo(9) = NaN;                      % inter-stimulus interval (samples)
    thisInfo(10) = NaN;                             % 50% decay time (samples)
    thisInfo(11) = NaN;                             % full decay time (samples)

    % Find the event number to add
    eventNoToAdd = eventNoPrev + 1;             % event number to add

    % Insert new event into eventInfo in the correct order
    %   Insert 8 (removed) in the corresponding place for eventClass
    %   Insert false in the corresponding place for isChecked
    %   Note: Mark event removed before deciding on actual class number
    if eventNoPrev == 0                     % if there's no previous event
        % Insert at the beginning
        eventInfo = [thisInfo; eventInfoOld];
        eventClass = [8; eventClassOld];
        isChecked = [false; isCheckedOld];
    elseif eventNoNext == nEventsOld + 1    % if there's no next event
        % Insert at the end
        eventInfo = [eventInfoOld; thisInfo];
        eventClass = [eventClassOld; 8];
        isChecked = [isCheckedOld; false];
    else                                    % if previous and next events exist
        % Insert in between previous and next events
        eventInfo = [eventInfoOld(1:eventNoPrev, :); ...
                     thisInfo; eventInfoOld(eventNoNext:end, :)];
        eventClass = [eventClassOld(1:eventNoPrev); ...
                     8; eventClassOld(eventNoNext:end)];
        isChecked = [isCheckedOld(1:eventNoPrev); ...
                     false; isCheckedOld(eventNoNext:end)];
    end

    % Update nEvents
    nEvents = nEventsOld + 1;                   % update total number of events
    set(dispNEvents, 'String', nEvents);        % update display of nEvents

    % If the currently examined event is after this new event, 
    %   update eventNoInAll to new value
    if eventNoInAllOld >= eventNoToAdd
        eventNoInAll = eventNoInAllOld + 1;
    end

    % If the last examined event is after this new event, 
    %   update eventNoInAll to new value
    if eventNoInAllLastOld >= eventNoToAdd
        eventNoInAllLast = eventNoInAllLastOld + 1;
    end

    % Update eventRank
    ButtonGroup1_SelectionChangedFcn(buttonGroup1);

    % Execute callback function for EditClassNo 
    %   (this updates thisClassInfo, rankThisClass & nThisClass)
    EditClassNo_Callback(editClassNo);

    % Remove previous plot of event breakpoints
    if ~isa(hPlot.eventBreaks, 'matlab.graphics.GraphicsPlaceholder')
        % Delete plot
        delete(hPlot.eventBreaks);
    
        % Reinitialize as GraphicsPlaceholder
        hPlot.eventBreaks = gobjects;
    end
    
    % Replot all event breakpoints as black crosses
    h = plot(swpT(eventInfo(:, 1)), eventInfo(:, 3), 'kx', ...
             'MarkerSize', markerSize, 'LineWidth', markerEdgeWidth, ...
             'DisplayName', 'Breakpoints', 'HitTest', 'off');
    if ~isempty(h)
        hPlot.eventBreaks = h;
    end

    % Zoom in on new event
    %   Note: this should update eventNoInAll to be this event
    set(editEventNo, 'String', num2str(eventNoToAdd));
    EditEventNo_Callback(editEventNo);
    if eventNoInAll ~= eventNoToAdd
        error('error with code!');
    end

    % Change class number for this event to classNoNew
    %   Note: Must follow change_class() by update_peak_markers()
    [classNoNew, prevClassNoOld, prevClassNoAuto, ...
        nextClassNoOld, nextClassNoAuto] = ...
        change_class(eventNoToAdd, classNoNewUser, 8);

    % Put all the changed classes together 
    %   and update the markers for those that are changed
    classesChanged = [8, classNoNew, ...
                         prevClassNoOld, prevClassNoAuto, ...
                         nextClassNoOld, nextClassNoAuto];
    update_peak_markers(classesChanged);

    % Set displayed class number to new value 
    % TODO: When to do this?
    % set(editClassNo, 'String', num2str(classNoNew));

    % Execute callback function for EditClassNo 
    %   (this updates thisClassInfo, rankThisClass & nThisClass and 
    %       update xlimits, checkbox, annotation if event is changed)
    EditClassNo_Callback(editClassNo);

    % Mark event checked
    mark_event_checked(eventNoToAdd, checkboxChecked);

    % Toggle button off
    set(hObject, 'Value', togmin);                      % toggle OFF
case togmin                 % if user toggles OFF
    % Prompt user on whether to cancel
    qString = 'Event addition incomplete, are you sure you want to cancel?';
    qTitle = 'Confirm cancellation';
    choice1 = 'YES, I give up!';
    choice2 = 'NO, I''m stilling trying ...';
    answer = questdlg(qString, qTitle, choice1, choice2, choice1);
    switch answer
    case choice1
        % Set abort to true to talk to the interrupted Callback
        abort = true;

        % Undo everything that was partially completed
        nEvents = nEventsOld;
        eventNoInAll = eventNoInAllOld;
        eventNoInAllLast = eventNoInAllLastOld;
        eventInfo = eventInfoOld;
        eventClass = eventClassOld;
        isChecked = isCheckedOld;

        % Execute callback function for EditClassNo 
        %   (this updates thisClassInfo, rankThisClass & nThisClass)
        EditClassNo_Callback(editClassNo);

        % Zoom in on last examined event
        set(editEventNo, 'String', num2str(eventNoInAll));
        EditEventNo_Callback(editEventNo);
    case choice2
        % Exit without resetting things
        return;
    otherwise
        error('Problem with code!');
    end
otherwise
    error('togstate unrecognised!\n\n');
end

%% Do the following if either the user successfully adds an event or 
%   the user cancels
if success || abort
    % Reset allowAddClick
    allowAddClick = false;

    % Remove cursor marks
    if isfield(hPlot, 'cursorMark')
        for k = 1:numel(hPlot.cursorMark)
            delete(hPlot.cursorMark{k});
        end
    end

    % Re-enable previously disabled UIControls
    % TODO: the following is only necessary if there are other enable/disables
    %       in the program
    % for k = 1:length(pLotherchildren)
    %     set(pLotherchildren(k), 'Enable', get(pLotherchildren(k), 'UserData'));    
    % end
    set(pLotherchildren, 'Enable', 'on');

    % Reset state, text and color of button
    set(hObject, 'String', 'Add Event');                    % reset button text
    set(hObject, 'BackgroundColor', rgb(colorButton));      % reset button color

    % Reset the Interruptible property
    set(hObject, 'Interruptible', 'on');

    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(editEventNo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditCorrRange_KeyPressFcn(hObject, callbackData)
%% Executes on key press in EditCorrRange

global isIncrDecr                                       % updated
global zoomWindowMs                                     % used

% If up or down arrow is pressed, increment or decrement, respectively. 
zWMs = zoomWindowMs;                    % window for examining events (ms)
switch callbackData.Key
case 'downarrow'
    increment_editbox(hObject, zWMs/100, zWMs, -zWMs/100, {});
case 'uparrow'
    increment_editbox(hObject, zWMs/100, zWMs, zWMs/100, {});
otherwise
    % Otherwise, do nothing
    return;
end

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditCorrRange
EditCorrRange_Callback(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditCorrRange_Callback(hObject, ~)
%% Executes on text change in EditCorrRange

global corrRangeMs                                      % updated
global isIncrDecr                                       % used, then updated
global zoomWindowMs                                     % used
global editEventNo                                      % used

% Save previous cursor correction range (ms)
corrRangeMsPrev = corrRangeMs;

% Update cursor correction range (ms)
corrRangeMs = str2double(get(hObject, 'String'));

% If cursor correction range not a number or not smaller than zoom window size, 
%   display error and choose a valid cursor correction range
if isnan(corrRangeMs) || corrRangeMs <= 0 || corrRangeMs > zoomWindowMs
    % Display error dialog box
    errormsg = sprintf(['Cursor correction range for adding events ', ...
                          'must be a positive number ', ...
                          'less than %g!'], zoomWindowMs);
    dlgname = 'Cursor Correction Range Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));

    % Choose a valid cursor correction range
    if isnan(corrRangeMs)                   % if correction range not a number
        % Use previous corrRangeMs
        corrRangeMs = corrRangeMsPrev;
    elseif corrRangeMs <= 0                 % if correction range too small
        % Use zoomWindowMs/100
        corrRangeMs = zoomWindowMs/100;
    elseif corrRangeMs > zoomWindowMs       % if correction range too large
        % Use zoomWindowMs
        corrRangeMs = zoomWindowMs;
    end

    % Update GUI
    set(hObject, 'String', num2str(corrRangeMs));   
end

% Give EditEventNo focus when not incrementing/decrementing
if isIncrDecr
    % Give focus to hObject (editCorrRange)
    uicontrol(hObject);

    % Restore isIncrDecr state
    isIncrDecr = false;
else
    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(editEventNo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditAddClass_KeyPressFcn(hObject, callbackData)
%% Executes on key press in EditCorrRange

global isIncrDecr                                       % updated
global nClass                                           % used

% If up or down arrow is pressed, increment or decrement, respectively. 
switch callbackData.Key
case 'downarrow'
    increment_editbox(hObject, 1, nClass, -1, {''});
case 'uparrow'
    increment_editbox(hObject, 1, nClass, 1, {''});
otherwise
    % Otherwise, do nothing
    return;
end

% Update isIncrDecr state
isIncrDecr = true;

% Execute callback function for EditAddClass
EditAddClass_Callback(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditAddClass_Callback(hObject, ~)
%% Executes on text change in EditAddClass

global classNoNewUser                                   % updated
global isIncrDecr                                       % used, then updated
global nClass                                           % used
global editEventNo                                      % used

% Update event class number to add or change
classNoNewUser = str2double(get(hObject, 'String'));

% If not within range, display error and change to valid value or ''
if isnan(classNoNewUser) || mod(classNoNewUser, 1) ~= 0 || ...
    classNoNewUser < 1 || classNoNewUser > nClass
    % Display error dialog box
    errormsg = sprintf(['New class number must be ', ...
                          'a positive integer between 1 and %d!'], ...
                          nClass);
    dlgname = 'Class Number Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));

    % Set classNoNewUser according to user input
    if isnan(classNoNewUser)                        % if not a number
        % Set classNoNewStr to ''
        classNoNewStr = '';
    else
        if mod(classNoNewUser, 1) ~= 0              % if not an integer
            % Find closest integer
            closestInteger = round(classNoNewUser);
            
            % Force classNoNewUser within range
            classNoNewUser = min(nClass, max(1, closestInteger));            
        elseif classNoNewUser < 1                   % if too small
            % Set classNoNewUser to minimum value
            classNoNewUser = 1;
        elseif classNoNewUser > nClass              % if too large
            % Set classNoNewUser to maximum value
            classNoNewUser = nClass;
        end
        % Update classNoNewStr
        classNoNewStr = num2str(classNoNewUser);
    end
    
    % Update GUI
    set(hObject, 'String', classNoNewStr);
end

% Give EditEventNo focus when not incrementing/decrementing
if isIncrDecr
    % Give focus to hObject (editAddClass)
    uicontrol(hObject);

    % Restore isIncrDecr state
    isIncrDecr = false;
else
    % Give EditEventNo focus (this highlights the text in the edit box)
    uicontrol(editEventNo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonRemove_Callback(~, ~)
%% Executes on press of ButtonRemove

global eventNoInAll                                     % used

% Do nothing if no event selected and doForAll not true
if isnan(eventNoInAll) && ~doForAll
    % Display error dialog box and return
    errormsg = 'No event selected!';
    dlgname = 'Remove Event Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
    return;
end

% Remove currently zoomed event
remove_event(eventNoInAll);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonChange_Callback(~, ~)
%% Executes on press of ButtonChange (changes class of an event)

global eventClass                                       % used
global eventNoInAll                                     % used
global classNoNewUser                                   % used
global editClassNo editEventNo                          % used
global doForAll toPrompt                                % used
global eventInfo                                        % used
global checkboxChecked checkboxNoPrompt                 % used

% Do nothing if not prepared or if no event selected and doForAll not true
if isnan(classNoNewUser) || (isnan(eventNoInAll) && ~doForAll)
    % Create error message
    if isnan(eventNoInAll) && ~doForAll % if no event selected
        errormsg = 'No event selected!';
    elseif isnan(classNoNewUser)        % if no class to change specified
        errormsg = 'Class number to change to not specified!';
    end

    % Display error dialog box and return
    dlgname = 'Class Change Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
    return;
end

% Retrieve the old class number of the chosen event
classNoOld = eventClass(eventNoInAll);

% Compare old class number with class number to change to
if classNoNewUser == classNoOld
    % Display error dialog box and return
    errormsg = 'New class number is same as old number!';
    dlgname = 'Class Change Error';
    uiwait(errordlg(errormsg, dlgname, 'modal'));
    return;
end

%{
% Display question dialogue box: 
answer = questdlg('Are you sure about changing the class of this event?', ...
                  'Prompt for Class Change', 'Yes');
switch answer
case 'Yes'
    % Proceed
case {'No', 'Cancel'}
    % Do nothing
    return;
otherwise
    error('Problem with code!');
end
%}

% Change class number for the chosen event or all events of same class
if doForAll
    % Do not prompt: %TODO: remove this when prompt is implemented with rezoomming
    set(checkboxNoPrompt, 'Value', get(checkboxNoPrompt, 'Max'));
    CheckboxNoPrompt_Callback(checkboxNoPrompt);

    % Get the event numbers of all events of this class
    eventNosToChange = find(eventClass == classNoOld);

    % Get the number of events to change
    nEventsToChange = length(eventNosToChange);

    % Change the class number for all events of this class
    classesChanged = [];
    for iEventsToChange = 1:nEventsToChange
        % Retrieve the number of this event in all events
        eventNoThis = eventNosToChange(iEventsToChange);

        % Change the class of this event and return the new class number
        %   Note: Must follow change_class() by update_peak_markers()
        [classNoNewThis, prevClassNoOld, prevClassNoAuto, ...
            nextClassNoOld, nextClassNoAuto] = ...
                change_class(eventNoThis, classNoNewUser, classNoOld);

        % Add any changed class to classesChanged
        classesChanged = [classesChanged, classNoOld, classNoNewThis, ...
                             prevClassNoOld, prevClassNoAuto, ...
                             nextClassNoOld, nextClassNoAuto];

        % TODO: Mark event checked?
        % mark_event_checked(eventNoThis, checkboxChecked);

        % Extract the classNoNew for the chosen event
        if eventNoThis == eventNoInAll
            classNoNew = classNoNewThis;
        end
    end

    % Update the markers for those classes that are changed
    update_peak_markers(classesChanged);
else
    % Change class number for this event to classNoNewUser,
    %   and return class number it was actually changed to (classNoNew)
    %   Note: Must follow change_class() by update_peak_markers()
    [classNoNew, prevClassNoOld, prevClassNoAuto, ...
        nextClassNoOld, nextClassNoAuto] = ...
        change_class(eventNoInAll, classNoNewUser, classNoOld);

    % Put all the changed classes together 
    %   and update the markers for those that are changed
    classesChanged = [classNoOld, classNoNew, ...
                         prevClassNoOld, prevClassNoAuto, ...
                         nextClassNoOld, nextClassNoAuto];
    update_peak_markers(classesChanged);
end

% If an event is examined and class change successful, 
%   mark event checked and change class number
if ~isnan(eventNoInAll)
    if classNoNew ~= classNoOld
        % Mark event checked
        mark_event_checked(eventNoInAll, checkboxChecked);

        % Set displayed class number to new class number
        % TODO: When to do this?
        % set(editClassNo, 'String', num2str(classNoNew));
    elseif ~toPrompt            % if the user wasn't prompt
        % Show a message box that tells user that the class couldn't be changed
        %   This message box don't need to be deleted before continuing
        message = 'Yuck! The class of the event wasn''t changed!';
        mTitle = 'Class not changed';
        % TODO: implement icon
        uiwait(msgbox(message, mTitle, 'non-modal'));
    else
        % If the user was already prompted, 
        %   no need to notify that the event class couldn't be changed
    end
end

% Execute callback function for EditClassNo 
%   (this updates thisClassInfo, rankThisClass & nThisClass and 
%       update xlimits, checkbox, annotation if event is changed)
EditClassNo_Callback(editClassNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckboxForAll_Callback(hObject, ~)
%% Executes on checking or unchecking CheckboxForAll

global doForAll                                         % updated
global editEventNo                                      % used

% Get checkbox value
valCheckboxForAll = get(hObject, 'Value');              % checkbox value

% Update doForAll status for this event
switch valCheckboxForAll
case get(hObject, 'Max')
    doForAll = true;
case get(hObject, 'Min')
    doForAll = false;
otherwise
    error('Error with code logic!');
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckboxNoPrompt_Callback(hObject, ~)
%% Executes on checking or unchecking CheckboxNoPrompt

global toPrompt                                           % updated
global editEventNo                                        % used

% Get checkbox value
valCheckboxNoPrompt = get(hObject, 'Value');              % checkbox value

% Update toPrompt status for this event
switch valCheckboxNoPrompt
case get(hObject, 'Max')
    toPrompt = false;
case get(hObject, 'Min')
    toPrompt = true;
otherwise
    error('Error with code logic!');
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckboxByForce_Callback(hObject, ~)
%% Executes on checking or unchecking CheckboxByForce

global byForce                                          % updated
global editEventNo                                      % updated

% Get checkbox value
valCheckboxByForce = get(hObject, 'Value');             % checkbox value

% Update byForce status for this event
switch valCheckboxForAll
case get(hObject, 'Max')
    byForce = true;
case get(hObject, 'Min')
    byForce = false;
otherwise
    error('Error with code logic!');
end

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SweepFigure_ButtonDownFcn(hObject, ~)
%% Executes upon mouse click in sweep figure

global idxNew                                           % updated
global allowAddClick                                    % used, then updated
global allowRemoveClick                                 % used
global abort                                            % used
global hPlot                                            % used
global swpT swpD eOrIFactor                             % used
global markerSize markerEdgeWidth                       % used
global bOrPFactor                                       % used
global corrRangeMs                                      % used
global eventInfo                                        % used
global siMs                                             % used

% If user is not ready to select a breakpoint or a peak for adding
%   or to remove an event, do nothing
if ~allowAddClick && ~allowRemoveClick
    return;
end

if allowAddClick || allowRemoveClick
    % The current dialog name for possible errors
    dlgname = 'Cursor Correction Error';

    % Do nothing if there is no cursor correction range defined
    %% TODO: Do not allow add or removal if not zoomed?
    if isnan(corrRangeMs)     % if no cursor correction range (ms) specified
        errormsg = 'Cursor correction range not specified!';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    else
        % Convert the cursor correction range (ms) to samples 
        corrRangeSamples = round(corrRangeMs / siMs);
    end

    % Get the cursor time position
    pos = get(hObject, 'CurrentPoint');
    timeAppr = pos(1, 1);
    if isempty(timeAppr)
        errormsg = 'Cursor time position cannot be obtained!';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    end

    % Get the index of the time vector corresponding to selected time
    idxAppr = find(swpT >= timeAppr, 1, 'first');
    if isempty(idxAppr)
        errormsg = 'Cannot find the nearest index from mouse click!';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    end
end

if allowAddClick
    % Get the direction factor
    %   Note: eOrIFactor is 1 for IPSC and -1 for EPSC
    %         bOrPFactor is 1 for peak and -1 for breakpoint
    %         directionFactor is 1 for finding maximum and -1 for finding minimum
    directionFactor = eOrIFactor * bOrPFactor;

    % Use adjust_peaks to find the nearest breakpoint/peak within corrRangeSamples
    [idxNew, ~] = adjust_peaks(swpD, idxAppr, directionFactor, corrRangeSamples);

    % Plot a marker on the found breakpoint/peak
    hPlot.cursorMark{2.5 + bOrPFactor * 1.5} = ...
        plot(swpT(idxNew), swpD(idxNew), 'x', 'MarkerSize', markerSize, ...
                                         'LineWidth', markerEdgeWidth, ...
                                         'Color', rgb('RoyalBlue'), ...
                                         'HitTest', 'off');
    hPlot.cursorMark{2.5 + bOrPFactor * 0.5} = ...
        plot(swpT(idxNew), swpD(idxNew), 'o', 'MarkerSize', markerSize, ...
                                         'LineWidth', markerEdgeWidth, ...
                                         'Color', rgb('RoyalBlue'), ...
                                         'HitTest', 'off');
    drawnow;

    % Disallow sweepFigure to respond to left mouse clicks
    allowAddClick = false;

    % Resume ButtonAdd_Callback()
    uiresume(gcbf);
end

if allowRemoveClick
    % Retrieve the event peak indices
    indEventPeaks = eventInfo(:, 2);
    
    % Find the first and last event numbers of all events with peak indices 
    %   within range (idxAppr - corrRangeSamples to idxAppr + corrRangeSamples)
    eventNoFirst = find(indEventPeaks > idxAppr - corrRangeSamples, 1, 'first');
    eventNoLast = find(indEventPeaks < idxAppr + corrRangeSamples, 1, 'last');

    % Don't remove anything if there are more than one events defined
    if eventNoFirst == eventNoLast
        % Choose this event
        eventNoToRemove = eventNoFirst;
    elseif eventNoFirst < eventNoLast
        % Create an error message
        % TODO: Don't show this when toPrompt == false?
        errormsg = {['There are more than one event within ', ...
                    'cursor correction range!'], ...
                    'No events removed!'};

        % Display an error dialog box and return
        dlgname = 'Remove Event Error';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    else
        % Create an error message
        % TODO: Don't show this when toPrompt == false?
        errormsg = {['There are no events found in the ', ...
                    'cursor correction range!'], ...        
                    'No events removed!'};

        % Display an error dialog box and return
        dlgname = 'Remove Event Error';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    end

    % Remove that event
    remove_event(eventNoToRemove);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ButtonDone_Callback(~, ~, buttonGroup2)
% Upon button press of DONE!

global finalCommand                                     % updated
global gui                                              % used
global abort                                            % used

% Store radiobutton selection
finalCommand = buttonGroup2.SelectedObject.Tag;

% Prompt user if abort not resolved
% TODO
if abort
    %TODO
end

% Close GUI
close(gui);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the uibuttongroup visible after creating child objects. 
% bg.Visible = 'on';

%% horizontal scroll bar.  Mark stole this from the web.
% figure('scroll', 'on')

%{
OLD CODE:

skip = 0;
if skip == 0
end

%% button group on the left. Mark stole this from the web.

plot(tVec, dataRaw, 'k'); 
plot(tVec, dataLowpass, 'b');
plot(eventInfo(:, 1), eventInfo(:, 3), 'x', ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 4);
plot(eventInfo(:, 2), eventInfo(:, 4), 'o', ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 4);
% TODO: Fix title
title(strcat(dataSubdirectory, ', ', {' '}, 'Trace/File', {' '}, ...
        num2str(cellNum)), 'interpreter', 'none', ...
        'FontSize', 20, 'FontWeight', 'bold');

posSweepFig = [0.13, 0.11, 0.775, 0.815];   % this is the Matlab default
%axes('Parent', pR, 'Position', posSweepFig);

%textFontUnitsDefault = 'points';
%textFontSizeDefault = 10;

set(gcf, 'WindowButtonDownFcn', 'xPOINTS = pickPointsMinis(gca);')

if exist('grayLine', 'var')
    delete(grayLine);
end
if exist('replot', 'var')
    delete(replot);
end

set(buttonRemove, 'FontSize', textFontSizeButton/2);
set(buttonChange, 'FontSize', textFontSizeButton/2);
'String', 'Remove Event', ...
'String', 'Change Class', ...

eventNoLast = str2double(get(editEventNo, 'String'));

% Save "previous class": the class that was last examined
lastClassInfo = thisClassInfo;
nLastClass = nThisClass;
classNoStrLast = classNoStr;

            % Update thisClassInfo, nThisClass using extract_class()
            [thisClassInfo, nThisClass] = extract_class(classNoStr);
            elseif eventNoNext ~= eventNoPrev + 1
                errormsg = 'The new peak cannot overlap an existing peak!';

% Find needed objects
gui = get(get(hObject, 'Parent'), 'Parent');              % find MainGUI
editZoomWindow = findall(gui, 'Tag', 'EditZoomWindow');   % find EditZoomWindow

% Find needed objects
gui = get(get(hObject, 'Parent'), 'Parent');              % find MainGUI
editEventNo = findall(gui, 'Tag', 'EditEventNo');         % find EditEventNo

% Find the number of this event in all events
eventNoInAll = find(eventInfo(:, 2) == idxPeak);

% Get event number of zoomed event
%   Note: cannot be 'all' because checkbox is disabled when eventNo == 'all'
eventNo = str2double(get(editEventNo, 'String'));       % number in this class

% Find needed objects
pL = get(hObject, 'Parent');                            % find panelLeft
editClassNo = findall(pL, 'Tag', 'EditClassNo');        % find EditClassNo

% Find needed objects
pL = get(hObject, 'Parent');                            % find panelLeft
editEventNo = findall(pL, 'Tag', 'EditEventNo');        % find EditEventNo
eventNoStrPrev = get(editEventNo, 'String');

%   Note: DispNThisClass will be changed as well, 
%         so need to follow by extract_class for current classNoStr

if ~(strcmpi(textEditClassNo, 'all') || ...
     isnan(classNo) || mod(classNo, 1) ~= 0 || ...
     classNo <= 1 || classNo > nClass)

if isnan(classNoPrev) || mod(classNoPrev, 1) ~= 0 || ...
    classNoPrev < 1 || classNoPrev > nClass
    return;    
end

classNoStrLast = classNoStr;
eventNoLast = NaN;
classNoStrLast = classNoStr;
eventNoLast = eventNo;

% Get lastClassInfo, nLastClass using extract_class()
%   Note: don't provide dispNThisClass to prevent from updating
[lastClassInfo, nLastClass] = extract_class(classNoStrLast);

idxPrevPeak = lastClassInfo(eventNoLast, 2);
            % Update the class number to last examined class
            set(hObject, 'String', classNoStrLast);

if ~isnan(eventNoLast) && mod(eventNoLast, 1) == 0 && ...
    eventNoLast >= 1 && eventNoLast <= nLastClass	% check if within range

if strcmpi(eventNoStr, 'all')               % to show all events
else                                        % to show a particular event

% Update thisClassInfo & nThisClass using extract_class()
[thisClassInfo, nThisClass] = extract_class(classNoStr, dispNThisClass);

function [tCInfo, nTC] = extract_class(thisClassNoStr, dispNThisClass)
%% Extract from eventInfo given thisClassNoStr

global eventInfo eventClass                             % used
global nClass                                           % used

    % Exit function
    return;
global classNoStr                                       % updated and used 
global eventNoStr                                           % updated

classNoStrPrev = ;        % previous text of EditClassNo
if strcmpi(classNoStrPrev, 'all')
end
    % If not a number, something is wrong in the code ()
    error('Code error: class number should always in range!\n')

% If previous text of EditEventNo is 'all', do nothing
eventNoStrPrev = get(hObject, 'String');        % previous text of EditEventNo
if strcmpi(eventNoStrPrev, 'all')
    return;    
end

    % Update eventNoInAll to NaN
    eventNoInAll = NaN;

%   Note: Assumes one of the following is true:
%       (1) eventNo is not a number
%       (2) eventNo is a number and thisClassInfo has at least eventNo rows
global thisClassInfo                                    % used
global eventNo                                          % used
if ~isnan(eventNo)                  % an event is selected
    idxPeak = thisClassInfo(eventNo, 2);
    idxBreak = thisClassInfo(eventNo, 1);       % index of event breakpoint
    fullDecayTime = thisClassInfo(eventNo, 11); % full decay time (samples)
       interStimulusInterval = thisClassInfo(eventNo, 9);
                              num2str(thisClassInfo(eventNo, 5), 3), ...
                              num2str(thisClassInfo(eventNo, 7) * siMs, 3), ...
                              num2str(thisClassInfo(eventNo, 10) * siMs, 3), ...
end

    % Execute callback function again
    EditClassNo_Callback(hObject, [], tVec, dataLowpass);

    % Execute callback function again;
    EditZoomWindow_Callback(hObject, [], tVec)

    % Execute callback function again
    EditEventNo_Callback(hObject, [], tVec, dataLowpass);

%   unless class number is not a number or at lower limit
if ~(isnan(classNoPrev) || classNoPrev <= 1)

% If previous class number is not a number, do nothing
if isnan(classNoPrev)
    return;    
end        % Locate the previous event
        eventNoPrev = find(eventInfoOld(:, 1) < idxBreak, 1, 'last');
        if isempty(eventNoPrev)
            eventNoPrev = 0;
        end


classNoStrPrev = get(editClassNo, 'String');    % previous class number string
classNoPrev = str2double(classNoStrPrev);       % previous class number
if strcmpi(classNoStrPrev, 'all')       % 'all' is higher than maximum
    % Set class number to maximum
    set(editClassNo, 'String', num2str(nClass));
elseif ~isnan(classNoPrev) && classNoPrev <= 1
    % Set class number string to 'all'
    set(editClassNo, 'String', 'all');
elseif ~isnan(classNoPrev)
    % Decrement class number
    set(editClassNo, 'String', num2str(classNoPrev - 1));
else
    error('Invalid class number!');
end

% Find needed objects
pL = get(hObject, 'Parent');                            % find panelLeft
editClassNo = findall(pL, 'Tag', 'EditClassNo');        % find EditClassNo

% If the key pressed is not the up or down arrow, do nothing
if ~strcmpi(callbackData.Key, 'uparrow') && ...
    ~strcmpi(callbackData.Key, 'downarrow')
    return;
end

% Get previous class number
classNoStrPrev = get(editClassNo, 'String');    % previous class number string
classNoPrev = str2double(classNoStrPrev);       % previous class number
    if classNoPrev <= 1
        return;
    else
        set(hObject, 'String', num2str(classNoPrev - 1));
    end    
    if classNoPrev >= nClass
        return;    
    else
        set(hObject, 'String', num2str(classNoPrev + 1));
    end

% If previous event number is not a number, do nothing
eventNoPrev = str2double(get(hObject, 'String'));
if isnan(eventNoPrev)
    return;    
end

% If up or down arrow is pressed and previous event number not at limits, 
%   increment or decrement the number, respectively. 
switch callbackData.Key
case 'downarrow'
    if eventNoPrev <= 1
        return;    
    else
        set(hObject, 'String', num2str(eventNoPrev - 1));
    end    
case 'uparrow'
    if eventNoPrev >= nThisClass
        return;    
    else
        set(hObject, 'String', num2str(eventNoPrev + 1));
    end
otherwise
    % This was checked before
end

% If previous zoom window size is not a number, do nothing
zoomWindowMsPrev = str2double(get(hObject, 'String'));
if isnan(zoomWindowMsPrev)
    return;    
end

% If up or down arrow is pressed and previous zoom window size not at limits, 
%   increment or decrement the number by 10, respectively. 
switch callbackData.Key
case 'downarrow'
    if zoomWindowMsPrev <= 10
        return;    
    else
        set(hObject, 'String', num2str(zoomWindowMsPrev - 10));
    end    
case 'uparrow'
    if zoomWindowMsPrev > nSamples * siMs - 10
        return;    
    else
        set(hObject, 'String', num2str(zoomWindowMsPrev + 10));
    end
otherwise
    % Do nothing
end

% Find needed objects
pL = get(hObject, 'Parent');                            % find panelLeft
editEventNo = findall(pL, 'Tag', 'EditEventNo');        % find EditEventNo

% Increment event number and execute callback function for EditEventNo
%   unless event number is not a number or at upper limit
eventNoPrev = str2double(get(editEventNo, 'String'));
if ~(isnan(eventNoPrev) || eventNoPrev >= nThisClass)
    % Increment event number
    set(editEventNo, 'String', num2str(eventNoPrev + 1));

    % Execute callback function for EditEventNo
    EditEventNo_Callback(editEventNo, [], tVec, dataLowpass);
end

% Decrement event number and execute callback function for EditEventNo
%   unless event number is not a number or at lower limit
eventNoPrev = str2double(get(editEventNo, 'String'));
if ~(isnan(eventNoPrev) || eventNoPrev <= 1)
    % Decrement event number
    set(editEventNo, 'String', num2str(eventNoPrev - 1));

    % Execute callback function for EditEventNo
    EditEventNo_Callback(editEventNo, [], tVec, dataLowpass);
end

% Increment zoom window size and execute callback function for EditZoomWindow
%   unless zoom window size is not a number or at upper limit
zoomWindowMsPrev = str2double(get(editZoomWindow, 'String'));
if ~(isnan(zoomWindowMsPrev) || zoomWindowMsPrev > nSamples * siMs - 10)
    % Increment zoom window size by 10
    set(editZoomWindow, 'String', num2str(zoomWindowMsPrev + 10));

    % Execute callback function for EditZoomWindow
    EditZoomWindow_Callback(editZoomWindow, [], tVec);
end

% Decrement zoom window size and execute callback function for EditZoomWindow
%   unless zoom window size is not a number or at lower limit
zoomWindowMsPrev = str2double(get(editZoomWindow, 'String'));
if ~(isnan(zoomWindowMsPrev) || zoomWindowMsPrev <= 10)
    % Decrement zoom window size by 10
    set(editZoomWindow, 'String', num2str(zoomWindowMsPrev - 10));

    % Execute callback function for EditZoomWindow
    EditZoomWindow_Callback(editZoomWindow, [], tVec);
end

% Find needed objects
sweepFigure = findall(gui, 'Tag', 'SweepFigure');  % find SweepFigure

% Find needed objects
gui = get(get(hObject, 'Parent'), 'Parent');            % find MainGUI
sweepFigure = findall(gui, 'Tag', 'SweepFigure');       % find SweepFigure
editEventNo = findall(gui, 'Tag', 'EditEventNo');       % find EditEventNo

% Find needed objects
gui = get(get(hObject, 'Parent'), 'Parent');            % find MainGUI
dispNThisClass = findall(gui, 'Tag', 'DispNThisClass'); % find DispNThisClass
editEventNo = findall(gui, 'Tag', 'EditEventNo');       % find EditEventNo

    %   and Type II PSCs are always followed by a Type III PSC
        eventClass(eventNoInAll - 1) == 1   % if previous event is a Type I PSC

    if classNoNew == 2 && ...               % if changed to a Type II PSC
        eventNoInAll ~= nEvents && ...      % and if not last event
        any(eventClass(eventNoInAll + 1) == [1, 2]) 
                                            % and if next event is a PSC
                                            %   but not Type III
        % The next PSC should really be Type III by truth (1)
        %   TODO
    elseif classNoNew == 3 && ...           % if changed to a Type III PSC
        eventNoInAll ~= 1 && ...            % and if not first event
        all(eventClass(eventNoInAll - 1) ~= [2, 3])   
                                            % and if previous event is 
                                            %   not a Type II or III PSC
        % The previous event should really be a Type II or III PSC by truth (2)
        %   TODO
    end
    if classNoOld == 2 && ...               % if changed from a Type II PSC
        classNoNew ~= 3 && ...              % and if not changed to a Type III PSC
        eventNoInAll ~= nEvents && ...      % if not last event
        eventClass(eventNoInAll + 1) == 3   % if next event is a Type III PSC
        if eventNoInAll ~= nEvents - 1 && ...   % if not second to last event
            eventClass(eventNoInAll + 2) == 3   % if next next event is 
                                                %   also a Type III PSC
            % The next PSC should really be Type II
            %   TODO
        else
            % The next PSC should really be Type I
            %   TODO
        end
    elseif classNoOld == 3 && ...           % if changed from a Type III PSC
        eventNoInAll ~= 1 && ...            % if not first event
        eventClass(eventNoInAll - 1) == 2   % if previous event is a Type II PSC
        % The previous PSC should really be Type I
        %   TODO        
    end

uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDone', ...
                     'String', 'DONE!', ...
                     'Position', posButtonDone, ...
                     'Callback', 'uiresume(gcbf)');    


%   Note: Initializing to gobjects(nEvents) takes a lot of time (4.3 seconds)
%           so try cell array instead. Not sure if this will be slower in the 
%           long term though
hPlot.peaksFilled = gobjects(nEvents);       % peaks that are checked and filled
    delete(hPlot.peaksFilled(eventNoInAll));
    hPlot.peaksFilled(eventNoInAll) = ...

        % Wait for user input on breakpoint time
        w = 1;                      % result for waitforbuttonpress
        while w ~= 0                % mouse button click
            w = waitforbuttonpress;
        end

% Remove cursor marks
if isfield(hPlot, 'cursorMark')
    delete(hPlot.cursorMark);
end

% Initialize zoomWindowMs
EditZoomWindow_Callback(editZoomWindow);

% Initialize corrRangeMs
EditCorrRange_Callback(editCorrRange);

% Initialize classNoNewUser
EditAddClass_Callback(editAddClass);

    hPlot.grayLine = line([swpT(idxPeak), swpT(idxPeak)], ...
                          [ylimits(1), ylimits(2)], ...
                          'Color', rgb('Gray'), 'LineStyle', '--', ...
                          'DisplayName', 'Current peak time');

global eventNoInAll                                     % used
global checkboxChecked                                  % used
% If it's the examined event, update displayed checked status
if eventNoInAll == thisEventNoInAll
    set(checkboxChecked, 'Value', get(checkboxChecked, 'Min'));
end

if ~isnan(eventNoInAll) && eventNoInAll ~= 0
    set(checkboxChecked, 'Value', get(checkboxChecked, 'Max'));
    CheckboxChecked_Callback(checkboxChecked, []);
end

    rise10Val = 0.1 * thisInfo(4);                  % 10% of event peak value
    rise90Val = 0.9 * thisInfo(4);                  % 90% of event peak value
    rise10Idx = find_custom(swpD(idxBreak:idxPeak) * eOrIFactor > ...
                            rise10Val * eOrIFactor, 1, 'first', ...
                            'ReturnNan', true) - 1;
    rise90Idx = find_custom(swpD(idxBreak:idxPeak) * eOrIFactor > ...
                            rise90Val * eOrIFactor, 1, 'first', ...
                            'ReturnNan', true) - 1;

if eventNoInAll ~= 1                    % the first event does not have prev
if eventNoInAll ~= nEvents              % the last event does not have next

    %   TODO: may need to change "previous event" to "previous PSC"
    %   TODO: may need to change to "PSC full decay time" for PSCs

classNoNew = verify_classNoNew(classNoNewUser);
function classNoNew = verify_classNoNew (classNoNewUser)
% Decide on actual class number to change to
global eventNoInAll                                     % used

% Change class number for current event to classNoNew
%   Note: assumes that eventNoInAll is a number in range
global eventNoInAll                                     % used
 
answer = questdlg(qString, qTitle, ...
                  choice1, choice2, choice3, choice1);
choice3 = 'Forget it! Don''t change!';            
case choice3
    % Do nothing
    return;
if toPrompt
else
    answer = choice1;
end

    %   Note: (1) EditClassNo_Callback will be called
    %         (2) event will be marked checked

    % Wait for the user to close the dialog box
    waitfor(errorAddEvent);
    % Wait for the user to close the dialog box
	waitfor(errorButtonAdd);
    % Wait for the user to close the dialog box
    waitfor(errorClassChange);
    % Wait for the user to close the dialog box
    waitfor(errorButtonRemove);
    % Wait for the user to close the dialog box
	waitfor(errorClassNo);
    % Wait for the user to close the dialog box
    waitfor(errorEventNo);
    % Wait for the user to close the dialog box
	waitfor(errorZoomWindow);
    % Wait for the user to close the dialog box
	waitfor(errorCorrRange);

    minLine = minD; maxLine = maxD;             % use the overall data range
    ylimits = get(gca, 'Ylim');                 % use the current y axis limits
    minLine = ylimits(1); maxLine = ylimits(2);
global minD maxD                            % used
global minD maxD                            % updated
minD = min(swpD);                           % maximum sweep data value (pA)
maxD = max(swpD);                           % minimum sweep data value (pA)

sweepFigure = ...
    axes('Parent', panelRight, 'Tag', 'SweepFigure', ...
         'ButtonDownFcn', @SweepFigure_ButtonDownFcn, ...
         'PickableParts', 'none', 'HitTest', 'off');
% Allow sweep figure to capture mouse clicks
set(sweepFigure, 'PickableParts', 'visible');
% Allow sweep figure to capture mouse clicks
set(sweepFigure, 'PickableParts', 'visible');
% Disallow sweep figure to capture mouse clicks
set(hObject, 'PickableParts', 'none');

% Allow sweep figure to respond to captured mouse clicks
set(sweepFigure, 'HitTest', 'on');
% Allow sweep figure to respond to captured mouse clicks
set(sweepFigure, 'HitTest', 'on');
% Disallow sweep figure to respond to captured mouse clicks
set(hObject, 'HitTest', 'off');

    % Save the current y axis limits
    ylimitsOld = get(gca, 'Ylim');              % use the current y axis limits
    ylimits = get(gca, 'Ylim');                 % use the current y axis limits
    minLine = ylimits(1); maxLine = ylimits(2);
    minLine = minY; maxLine = maxY;             % use the overall y axis range
    % Update y limits
    set(gca, 'YLim', ylimitsOld);

global hPlot                                            % used, then updated

% Remove previous gray dashed line if exists
if isfield(hPlot, 'grayLine')
    delete(hPlot.grayLine);
    grayLineExists = true;
else
    grayLineExists = false;
end

% Replot gray dashed line if originally present
if grayLineExists
    % Draw a new gray dashed line with the new y axis limits
    draw_gray_line(idxPeak);
end

% Use the new y axis limits for the gray dashed line
ylimits = get(gca, 'Ylim');                 % use the current y axis limits
minLine = ylimits(1); maxLine = ylimits(2);

% Change the x-axis limits and allow y axis limits to change automatically
xlim(xlimits);

% Save old ylimits
ylimitsOld = get(gca, 'Ylim');

ylim([minDataNow, maxDataNow]);

posButtonDown2Zoom = [0.10, 0.560, 0.20, 0.04];
posButtonUp2Zoom   = [0.70, 0.560, 0.20, 0.04];
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonDown2Zoom', ...
                     'String', '-100', ...
                     'Position', posButtonDown2Zoom, ...
                     'Callback', @ButtonDown2Zoom_Callback);
uicontrol(panelLeft, 'Style', 'pushbutton', 'Tag', 'ButtonUp2Zoom', ...
                     'String', '+100', ...
                     'Position', posButtonUp2Zoom, ...
                     'Callback', @ButtonUp2Zoom_Callback);
increment_editbox(editZoomWindow, 1, nSamples * siMs, -100, {});
increment_editbox(editZoomWindow, 1, nSamples * siMs, 100, {});

    % Change text on the button to Adding Event ... and change the color to red
    set(hObject, 'String', 'Adding Event ...');         % change button text

pLotherchildren = setdiff(allUiControl, hObject);
    % Convert to cursor correction range (ms) to samples 
    corrRangeSamples = round(corrRangeMs / siMs);
global bOrPFactor corrRangeSamples                      % updated
global bOrPFactor corrRangeSamples                      % used

        % Create an error message
        if isnan(corrRangeMs)  
        end

        % Display an error dialog box and return

pLotherchildren = setdiff(allUiControl, hObject, corrRangeMs);

% Give EditZoomWindow focus (this highlights the text in the edit box)
uicontrol(editZoomWindow);

% Set class number string to 'all' and execute callback function for EditClassNo
set(editClassNo, 'String', 'all');
EditClassNo_Callback(editClassNo);

% Give EditClassNo focus (this highlights the text in the edit box)
uicontrol(editClassNo);

% Give EditZoomWindow focus (this highlights the text in the edit box)
uicontrol(editZoomWindow);

% Give EditEventNo focus (this highlights the text in the edit box)
uicontrol(editEventNo);

% Give SweepFigure focus
axes(sweepFigure);

    % Update thisClassInfo
    thisClassInfo = eventInfo;

    % Update number of events in this class
    nThisClass = nEvents;                   % number of events of this class

    % Update rank of events in this class
    ButtonGroup1_SelectionChangedFcn(buttonGroup1);

    % Update dispNThisClass
    set(dispNThisClass, 'String', num2str(nThisClass));

global thisClassInfo nThisClass                         % used
    for iThisClass = 1:nThisClass
        eventNoThis = find(eventInfo(:, 2) == thisClassInfo(iThisClass, 2), 1);
        if ~isnan(eventNo) && iThisClass == eventNo

% Change class number for the examined event or all events of same class

if ~isnan(eventNoInAll)
    classNoOld = eventClass(eventNoInAll);
elseif doForAll
    classNoOld = str2double(get(editClassNo, 'String'));
    if isnan(classNoOld)
        % Display error dialog box and return
        errormsg = 'No class number to change from!';
        dlgname = 'Class Change Error';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    end
end

global thisClassInfo nThisClass                         % used
    % Get the event number in this class
    eventNoStr = get(editEventNo, 'String');
    eventNo = str2double(eventNoStr);
    %   Note: classNoNew will not be returned if isnan(eventNo) == true
    for iThisClass = 1:nThisClass
        % Retrieve the number of this event in all events
        eventNoThis = find(eventInfo(:, 2) == thisClassInfo(iThisClass, 2), 1);

        % Retrieve the number of this event in all events
        eventNoThis = find(eventInfo(:, 2) == thisClassInfo(iThisClass, 2), 1);

        % Extract the classNoNew for the examined event
        if ~isnan(eventNo) && iThisClass == eventNo
            classNoNew = classNoNewThis;
        end
    end

if ~isnan(eventNoToRemove)
    classNoOld = eventClass(eventNoToRemove);
elseif doForAll
    classNoOld = str2double(get(editClassNo, 'String'));
    if isnan(classNoOld)
        % Display error dialog box and return
        errormsg = 'No class number to change from!';
        dlgname = 'Remove Event Error';
        uiwait(errordlg(errormsg, dlgname, 'modal'));
        return;
    end
end
global thisClassInfo nThisClass                         % used
global thisClassInfo nThisClass                         % updated
global dispNThisClass                                   % used

    % First change class number to 'all' and call EditClassNo_Callback
    %   Note: this will also update eventNoInAll to be last examined event
    %         and this will use eventInfo and eventClass, so this must 
    %         be done before eventNoInAll, eventInfo & eventClass are updated

    % The currently examined event will be the last examined event after
    %   event addition
    eventNoInAllLast = eventNoInAll;

    % If last examined event is after this new event, update eventNoInAllLast
    if eventNoInAllLastOld >= eventNoInAll
        eventNoInAllLast = eventNoInAllLastOld + 1;
    end

pLotherchildren = setdiff(setdiff(allUiControl, hObject), editCorrRange);

% Find all other buttons, sliders, checkboxes & editable text fields
%   except this button and cursor correction range

        % Locate the previous event
        eventNoPrev = find(eventInfoOld(:, 1) < idxBreak, 1, 'last');
        if isempty(eventNoPrev)
            eventNoPrev = 0;
        end

        if idxBreak >= idxPeak || eventNoNext ~= eventNoPrev + 1

            elseif eventNoNext ~= eventNoPrev + 1
                errormsg = 'The new peak cannot overlap an existing peak!';

    if ~success
        return;
    end

% Reset abort
%   Note: It is necessary to put this here because abort needs time
%           to terminate any remnant ButtonAdd_Callbacks
abort = false;

% If user cancels, do nothing and resume ButtonAdd_Callback()
if abort
    uiresume(gcbf);
    return;
end

legend(hAll(validPlots));

%}

