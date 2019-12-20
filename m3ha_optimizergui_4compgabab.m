function [hfig] = m3ha_optimizergui_4compgabab (realdata_cpr, realData, hfig, varargin)
% OPTIMIZERGUI  The GUI interface for OPTIMIZER.m, which runs NEURON
% simulations and can (in auto mode) perform optimization of parameters. 
% Usage: [hfig] = m3ha_optimizergui_4compgabab (realdata_cpr, realData, hfig, varargin)
%
% Asks for data file names, loads data, generates GUI. Operates in two
% modes: "manual mode" with user-controllable sliders, and "auto mode" with
% a simplex optimization procedure. 
%
% Requires:
%       cd/m3ha_optimizer_4compgabab.m
%       /home/Matlab/Adams_Functions/set_fields_zero.m
%       /home/Matlab/Adams_Functions/restore_fields.m
%       /home/Matlab/Adams_Functions/find_in_strings.m
%       /home/Matlab/Adams_Functions/my_closereq.m
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:
%       cd/singleneuronfitting2.m and later versions
%
% By Christine K. Lee 2011-01-14    last modified 2014-04
% 2016-07-15 - Added MANUAL WITH JITTER
% 2016-07-19 - Improved uigetfile
% 2016-07-20 - Modified JITTER mode
% 2016-07-20 - Modified swprightedge (Right edge of sweep of interest)
% 2016-07-20 - Determine holding potential from real data, also used as baseline for LTS analysis
% 2016-10-05 - Added current pulse response
% 2016-10-06 - Started to reorganize code
% 2016-10-06 - Renamed figure handles so that they are now all in a structure 
%               hfig that is passed to and from functions
% 2016-10-07 - outParams.currpulse(k) is now converted to nA
% 2016-10-21 - moved files up one level to /media/adamX/m3ha/optimizer4gabab
% 2016-10-21 - added functionsdirectory, parentDirectory, infolder
% 2016-10-21 - updated locations of folders and files
% 2017-01-14 - Added numbuildparams and added build parameters to neuronparamnames, 
%               neuronparams_max, neuronparams_min, paraminit, neuronparamislog
% 2017-01-15 - Shortened cprwin_orig from [100, 500] to [100, 260] to be 
%               consistent with dclampPassiveFitter.m
% 2017-01-16 - Added rowmode so that each pharm, gincr pair is a row for rowmode == 2
% 2017-01-17 - Make cm and Ra fixed values
% 2017-01-18 - Changed 'MaxIter' & 'MaxFunEvals' to 200
% 2017-01-25 - Added outParams.numIC
% 2017-01-27 - Added outParams.plotMarkflag to plot figures for Mark's talk
% 2017-02-06 - Changed range of epas from [-120 -30] to [-90 -50]
% 2017-04-19 - Replaced across_trials & across_cells with colmode
% 2017-04-19 - Added plotwin
% 2017-04-19 - Changed positions in panels to normalized units and made things always fit 
%               without overlapping in the GUI
% 2017-04-19 - Renamed runbuttonpanel -> p6
% 2017-04-19 - Renamed fitmode -> datamode; 
%               fitmode now goes from 1~4 and mean MANUAL, AUTO, JITTER, AUTO_WITH_JITTER, resp.
% 2017-04-20 - Removed update_sliderposition() and renamed update_sliderbounds -> update_sliders()
% 2017-04-20 - Now updates neuronparams_use whenever check box is selected or non-selected
% 2017-04-21 - Added outfoldername to file names
% 2017-04-22 - Added callbacks for all check boxes and editable text fields
% 2017-04-22 - outParams.ltserrorflag is now dependent on default_ltsuse and default_ltsw
% 2017-04-24 - Removed the handles structure and cleaned code
% 2017-04-25 - Added my_closereq.m
% 2017-05-01 - Fixed across_cells attempt #4 to include all g incr 
% 2017-05-01 - Moved code to subfunction select_cellstofit()
% 2017-05-01 - Added fiti & fiti_start to iterate over all cells to fit
% 2017-05-01 - Moved code to subfunction m3ha_select_raw_traces()
% 2017-05-01 - Moved code to subfunction m3ha_import_raw_traces()
% 2017-05-01 - Created log_arraytext.m
% 2017-05-01 - Added autoparams to Panel 5
% 2017-05-12 - Move everything unrelated to GUI to singleneuronfitting.m
% 2017-05-13 - Now updates outParams.runnumtotal in m3ha_optimizer_4compgabab.m
% 2017-05-13 - Added outParams.cellname
% 2017-05-13 - Now gets parentDirectory from outParams
% 2017-05-15 - Made build params and other passive params a different background color on GUI
% 2017-05-17 - Renamed fitmode -> runmode
% 2017-05-22 - Changed line width and indentation
% 2017-05-23 - Fixed the placement of the rgb() function to be after addpath()
% 2017-05-23 - Removed modeselected from outParams and replaced with updated outParams.runMode
% 2017-05-23 - Added otherwise to all switch statements
% 2017-07-26 - Added scroll bar for NEURON parameters Panel (now Panel 8)
% 2017-08-29 - Only show scroll bar for Panel 8 if necessary
% 2018-01-24 - Added isdeployed
% 2018-08-10 - Changed swpedges to fitwin (ms) and added fitwinCpr (ms)
%% TODO: Add 'ColMode' as a parameter-value pair

global outParams

%% Parameters used for GUI construction
% General GUI appearance
doublecolumnflag = 1;           % whether two columns are needed for sliders
bgcolor0 = [0.8, 1, 0.8];       % background color of Panels 1, 2 & 6
bgcolor1 = 'Plum';              % background color of build parameters
bgcolor2 = 'LightBlue';         % background color of global passive parameters
bgcolor3 = 'PaleGreen';         % background color of pas parameters
bgcolor4 = 'LemonChiffon';      % background color of IT parameters
bgcolor5 = 'Salmon';            % background color of Ih parameters
bgcolor6 = 'Honeydew';          % background color of IKir parameters
bgcolor7 = 'Coral';             % background color of IA parameters
bgcolor8 = 'Bisque';            % background color of INaP parameters
bgcolor9 = 'Turqoise';          % background color of Iahp parameters
bgcolor10 = 'Tan';              % background color of IKCa parameters
bgcolor11 = 'Silver';           % background color of IL parameters
bgcolor12 = 'Lavender';         % background color of cad parameters
bgcolor13 = 'Khaki';            % background color of sliders
buttoncolor_OFF = [0, 1, 0];    % button color when it's not pressed
buttoncolor_ON = [1, 0, 0];     % button color when it's pressed
rowh_text = 20;                 % row height for pure texts (pixels)
rowh_swp = 23.5;                % row height for sweeps (pixels)
rowh_lts = 20;                  % row height for LTS properties (pixels)
rowh_par = 27.5;                % row height of NEURON parameters (pixels)

% Position of main GUI figure in pixels, [left bottom width height]
foutpos_singcol = [13, 52, 850, 909];
foutpos_doubcol = [13, 52, 1350, 909];

% Position of panels in normalized units, [left bottom width height]
p1pos_singcol = [0.02 0.05 0.41 0.90];      % position of Panel 1: 'Error Function Maker' 
                                            %   when doublecolumnflag == 0
p2pos_singcol = [0.45 0.05 0.52 0.90];      % position of Panel 2: 'NEURON Controller'
                                            %   when doublecolumnflag == 0
p1pos_doubcol = [0.02 0.05 0.26 0.90];      % position of Panel 1: 'Error Function Maker'
                                            %   when doublecolumnflag == 1
p2pos_doubcol = [0.30 0.05 0.67 0.90];      % position of Panel 2: 'NEURON Controller'
                                            %   when doublecolumnflag == 1
p3pos         = [0.05 0.91 0.60 0.07];      % position of Button Group 3: 'Mode Chooser' in Panel 2
p4pos         = [0.05 0.22 0.90 0.69];      % position of Panel 4: 'NEURON parameters' in Panel 2
p5pos         = [0.05 0.03 0.90 0.19];      % position of Panel 5: 'AUTO mode' in Panel 2
p6pos         = [0.70 0.91 0.30 0.07];      % position of Panel 6: 'RUN Button' in Panel 2
p7pos         = [0.00 0.09 0.98 0.91];      % position of Panel 7 in Panel 4
        
% Position of contents in normalized units, [left bottom width height]
pos_pfileload     = [0.02 0.02 0.1 0.05];   % position of 'Load P-file' button in Panel 4
pos_pfilesave     = [0.125 0.02 0.1 0.05];  % position of 'Save P-file' button in Panel 4
pos_chbounds      = [0.25 0.02 0.2 0.05];   % position of 'Change Bounds' button in Panel 4
pos_autologger    = [0 0.8 1 0.2];          % position of autologger in Panel 5
pos_autotable1    = [0 0 0.3 0.8];          % position of autotable1 in Panel 5
pos_autotable2    = [0.3 0 0.6 0.8];        % position of autotable2 in Panel 5
pos_runbutton     = [0 0 0.3 1];            % position of 'RUN' button in Panel 6
pos_twiddlebutton = [0.4 0 0.5 1];          % position of 'TWIDDLE' button in Panel 6
pos_p8slider      = [0.98 0.09 0.02 0.91];  % position of slider for Panel 8 in Panel 4

% Other positions of contents in normalized units
xpos_param        = [0, 0.17, 0.20, 0.39];  % left position of name, checkbox, slider, value 
                                            %   of NEURON parameters in Panel 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath_custom(fullfile(functionsdirectory, '/Adams_Functions/'));        
                                    % for set_fields_zero.m, restore_fields.m,
                                    % find_in_strings.m, my_closereq.m
    addpath_custom(fullfile(functionsdirectory, '/Downloaded_Functions/'));        
                                    % for rgb.m
end

%% Extract default values from outParams
parentDirectory = outParams.parentDirectory;    % parent directory for input/output folders
cellname = outParams.cellname;              % name of current neuron to be fitted
debugflag = outParams.debugflag;            % whether in debug mode
runmode = outParams.runMode;                % run mode
filenames = outParams.filenames;            %% TODO: copy description from singleneuronfitting.m
numswps = outParams.numswps;                %% TODO: copy description from singleneuronfitting.m
default_swpuse = outParams.default_swpuse;        %% TODO: copy description from singleneuronfitting.m
default_swpw = outParams.default_swpw;            %% TODO: copy description from singleneuronfitting.m
fitwin = outParams.fitwin;                %% TODO: copy description from singleneuronfitting.m
fitwinCpr = outParams.fitwinCpr;                %% TODO: copy description from singleneuronfitting.m
ltsProperties = outParams.ltsProperties;        %% TODO: copy description from singleneuronfitting.m
default_ltsuse = outParams.default_ltsuse;        %% TODO: copy description from singleneuronfitting.m
default_ltsw = outParams.default_ltsw;            %% TODO: copy description from singleneuronfitting.m
default_errratio = outParams.lts_to_swp_errratio;    %% TODO: copy description from singleneuronfitting.m
numparams = outParams.numparams;            %% TODO: copy description from singleneuronfitting.m
neuronparams_val = outParams.neuronparams;        %% TODO: copy description from singleneuronfitting.m
neuronparams_use = outParams.neuronparams_use;        %% TODO: copy description from singleneuronfitting.m
neuronparams_min = outParams.neuronparams_min;        %% TODO: copy description from singleneuronfitting.m
neuronparams_max = outParams.neuronparams_max;        %% TODO: copy description from singleneuronfitting.m
neuronparamnames = outParams.neuronparamnames;        %% TODO: copy description from singleneuronfitting.m
neuronparamclass = outParams.neuronparamclass;    %% TODO: copy description from singleneuronfitting.m
neuronparamislog = outParams.neuronparamislog;        %% TODO: copy description from singleneuronfitting.m
neuronparamispas = outParams.neuronparamispas;        %% TODO: copy description from singleneuronfitting.m
neuronparams_jit = outParams.neuronparams_jit;        %% TODO: copy description from singleneuronfitting.m
simplexparamnames = outParams.simplexparamnames;    %% TODO: copy description from singleneuronfitting.m
simplexparamsinit = outParams.simplexparams;        %% TODO: copy description from singleneuronfitting.m
autoparamnames = outParams.autoparamnames;        %% TODO: copy description from singleneuronfitting.m
autoparamsinit = outParams.autoparams;            %% TODO: copy description from singleneuronfitting.m
autologgerflag = outParams.autologgerflag;        %% TODO: copy description from singleneuronfitting.m

%% CREATE GUI
% Initialize and hide the GUI as it is being constructed
% For figures, default 'Units' is 'pixels'; for uipanels, default 'Units' is 'normalized'
%     Format: [left bottom width height]
hfig.gui = figure('Visible', 'off', 'Name', ['OPTIMIZERGUI (c) M3HA 2017 for ', cellname], ...
        'Tag', 'mainGUI', 'CloseRequestFcn', {@my_closereq, 'No', 'OPTIMIZERGUI'});
if doublecolumnflag == 0
    set(hfig.gui, 'OuterPosition', foutpos_singcol);
    p1pos = p1pos_singcol;                  % location and size of Panel 1: 'Error Function Maker'
    p2pos = p2pos_singcol;                  % location and size of Panel 2: 'NEURON Controller'
else
    set(hfig.gui, 'OuterPosition', foutpos_doubcol);
    p1pos = p1pos_doubcol;
    p2pos = p2pos_doubcol;
end
fpos = get(hfig.gui, 'Position');           % location and size of drawable area of 
                                            % figure 'OPTIMIZERGUI (c) M3HA 2017' in pixels

% Define row heights and bottom positions of rows, some based on number of sweeps and numel(ltsProperties)
numlts = numel(ltsProperties);
p1trh = rowh_text/(fpos(4)*p1pos(4));                   % height of each text row of Panel 1
p1lrh = rowh_lts/(fpos(4)*p1pos(4));                    % height of each LTS property row of Panel 1
p1tr1b = 1 - p1trh;                                     % bottom of 1st text row of Panel 1
p1srsp = (1 - p1trh*6 - p1lrh * numlts)/(numswps+2);    % spacing of each sweep row of Panel 1
p1srh = min(rowh_swp/(fpos(4)*p1pos(4)), p1srsp);       % height of each sweep row of Panel 1
p1tr2b = p1tr1b - p1srsp * numswps - p1srsp - p1trh;    % bottom of 2nd text row of Panel 1
p1tr3b = p1tr2b - p1srsp - p1trh;                       % bottom of 3rd text row of Panel 1
p1tr4b = p1tr3b - p1lrh * numlts - p1trh;               % bottom of 4th text row of Panel 1
p1tr5b = p1tr4b - p1trh;                                % bottom of 5th text row of Panel 1
p1tr6b = p1tr5b - p1trh;                                % bottom of 6th text row of Panel 1
p7prh = rowh_par/(fpos(4)*p2pos(4)*p4pos(4)*p7pos(4));  % height of each parameter row relative to Panel 7

% Construct Panel 1: 'Error Function Maker'
p1 = uipanel(hfig.gui, 'Title', 'Error Function Maker', ...
    'Position', p1pos, 'Tag', 'p1', 'BackgroundColor', bgcolor0);

% 1st text row: Header for sweep errors, all static text fields
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr1b 0.4 p1trh], 'String', 'Sweep name', ...
    'Style', 'text', 'Tag', 'errfuntext1', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.4 p1tr1b 0.1 p1trh], 'String', {'Use?'}, ...
    'Style', 'text', 'Tag', 'errfuntext2', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.475 p1tr1b 0.1 p1trh], 'String', {'w'}, ...
    'Style', 'text', 'Tag', 'errfuntext3', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.575 p1tr1b 0.125 p1trh], 'String', {'Ledge'}, ...
    'Style', 'text', 'Tag', 'errfuntext4', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.7 p1tr1b 0.125 p1trh], 'String', {'Redge'}, ...
    'Style', 'text', 'Tag', 'errfuntext5', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.825 p1tr1b 0.175 p1trh], 'String', {'RMSE'}, ...
    'Style', 'text', 'Tag', 'errfuntext6', 'BackgroundColor', bgcolor0);

% Sweep rows
for k = 1:numswps
    % Static text field for name of sweep
    uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr1b-k*p1srsp 0.425 p1srh], ...
            'String', filenames{k}, ... 
            'Style', 'text', 'Tag', ['swptext', num2str(k)], 'BackgroundColor', bgcolor0);
    % Check box for whether to use sweep; default == checked
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.425 p1tr1b-k*p1srsp 0.05 p1srh], ...
            'String', {''}, 'Style', 'checkbox', 'Tag', ['swpuse', num2str(k)], ...
            'Value', default_swpuse, 'Callback', {@update_swpuse, k});
    % Editable text field for weight of sweep; default == 1
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.475 p1tr1b-k*p1srsp 0.1 p1srh], ...
            'String', num2str(default_swpw), ...
            'Style', 'edit', 'Tag', ['swpw', num2str(k)], ...
            'Callback', {@update_swpw, k});
    % Editable text field for left edge of sweep; default: plotwin(1)
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.575 p1tr1b-k*p1srsp 0.125 p1srh], ...
            'String', {num2str(fitwin(k, 1))}, ...
            'Style', 'edit', 'Tag', ['swpleftedge', num2str(k)], ...
            'Callback', {@update_swpleftedge, k});
    % Editable text field for right edge of sweep; default: plotwin(2)
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.7 p1tr1b-k*p1srsp 0.125 p1srh], ...
            'String', {num2str(fitwin(k, 2))}, ...
            'Style', 'edit', 'Tag', ['swprightedge', num2str(k)], ...
            'Callback', {@update_swprightedge, k});
    % Static text field for sweep error; to be computed
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.825 p1tr1b-k*p1srsp 0.175 p1srh], ...
            'String', {'-'}, 'Style', 'text', 'Tag', ['swperrtext', num2str(k)], 'BackgroundColor', bgcolor0);
        
end

% 2nd text row: Weighted Total Sweep Error, all static text fields
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr2b 0.8 p1trh], ...
    'String', {'Weighted Total Sweep Error = '}, ...
    'Style', 'text', 'Tag', 'totswperrlabel', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.8 p1tr2b 0.2 p1trh], ...
    'String', {'-'}, ...
    'Style', 'text', 'Tag', 'totswperrtext', 'BackgroundColor', bgcolor0);

% 3rd text row: Header for LTS property errors, all static text fields
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr3b 0.4 p1trh], 'String', 'LTS Property', ...
    'Style', 'text', 'Tag', 'ltstext', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.4 p1tr3b 0.1 p1trh], 'String', {'Use?'}, ...
    'Style', 'text', 'Tag', 'ltstext', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.5 p1tr3b 0.15 p1trh], 'String', {'w'}, ...
    'Style', 'text', 'Tag', 'ltstext', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.65 p1tr3b 0.35 p1trh], 'String', {'Error'}, ...
    'Style', 'text', 'Tag', 'ltstext', 'BackgroundColor', bgcolor0);

% LTS property error rows
for k = 1:numel(ltsProperties)
    % Static text field for name of LTS property
    uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr3b-k*p1lrh 0.4 p1trh], ...
        'String', ltsProperties{k}, ...
        'Style', 'text', 'Tag', ['ltsnametext', num2str(k)], 'BackgroundColor', bgcolor0);
    % Check box for whether to use LTS property; default == default_ltsuse
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.425 p1tr3b-k*p1lrh 0.05 p1trh], ...
        'String', {''}, 'Style', 'checkbox', 'Tag', ['ltsuse', num2str(k)], ...
        'Value', default_ltsuse(k), 'Callback', {@update_ltsuse, k});
    % Editable text field for weight of LTS property; default == default_ltsw
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.5 p1tr3b-k*p1lrh 0.15 p1trh], ...
        'String', {num2str(default_ltsw(k))}, 'Style', 'edit', 'Tag', ['ltsw', num2str(k)], ...
        'Callback', {@update_ltsw, k});
    % Static text field for LTS property error; to be computed
    uicontrol(p1, 'Units', 'normalized', 'Position', [0.65 p1tr3b-k*p1lrh 0.35 p1trh], ...
        'String', {'-'}, 'Style', 'text', 'Tag', ['ltserrtext', num2str(k)], ...
        'BackgroundColor', bgcolor0);
end

% 4th text row: Total LTS property error, all static text fields
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr4b 0.8 p1trh], ...
    'String', {'Total LTS Err = '}, ...
    'Style', 'text', 'Tag', 'totltserrlabel', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.8 p1tr4b 0.2 p1trh], ...
    'String', {'-'}, ...
    'Style', 'text', 'Tag', 'totltserrtext', 'BackgroundColor', bgcolor0);

% 5th text row: LTS/TotSwpErr weight ratio, one static and one editable text field
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr5b 0.8 p1trh], ...
    'String', 'LTS/TotSwpErr weight ratio = ', ...
    'Style', 'text', 'Tag', 'errratiolabel', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.8 p1tr5b 0.2 p1trh], ...
    'String', num2str(default_errratio), 'Style', 'edit', 'Tag', 'errratio', ...
    'Callback', @update_errratio);

% 6th text row: Weighted Total Error, all static text fields
uicontrol(p1, 'Units', 'normalized', 'Position', [0 p1tr6b 0.8 p1trh], ...
    'String', {'Weighted Total Error = '}, ...
    'Style', 'text', 'Tag', 'toterrlabel', 'BackgroundColor', bgcolor0);
uicontrol(p1, 'Units', 'normalized', 'Position', [0.8 p1tr6b 0.2 p1trh], ...
    'String', {'-'}, ...
    'Style', 'text', 'Tag', 'toterrtext', 'BackgroundColor', bgcolor0);

% Construct Panel 2: 'NEURON Controller'
p2 = uipanel(hfig.gui, 'Title', 'NEURON Controller', ...
    'Position', p2pos, 'Tag', 'p2', 'BackgroundColor', bgcolor0);

% Construct Button Group 3: 'Mode Chooser'
p3 = uibuttongroup(p2, 'Position', p3pos, 'Tag', 'p3', ...
    'SelectionChangeFcn', @modebutton_selectionchange);
uicontrol(p3, 'Units', 'normalized', 'Position', [0 0 0.2 1], ...
    'String', ' MANUAL', 'Style', 'radiobutton', 'Tag', 'modebutton_manual');
                                % Button for MANUAL mode
uicontrol(p3, 'Units', 'normalized', 'Position', [0.2 0 0.4 1], ...
    'String', ' AUTO', 'Style', 'radiobutton', 'Tag', 'modebutton_auto');
                                % Button for AUTO mode
uicontrol(p3, 'Units', 'normalized', 'Position', [0.4 0 0.6 1], ...
    'String', ' JITTER', 'Style', 'radiobutton', 'Tag', 'modebutton_jitter');
                                % Button for JITTER mode
uicontrol(p3, 'Units', 'normalized', 'Position', [0.6 0 0.4 1], ...
    'String', ' AUTO WITH JITTER', 'Style', 'radiobutton', 'Tag', 'modebutton_auto_w_jitter');
                                % Button for AUTO WITH JITTER mode

% Construct Panel 4: 'NEURON parameters'
p4 = uipanel(p2, 'Title', 'NEURON parameters', ...
    'Position', p4pos, 'Tag', 'p4');

% Find the excess height of Panel 8 in Panel 7 and position vector
halfnumparams = ceil(numparams/2);      % about half of the number of NEURON parameters
p8Height = (p7prh*halfnumparams)/0.96;  % required height of Panel 8 in Panel 7
excessHeight = p8Height - 1;            % excess height of Panel 8 in Panel 7
p8pos = [0 -excessHeight 1 p8Height];   % position of Panel 8 in Panel 7
p8prh = rowh_par/(fpos(4)*p2pos(4)*p4pos(4)*p7pos(4)*p8Height);
                                        % height of each parameter row relative to Panel 8

% Construct Panel 7 for restricting the NEURON parameters
p7 = uipanel(p4, 'Position', p7pos, 'Tag', 'p7');

% Construct Panel 8 for placing the NEURON parameters
p8 = uipanel(p7, 'Position', p8pos, 'Tag', 'p8');

% Slider (if necessary) for moving Panel 7 through Panel 4
if excessHeight > 0
    uicontrol(p4, 'Units', 'normalized', 'Position', pos_p8slider, ...
                'Style', 'slider', 'Tag', 'p8slider', ...
                'Min', 0, 'Max', excessHeight, 'Value', excessHeight, ...
                'Callback', {@update_p8pos, p8});
end

% Prepare values for NEURON parameter rows
slider_min = neuronparams_min;      % minimum position of NEURON parameters on the slider
slider_max = neuronparams_max;      % maximum position of NEURON parameters on the slider
slider_val = neuronparams_val;      % initial position of NEURON parameters on the slider
for k = 1:numparams
    % Change positions on slider to log scale for designated parameters
    if neuronparamislog(k)
        slider_min(k) = log10(neuronparams_min(k));
        slider_max(k) = log10(neuronparams_max(k));
        slider_val(k) = log10(neuronparams_val(k));
    end
    if slider_min(k) > slider_max(k) || ...
        slider_min(k) > slider_val(k) || ...
        slider_max(k) < slider_val(k)
        error(['Initial slider value or bounds for the %dth ', ...
            'NEURON parameter is incorrect!'], k);
    end
end

% NEURON parameter rows
for k = 1:numparams
    switch neuronparamclass(k)
    case 1
        bgcolor = rgb(bgcolor1);
    case 2
        bgcolor = rgb(bgcolor2);
    case 3
        bgcolor = rgb(bgcolor3);
    case 4
        bgcolor = rgb(bgcolor4);
    case 5
        bgcolor = rgb(bgcolor5);
    case 6
        bgcolor = rgb(bgcolor6);
    case 7
        bgcolor = rgb(bgcolor7);
    case 8
        bgcolor = rgb(bgcolor8);
    case 9
        bgcolor = rgb(bgcolor9);
    case 10
        bgcolor = rgb(bgcolor10);
    case 11
        bgcolor = rgb(bgcolor11);
    otherwise
        bgcolor = rgb(bgcolor12);
    end

    % Find left and bottom in normalized units of name, checkbox, slider, value, respectively
    if doublecolumnflag == 1 && k > halfnumparams
        xshift = 0.5;           % shift of left position in normalized units
        ordercol = k - halfnumparams;   % order of parameter in this column
    else
        xshift = 0;
        ordercol = k;
    end
    xpos = xpos_param + xshift; % left position in normalized units of name, checkbox, slider, value
    ypos = 0.98 - p8prh .* [ordercol, (ordercol - 0.15), (ordercol - 0.3), ordercol];
                        % bottom position in normalized units of name, checkbox, slider, value
    % Static text field for name of the NEURON parameter
    uicontrol(p8, 'Units', 'normalized', 'Position', [xpos(1) ypos(1) xpos(2)-xpos(1) p8prh], ...
                'String', neuronparamnames{k}, 'Style', 'text', ...
                'Tag', ['paramtext', num2str(k)], 'BackgroundColor', bgcolor);
    % Check box for whether to use the NEURON parameter
    uicontrol(p8, 'Units', 'normalized', 'Position', [xpos(2) ypos(2) xpos(3)-xpos(2) p8prh], ...
                'String', {''}, 'Style', 'checkbox', 'Tag', ['paramuse', num2str(k)], ...
                'Value', neuronparams_use(k), ...
                'Callback', {@update_neuronparams_use, k});
    % Slider for tuning the value of the NEURON parameter
    uicontrol(p8, 'Units', 'normalized', 'Position', [xpos(3) ypos(3) xpos(4)-xpos(3) 0.75*p8prh], ...
                'String', {''}, 'Style', 'slider', 'Tag', ['paramslider', num2str(k)], ...
                'Min', slider_min(k), 'Max', slider_max(k), ...
                'Value', slider_val(k), 'BackgroundColor', bgcolor, ...
                'Callback', {@update_neuronparams, k});
    % Static text field for value of the NEURON parameter
    uicontrol(p8, 'Units', 'normalized', 'Position', [xpos(4) ypos(4) (xshift+0.5)-xpos(4) p8prh], ...
                'String', num2str(neuronparams_val(k)), 'Style', 'text', ...
                'Tag', ['paramvaltext', num2str(k)]);

end

% Button for 'Load P-file'
uicontrol(p4, 'Units', 'normalized', 'Position', pos_pfileload, ...
    'String', 'Load P-file', 'Style', 'togglebutton', 'Tag', 'parameterfileloadbutton', ...
    'BackgroundColor', buttoncolor_OFF, ...
    'Callback', {@parameterfileloadbutton_toggle, parentDirectory, buttoncolor_ON, buttoncolor_OFF});

% Button for 'Save P-file'
uicontrol(p4, 'Units', 'normalized', 'Position', pos_pfilesave, ...
    'String', 'Save P-file', 'Style', 'togglebutton', 'Tag', 'parameterfilesavebutton', ...
    'BackgroundColor', buttoncolor_OFF, ...
    'Callback', {@parameterfilesavebutton_toggle, parentDirectory, buttoncolor_ON, buttoncolor_OFF});

% Button for 'Change Bounds'
uicontrol(p4, 'Units', 'normalized', 'Position', pos_chbounds, ...
    'String', 'Change Bounds', 'Style', 'togglebutton', 'Tag', 'boundsbutton', ...
    'BackgroundColor', buttoncolor_OFF, ...
    'Callback', {@boundsbutton_toggle, buttoncolor_ON, buttoncolor_OFF});

% Construct Panel 5: 'AUTO mode'
p5 = uipanel(p2, 'Title', 'AUTO mode', ...
    'Position', p5pos, 'Tag', 'p5');

% Check box for whether to log simplex errors and parameters
uicontrol('Parent', p5, 'Units', 'normalized', 'Position', pos_autologger, ...
        'String', {' Log simplex errors and params'}, ...
        'Style', 'checkbox', 'Tag', 'autologger', 'Value', autologgerflag, ...
        'Callback', @update_autologger);

% Table of simplex parameters
uitable('Parent', p5, 'Units', 'normalized', 'Position', pos_autotable1, ...
        'Data', [simplexparamnames, num2cell(simplexparamsinit)], ...
        'Tag', 'autotable1', 'ColumnWidth', {100, 100}, 'ColumnEditable', [false, true], ...
        'CellEditCallback', @update_simplexparams);

% Table of other auto mode parameters
uitable('Parent', p5, 'Units', 'normalized', 'Position', pos_autotable2, ...
        'Data', [autoparamnames, num2cell(autoparamsinit)], ...
        'Tag', 'autotable2', 'ColumnWidth', {100, 100}, 'ColumnEditable', [false, true], ...
        'CellEditCallback', @update_autoparams);

% Construct Run Button Panel
p6 = uipanel(p2, 'Position', p6pos, 'Tag', 'runbuttonpanel', ...
        'BackgroundColor', bgcolor0, 'BorderType', 'none');

% Button for 'RUN'
runbutton = uicontrol(p6, 'Units', 'normalized', 'Position', pos_runbutton, ...
    'String', 'RUN', 'Style', 'togglebutton', 'Tag', 'runbutton', ...
    'BackgroundColor', buttoncolor_OFF, ...
    'Callback', {@runbutton_toggle, hfig, realdata_cpr, realData, ...
            ltsProperties, buttoncolor_ON, buttoncolor_OFF, debugflag});

% Button for 'TWIDDLE'
uicontrol(p6, 'Units', 'normalized', 'Position', pos_twiddlebutton, ...
    'String', 'TWIDDLE', 'Style', 'togglebutton', 'Tag', 'twiddlebutton', ...
    'BackgroundColor', buttoncolor_OFF, ...
    'Callback', {@twiddlebutton_toggle, hfig, realdata_cpr, realData, ...
            buttoncolor_ON, buttoncolor_OFF});

%% INITIATE GUI

% Set default button to be selected in Button Group 3 then execute modebutton_selectionchange
switch runmode
case 1
    set(p3, 'SelectedObject', findall(p3, 'Tag', 'modebutton_manual'));
case 2
    set(p3, 'SelectedObject', findall(p3, 'Tag', 'modebutton_auto'));
case 3
    set(p3, 'SelectedObject', findall(p3, 'Tag', 'modebutton_jitter'));
case 4
    set(p3, 'SelectedObject', findall(p3, 'Tag', 'modebutton_auto_w_jitter'));
otherwise
    set(p3, 'SelectedObject', findall(p3, 'Tag', 'modebutton_manual'));
    outParams.runMode = 1;
    fprintf('Warning: run mode out of range, changed to 1!\n\n');
end
modebutton_selectionchange(p3);

% Make GUI visible
set(hfig.gui, 'Visible', 'on');
fprintf('OPTIMIZER GUI for %s is ready!\n\n', cellname);
drawnow;                % Update GUI

if outParams.autorunflag
    %% RUN ONCE
    fprintf('RUNNING!\n\n');
    set(runbutton, 'Value', runbutton.Max);
    drawnow;            % Update GUI
    runbutton_toggle (runbutton, [], hfig, realdata_cpr, realData, ...
            ltsProperties, buttoncolor_ON, buttoncolor_OFF, debugflag);
    drawnow;            % Update GUI
end

fprintf('Close GUI when done ... \n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Callbacks for OPTIMIZERGUI  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_flags
%% Update flags when outParams.ltsw is changed

global outParams

% Determine whether LTS error needs to be computed
if all(outParams.ltsw == 0)
    % Not computed when no LTS property is checked or all have zero weights
    outParams.ltserrorflag = 0;     

    % Use user's findltsflag in this case
    outParams.findltsflag = outParams.findltsflag_user;
else
    % Computed when at least one LTS property is checked and has non-zero weight
    outParams.ltserrorflag = 1;   

    % Force find LTS in this case
    if ~outParams.findltsflag_user
        outParams.findltsflag = 1;
        fprintf('Warning: findltsflag changed to true because ltsw is not all 0!\n\n');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_sliders
%% On changing bounds/initial values, update sliders

global outParams

% Initialize new slider bounds and values
slidermin = outParams.neuronparams_min;
slidermax = outParams.neuronparams_max;
sliderval = outParams.neuronparams;

% Update slider bounds and values for each parameter
for k = 1:numel(outParams.neuronparams)
    % Change to log scale if designated
    if outParams.neuronparamislog(k)
        slidermin(k) = log10(outParams.neuronparams_min(k));
        slidermax(k) = log10(outParams.neuronparams_max(k));
        sliderval(k) = log10(outParams.neuronparams(k));
    end
    
    % Adapt current parameter value to new bounds
    if sliderval(k) < slidermin(k)        %if trying to set lower bound such that current param val is out of bounds
        % Force set current parameter value to lower bound value
        sliderval(k) = slidermin(k);
        outParams.neuronparams(k) = outParams.neuronparams_min(k);  
    elseif sliderval(k) > slidermax(k)
        % Force set current parameter value to upper bound value
        sliderval(k) = slidermax(k);
        outParams.neuronparams(k) = outParams.neuronparams_max(k);
    end

    % Update slider bounds and values
    paramslider_k = findobj('Tag', ['paramslider', num2str(k)]);
    paramvaltext_k = findobj('Tag', ['paramvaltext', num2str(k)]);
    set(paramslider_k, 'Min', slidermin(k));
    set(paramslider_k, 'Max', slidermax(k));
    set(paramslider_k, 'Value', sliderval(k));
    if outParams.neuronparamislog(k)
        set(paramvaltext_k, 'String', num2str(outParams.neuronparams(k), '%5.2e'));
    else
        set(paramvaltext_k, 'String', num2str(outParams.neuronparams(k), 3));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_errortext (err, ltsProperties)
%% On finishing a simulation/optimization, update error texts on GUI

global outParams

for k = 1:outParams.numswps
    set(findobj('Tag', ['swperrtext', num2str(k)]), 'String', num2str(err.swperr(k)));
end
set(findobj('Tag', 'totswperrtext'), 'String', num2str(err.avgswperr, 3));
amp_k = find_in_strings('amp', ltsProperties);
set(findobj('Tag', ['ltserrtext', num2str(amp_k)]), ...
            'String', num2str(err.avgltsverr, 3));
time_k = find_in_strings('time', ltsProperties);
set(findobj('Tag', ['ltserrtext', num2str(time_k)]), ...
            'String', num2str(err.avgltsterr, 3));
slope_k = find_in_strings('slope', ltsProperties);
set(findobj('Tag', ['ltserrtext', num2str(slope_k)]), ...
            'String', num2str(err.avgltsdvdtverr, 3));
set(findobj('Tag', 'totltserrtext'), 'String', num2str(err.avgltserr, 3));
set(findobj('Tag', 'toterrtext'), 'String', num2str(err.toterr, 3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_swpuse (hObject, ~, k)
%% On check box click of sweep use, update sweep weight in outParams

global outParams

% Find the corresponding sweep weight text field
p1 = get(hObject, 'Parent');
swpw_k = findall(p1, 'Tag', ['swpw', num2str(k)]);

% Update sweep weight in outParams
outParams.swpw(k) = get(hObject, 'Value') * str2double(get(swpw_k, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_swpw (hObject, ~, k)
%% On text field edit of sweep weight, update sweep weight in outParams

global outParams 

% Find the corresponding sweep use check box
p1 = get(hObject, 'Parent');
swpuse_k = findall(p1, 'Tag', ['swpuse', num2str(k)]);

% Update sweep weight in outParams
outParams.swpw(k) = get(swpuse_k, 'Value') * str2double(get(hObject, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_swpleftedge (hObject, ~, k)
%% On text field edit, update sweep left edge

global outParams 

outParams.fitwin(k, 1) = str2double(get(hObject, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_swprightedge (hObject, ~, k)
%% On text field edit, update sweep left edge

global outParams 

outParams.fitwin(k, 2) = str2double(get(hObject, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_ltsuse (hObject, ~, k)
%% On check box click of LTS use, update LTS weight in outParams

global outParams

% Find the corresponding LTS weight text field
p1 = get(hObject, 'Parent');
ltsw_k = findall(p1, 'Tag', ['ltsw', num2str(k)]);

% Update LTS weight in outParams
outParams.ltsw(k) = get(hObject, 'Value') * str2double(get(ltsw_k, 'String'));

% Update flags if needed
update_flags;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_ltsw (hObject, ~, k)
%% On text field edit of LTS weight, update LTS weight in outParams

global outParams

% Find the corresponding LTS use check box
p1 = get(hObject, 'Parent');
ltsuse_k = findall(p1, 'Tag', ['ltsuse', num2str(k)]);

% Update LTS weight in outParams
outParams.ltsw(k) = get(ltsuse_k, 'Value') * str2double(get(hObject, 'String'));

% Update flags if needed
update_flags;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_errratio (hObject, ~)
%% On text field edit, update LTS to sweep error ratio

global outParams

outParams.lts_to_swp_errratio = str2double(get(hObject, 'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_neuronparams (hObject, ~, k)
%% On slider movement, get value, update slider text, and update NEURON parameter value

global outParams

% Get corresponding text field
p4 = get(hObject, 'Parent');
paramvaltext_k = findall(p4, 'Tag', ['paramvaltext', num2str(k)]);

% Get slider value
sliderval = get(hObject, 'Value');    % the slider value may be log scaled

% Update slider text and value
if outParams.neuronparamislog(k)
    trueval = 10^sliderval;        % the true value 
    set(paramvaltext_k, 'String', num2str(trueval, '%5.2e'));
else
    trueval = sliderval;
    set(paramvaltext_k, 'String', num2str(trueval, 3));
end

% Update NEURON parameter value
outParams.neuronparams(k) = trueval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_neuronparams_use (hObject, ~, k)
%% On check box click, update whether NEURON parameter needs to be used

global outParams

outParams.neuronparams_use(k) = get(hObject, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_autologger (hObject, ~)
%% On check box click, update whether to log simplex errors and parameters

global outParams

outParams.autologgerflag = get(hObject, 'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_simplexparams (hObject, ~)
%% On table edit, update simplex parameters to be used

global outParams

tabdata = get(hObject, 'Data');
outParams.simplexparams = cell2mat(tabdata(:, 2));        % update simplex parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_autoparams (hObject, ~)
%% On table edit, update other auto mode parameters to be used

global outParams

tabdata = get(hObject, 'Data');
outParams.autoparams = cell2mat(tabdata(:, 2));            % update other auto mode parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameterfileloadbutton_toggle (hObject, ~, parentDirectory, buttoncolor_ON, buttoncolor_OFF)
%% On toggle of parameterfileloadbutton, open p-file

global outParams

togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
if togstate == togmax           % if toggled ON
    % Make button red
    set(hObject, 'BackgroundColor', buttoncolor_ON);

    % Load file with selection box
    % File must be a 4 column matrix: parameter names, values, lower bounds, upper bounds
    [pfilename, pfilepath] = uigetfile(fullfile(parentDirectory, '/optimizer4gabab/pfiles/*.p')); 
    if isequal(pfilename, 0)
        disp('User selected Cancel');
    else
        pfiledata = importdata(fullfile(pfilepath, pfilename));
        for k = 1:size(pfiledata, 1)
            pn = find(strcmp(outParams.neuronparamnames, pfiledata{k, 1}) == 1);
            if isempty(pn)                  % if parameter is not fitted anymore
                fprintf('The parameter %s is not fitted anymore!\n\n', pfiledata{k, 1});
            else
                outParams.neuronparams(pn) = pfiledata{k, 2};
                outParams.neuronparams_min(pn) = pfiledata{k, 3};
                outParams.neuronparams_max(pn) = pfiledata{k, 4};
            end
        end
        update_sliders;         % update slider bounds
    end
end
set(hObject, 'Value', togmin);                      % untoggle
set(hObject, 'BackgroundColor', buttoncolor_OFF);   % make button green

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parameterfilesavebutton_toggle (hObject, ~, parentDirectory, buttoncolor_ON, buttoncolor_OFF)
%% On toggle of parameterfilesavebutton, open selection box to save p-file

global outParams

togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
if togstate == togmax           % if toggled ON
    % Make button red
    set(hObject, 'BackgroundColor', buttoncolor_ON);

    % Load file with selection box
    [pfilename, pfilepath] = uiputfile(fullfile(parentDirectory, '/optimizer4gabab/pfiles/*.p'));
    if isequal(pfilename, 0)
        disp('User selected Cancel');
    else
        % Save a 4 column matrix: parameter names, values, lower bounds, upper bounds
        tosave = [outParams.neuronparamnames', num2cell(outParams.neuronparams'), ...
            num2cell(outParams.neuronparams_min'), num2cell(outParams.neuronparams_max')];
        save(fullfile(pfilepath, pfilename), 'tosave');
    end
end
set(hObject, 'Value', togmin);                      % untoggle
set(hObject, 'BackgroundColor', buttoncolor_OFF);   % make button green

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function boundsbutton_toggle (hObject, ~, buttoncolor_ON, buttoncolor_OFF)
%% On toggle of boundsbutton, open little GUI to update bounds and initial values

global outParams

togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
switch togstate
case togmax        % if toggled ON
    % Make button red and change text
    set(hObject, 'BackgroundColor', buttoncolor_ON);
    set(hObject, 'String', 'Changing bounds ...');

    % Open little GUI
    fbounds = figure('Visible', 'on', 'Name', 'Choose bounds and initial values', ...
        'Position', [500, 600, 450, 400], 'Tag', 'fbounds');

    % Table of parameter values and bounds
    uitable(fbounds, 'Data', ...
        [outParams.neuronparamnames', num2cell(outParams.neuronparams_min'), ...
        num2cell(outParams.neuronparams_max'), num2cell(outParams.neuronparams')], ...
        'Units', 'normalized', 'Position', [0.02, 0.275, 0.96, 0.7], 'Tag', 'boundstable', ...
        'ColumnWidth', {100, 90, 90, 90}, 'ColumnEditable', logical([0, 1, 1, 1]), ...
        'ColumnName', {'Parameter', 'Lower bound', 'Upper bound', 'Initial value'});

    % Button for 'Save & close'
    uicontrol(fbounds, 'Units', 'normalized', 'Position', [0.2, 0.08, 0.2, 0.1], ...
        'String', 'Save & close', 'Style', 'pushbutton', 'Tag', 'boundsdonebutton', ...
        'BackgroundColor', buttoncolor_OFF, ...
        'Callback', {@boundsdonebutton_push, hObject, buttoncolor_OFF});
    % Button for 'Cancel'
    uicontrol(fbounds, 'Units', 'normalized', 'Position', [0.6, 0.08, 0.2, 0.1], ...
        'String', 'Cancel', 'Style', 'pushbutton', 'Tag', 'boundscancelbutton', ...
        'BackgroundColor', buttoncolor_OFF, ...
        'Callback', {@boundscancelbutton_push, hObject, buttoncolor_OFF});
case togmin        % if toggled OFF
    fbounds = findobj('Tag', 'fbounds');
    if ~isempty(fbounds)                            % if little GUI still exists
        set(hObject, 'Value', 1);                   % toggle back ON
        warndlg('Bounds still being changed!', ...
            'Button press warning', 'modal');       % open warning dialog box
        figure(fbounds);                            % make current figure
    else
        set(hObject, 'String', 'Change Bounds')             % restore text
        set(hObject, 'BackgroundColor', buttoncolor_OFF);   % make button green
    end
otherwise
    error('togstate unrecognised!\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function boundsdonebutton_push (hObject, ~, boundsbutton, buttoncolor_OFF)
% Save new bounds and initial values upon pressing 'Save & close' of 'Choose bounds and initial values' GUI

global outParams

% Retrieve bounds and initial values set by user
tabdata = get(findobj('Tag', 'boundstable'), 'Data');

% Check bounds and initial values
for p = 1:outParams.numparams
    if ~(tabdata{p, 2} <= tabdata{p, 4} && tabdata{p, 4} <= tabdata{p, 3})
        % Pop up warning dialog box asking to change value for this parameter
        warndlg('Initial values outside of bounds!', ...
            'Invalid values warning', 'modal');             % open warning dialog box
        return;
    end
end

% Update bounds and initial values and slider bounds and positions
outParams.neuronparams_min = reshape(cell2mat(tabdata(:, 2)), 1, []);   % update lower bounds
outParams.neuronparams_max = reshape(cell2mat(tabdata(:, 3)), 1, []);   % update upper bounds
outParams.neuronparams = reshape(cell2mat(tabdata(:, 4)), 1, []);       % update initial vals
update_sliders;                                             % update slider bounds

% Reset 'Change Bounds' button
set(boundsbutton, 'Value', 0);                              % toggle OFF
set(boundsbutton, 'String', 'Change Bounds');               % restore text
set(boundsbutton, 'BackgroundColor', buttoncolor_OFF);      % make button green
close(get(hObject, 'Parent'));                              % close GUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function boundscancelbutton_push (hObject, ~, boundsbutton, buttoncolor_OFF)
% Upon pressing 'Cancel' of 'Choose bounds and initial values' GUI, reset button

% Reset 'Change Bounds' button
set(boundsbutton, 'Value', 0);                              % toggle OFF
set(boundsbutton, 'String', 'Change Bounds');               % restore text
set(boundsbutton, 'BackgroundColor', buttoncolor_OFF);      % make button green
close(get(hObject, 'Parent'));                              % close GUI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modebutton_selectionchange(hObject, ~)
%% On change of Button Group 3, update outParams.runMode and enable/disable things in Panel 4/5 (NEURON parameters/AUTO)

global outParams

% Get the selected button and store in modeselected
modeselected = get(get(hObject, 'SelectedObject'), 'Tag');

% Find Panels 4 and 5 from Panel 2, which is parent of Button Group 3
p2 = get(hObject, 'Parent');
p4 = findall(p2, 'Tag', 'p4');
p5 = findall(p2, 'Tag', 'p5');

% Enable/disable things in Panel 4/5 (NEURON parameters/AUTO) according to mode selected
p4children = allchild(p4);                                          % all children objects of p4
p4children_sliders = findall(p4children, 'Style', 'slider');        % all sliders of p4        
p4children_checkboxes = findall(p4children, 'Style', 'checkbox');   % all check boxes of p4
p5children = allchild(p5);                                          % all children objects of p5
switch modeselected
case 'modebutton_manual'                            % for MANUAL mode
    set(p4children_sliders, 'Enable', 'on');            % can slide parameters to tweak them
    set(p4children_checkboxes, 'Enable', 'off');        % cannot select parameters; gray out check boxes
    set(p5children, 'Enable', 'off');                   % disable and gray out AUTO panel
    outParams.runMode = 1;                              % update runmode
case 'modebutton_jitter'                            % for JITTER mode
    set(p4children_sliders, 'Enable', 'inactive');      % cannot slide parameters during jittering
    set(p4children_checkboxes, 'Enable', 'on');         % can select parameters to jitter
    set(p5children, 'Enable', 'off');                   % disable and gray out AUTO panel
    outParams.runMode = 3;                              % update runmode
case {'modebutton_auto', 'modebutton_auto_w_jitter'}% for both AUTO & AUTO WITH JITTER modes
    set(p4children_sliders, 'Enable', 'inactive');      % cannot slide parameters during 
    set(p4children_checkboxes, 'Enable', 'on');         % can select parameters to optimize
    set(p5children, 'Enable', 'on');                    % enable AUTO panel
    switch modeselected
    case 'modebutton_auto'
        outParams.runMode = 2;                          % update runmode
    case 'modebutton_auto_w_jitter'
        outParams.runMode = 4;                          % update runmode
    end
otherwise
    error('modeselected unrecognised!\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_p8pos (hObject, ~, p8)
%% On movement of p8slider, update Panel 8 position

% Get previous position of Panel 8 in Panel 7
p8pos = get(p8, 'Position');

% Change the bottom of Panel 8 position
p8pos(2) = -get(hObject, 'value');

% Set new position of Panel 8 in Panel 7
set(p8, 'Position', p8pos);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runbutton_toggle (hObject, ~, hfig, realdata_cpr, realData, ltsProperties, buttoncolor_ON, buttoncolor_OFF, debugflag)
%% On toggle of runbutton: ALL SYSTEMS GO (OR KILL)
%% TODO: Isolate run part to bypass GUI

global outParams

% Find all other buttons, sliders, checkboxes & editable text fields
p1 = findobj('Tag', 'p1');
p3 = findobj('Tag', 'p3');
p4 = findobj('Tag', 'p4');
p5 = findobj('Tag', 'p5');
p6 = findobj('Tag', 'p6');
p7 = findobj('Tag', 'p7');
p8 = findobj('Tag', 'p8');

% Find all other UIControls
p1children = allchild(p1);
p3children = allchild(p3);
p4children = allchild(p4);
p4nonpanelchildren = setdiff(allchild(p4), [p7; p8]);
p5children = allchild(p5);
p6otherchildren = setdiff(allchild(p6), hObject);
allotherUIC = vertcat(p1children, p3children, p4nonpanelchildren, p5children, p6otherchildren);

% Store enable/disable UIControls in UserData
for k = 1:length(allotherUIC)
    set(allotherUIC(k), 'UserData', get(allotherUIC(k), 'Enable'));
end

% Get toggle state of run button and act accordingly
togstate = get(hObject, 'Value');
togmax = get(hObject, 'Max');
togmin = get(hObject, 'Min');
switch togstate
case togmax                 % if user toggles ON then send parameters to m3ha_optimizer_4compgabab.m
                            % (which will call NEURON and run simulations in auto or manual mode;
                            % when m3ha_optimizer_4compgabab.m is done, runbutton will be toggled off)

    % Change text on the button to RUNNING and change the color to red
    set(hObject, 'String', 'RUNNING');                      % set button text to RUNNING
    set(hObject, 'BackgroundColor', buttoncolor_ON);        % set button color to red

    % Disable all other UIControls
    set(p1children, 'Enable', 'off');    
    set(p3children, 'Enable', 'off');    
    set(p4nonpanelchildren, 'Enable', 'inactive');    
    set(p5children, 'Enable', 'off');    
    set(p6otherchildren, 'Enable', 'off');    

    %##########
    %##############

    % Call m3ha_optimizer_4compgabab.m to SIMULATE or OPTIMIZE
    [done, outParams, hfig] = m3ha_optimizer_4compgabab(realdata_cpr, realData, outParams, hfig); 

    %##############
    %##########

    % If on AUTO or AUTO WITH JITTER modes, update slider position based on optimized parameters
    if outParams.runMode == 2 || outParams.runMode == 4
        update_sliders;                                     % update slider bounds
    end

    % Update error text on GUI
    update_errortext(outParams.err{outParams.runnumtotal}, ltsProperties);

    % Reset flags and button if done
    if done
        % Re-enable previously disabled buttons, sliders and checkboxes
        for k = 1:length(allotherUIC)
            set(allotherUIC(k), 'Enable', get(allotherUIC(k), 'UserData'));    
        end

        % Reset state, text and color of button
        set(hObject, 'Value', togmin);                      % toggle OFF
        set(hObject, 'String', 'RUN');                      % reset button text to RUN
        set(hObject, 'BackgroundColor', buttoncolor_OFF);   % reset button color to green
    else
        fprintf('Optimizer terminated with error!\n\n');
    end
case togmin                 % if user toggles OFF then toggle it back on, 
                            %   unless in debug mode, when the button is reset
    %%% TODO: Create a warning box for the user displaying the following message:
    %     'Simulation/optimization still running, are you sure you want to cancel?'
    %    with two buttons: 'Return', 'Cancel'
    %     if user hits 'Return', the warning box closes and the RUN button is toggled back on
    %     if user hits 'Cancel', the warning box closes, 
    %            m3ha_optimizer_4compgabab.m stops,                  %%%%% NOT SURE HOW TO DO THIS YET
    %            the last optimized parameters are returned,    %%%%% NOT SURE HOW TO DO THIS YET
    %            and the RUN button is reset
    %                                

    if ~debugflag    
        set(hObject, 'Value', togmax);                      % toggle back ON        
    else                % only allow resetting in DEBUG mode, lest one accidentally presses the button
        % Re-enable previously disabled UIControls
        for k = 1:length(allotherUIC)
            set(allotherUIC(k), 'Enable', get(allotherUIC(k), 'UserData'));    
        end

        % Reset text and color of button
        set(hObject, 'String', 'RUN');                      % reset button text to RUN
        set(hObject, 'BackgroundColor', buttoncolor_OFF);   % reset button color to green
    end
otherwise
    error('togstate unrecognised!\n\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function twiddlebutton_toggle (hObject, ~, hfig, realdata_cpr, realData, buttoncolor_ON, buttoncolor_OFF)
% Knob twiddling
% TODO: FIGURE OUT WHAT THIS DOES

global outParams

%[done err outParams] = optimizer(outfoldername,d, outParams,handles); % the business end of things
%outParams.err{outParams.runnum_auto + outParams.runnum_manual} = err;

%neuronparamnames = {'pcabar1', 'pcabar2', 'actshift_itGHK', 'shift_itGHK', 'ghbar1', 'ghbar2', 'shift_ih', 'gpas', 'epas'};

% twiddleparams{1} = [1e-5 1e-6 1e-7];
% twiddleparams{2} = [1e-5 1e-6 1e-7];
% twiddleparams{3} = [-10 -5 0 5 10];
% twiddleparams{4} = [-10 -5 0 5 10];
% twiddleparams{5} = [1e-5 1e-6 1e-7];
% twiddleparams{6} = [1e-5 1e-6 1e-7];
% twiddleparams{7} = [-10 -5 0 5 10];
% twiddleparams{8} = [1e-5 1e-6 1e-7];
% twiddleparams{9} = [-75 -70 -65 -60 -55];

% Don't plot figures if knobtwiddler is used
[outParams] = set_fields_zero(outParams, ...
            'plotzoomedflag', 'plotoverlappedflag', 'plotconductanceflag', 'plotcurrentflag', ...
            'plotipeakflag', 'plotLTSflag', 'plotstatisticsflag');

twiddleparams{1} = 5e-7;
twiddleparams{2} = 5e-7;
twiddleparams{3} = 4;
twiddleparams{4} = -4;
twiddleparams{5} = [5e-6 5e-5 5e-4];
twiddleparams{6} = [5e-6 5e-5 5e-4];
twiddleparams{7} = [-20 0 20];
twiddleparams{8} = [5e-6 5e-5 5e-4];
twiddleparams{9} = [-75 -60 -45];

initparams = zeros(1,9);
for k = 1:9
    initparams(k) = twiddleparams{k}(1);
end
nowparams = initparams;

%allerr = [];
for tw5 = 1:numel(twiddleparams{5})
    nowparams(5) = twiddleparams{5}(tw5);
    for tw6 = 1:numel(twiddleparams{6})
        nowparams(6) = twiddleparams{6}(tw6);
        for tw7 = 1:numel(twiddleparams{7})
            nowparams(7) = twiddleparams{7}(tw7);
            for tw8 = 1:numel(twiddleparams{8})
                nowparams(8) = twiddleparams{8}(tw8);
                for tw9 = 1:numel(twiddleparams{9})
                    nowparams(9) = twiddleparams{9}(tw9);
                    outParams.neuronparams = nowparams;
                                        
%                     [done err outParams] = optimizer_comp(outfoldername,d, outParams,handles); 
                     
                    %allerr = [allerr; err];
                end
            end
        end
    end
end

% Restore plot flags
[outParams] = restore_fields(outParams, ...
        'plotzoomedflag', 'plotoverlappedflag', 'plotconductanceflag', 'plotcurrentflag', ...
        'plotipeakflag', 'plotLTSflag', 'plotstatisticsflag');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
function local_CreateFcn (hObject, eventdata, createfcn, appdata)
% Utility function for GUI    %%% TODO: Understand what this can do

if ~isempty(appdata)
    names = fieldnames(appdata);
    for i = 1:length(names)
        name = char(names(i));
        setappdata(hObject, name, getfield(appdata, name));
    end
end

if ~isempty(createfcn)
    if isa(createfcn, 'function_handle')
        createfcn(hObject, eventdata);
    else
        eval(createfcn);
    end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

bgcolor1 = [0.9, 0.9, 0.9];        % background color of sliders

    if strcmp(outParams.modeselected, 'modebutton_auto') == 1 ...
        || strcmp(outParams.modeselected, 'modebutton_auto_w_jitter') == 1

    if neuronparamisbuild(k)
        bgcolor = rgb(bgcolor1);
    elseif neuronparamispas(k)
        bgcolor = rgb(bgcolor2);
    else
        bgcolor = rgb(bgcolor3);
    end
neuronparamisbuild = outParams.neuronparamisbuild;    %% TODO: copy description from singleneuronfitting.m

%}
