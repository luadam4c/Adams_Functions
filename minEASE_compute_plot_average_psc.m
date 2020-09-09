function minEASE_compute_plot_average_psc (eventInfo, eventClass, ...
                data, siMs, traceLengthMs, beforePeakMs, ...
                dealWithTooShort, outputDirectory, outputLabel, varargin)
%% Computes and plots averaged Type I, II & III PSC traces
% Usage: minEASE_compute_plot_average_psc (eventInfo, eventClass, ...
%               data, siMs, traceLengthMs, beforePeakMs, ...
%               dealWithTooShort, outputDirectory, outputLabel, varargin)
% Outputs:
%   TODO
%
% Arguments:
%   TODO
%       varargin    - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'FigTypes': figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/compute_average_psc_trace.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%       cd/print_or_show_message.m
%
% Used by:
%       cd/minEASE.m
%       cd/minEASE_detect_gapfree_events.m

% File History:
% 2017-07-25 AL - Moved from minEASE_detect_gapfree_events.m
% 2018-02-02 AL - Added messageMode as an optional parameter-value pair argument
% 2018-02-02 AL - Now uses print_or_show_message.m for output
% 2018-02-07 MD - Changed usage of print_or_show_message()
% 2018-02-27 AL - Changed messageModes to messageMode with possible values:
%                   'wait', 'show', 'none'
% 2018-03-02 MD - Defined verbose parameter for print_or_show_message
% 2018-08-03 AL - Renamed sweepLabel -> outputLabel
% 2018-09-17 AL - Now uses save_all_figtypes.m
%

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};

%% Default values for optional arguments
messageModeDefault = 'none';    % print to standard output by default
verboseDefault = false;             % default: Program does not print message
                                    %   even if message box is shown
figTypesDefault = 'png';        % default figure type(s) for saving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) any(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, varargin{:});
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);
verbose = iP.Results.Verbose;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Compute the "average Type I PSC trace"
typeOneInfo = eventInfo(eventClass == 1, :);
if ~isempty(typeOneInfo)
    [averageTypeOneTrace, allTypeOneTraces] = ...
        compute_average_psc_trace(typeOneInfo, data, 'SiMs', siMs, ...
                                    'TraceLengthMs', traceLengthMs, ...
                                    'BeforePeakMs', beforePeakMs, ...
                                    'DealWithTooShort', dealWithTooShort, ...
                                    'MessageMode', messageMode);
end

% Prepare for plotting
message = '';                                   % no warning message at start
nSamplesAvgTrace = length(averageTypeOneTrace); % number of sample points
                                                %   for average PSC trace
tVecAvgTrace = linspace(siMs, nSamplesAvgTrace * siMs, nSamplesAvgTrace);
                                                % time vector (ms)
                                                %   for average PSC trace

% Plot "average Type I PSC trace"
if isempty(typeOneInfo) || isempty(allTypeOneTraces)
    % Add to warning message
    message = {message, 'No Type I PSCs satisfying criteria are found!'};
else
    % Plot averaged Type I PSC trace
    fig3 = figure(3);
    clf(fig3);
    title(['Average Type I PSC Trace (''', dealWithTooShort, ''' mode)']);
    hold on;

    % Plot all PSC traces to be averaged
    plot(tVecAvgTrace, allTypeOneTraces);

    % Plot average PSC trace in cyan bolded
    plot(tVecAvgTrace, averageTypeOneTrace, '-c', 'LineWidth', 4);
 
    % Change x limits
    xlim([0, traceLengthMs]);

    % Save figure
    figName3 = fullfile(outputDirectory, ...
                ['Average_Type_I_PSC_Trace_', outputLabel, ...
                 '_mode_', dealWithTooShort]);
    save_all_figtypes(fig3, figName3, figTypes);

    % Close figure
    % close(fig3);
end 

% Compute the "average Type II PSC trace"
typeTwoInfo = eventInfo(eventClass == 2, :);
if ~isempty(typeTwoInfo)
    [averageTypeTwoTrace, allTypeTwoTraces] = ...
        compute_average_psc_trace(typeTwoInfo, data, 'SiMs', siMs, ...
                                    'TraceLengthMs', traceLengthMs, ...
                                    'BeforePeakMs', beforePeakMs, ...
                                    'DealWithTooShort', dealWithTooShort, ...
                                    'MessageMode', messageMode);
end

% Plot "average Type II PSC trace"
if isempty(typeTwoInfo) || isempty(allTypeTwoTraces)
    % Add to warning message
    message = {message, 'No Type II PSCs satisfying criteria are found!'};
else
    % Plot averaged Type II PSC trace
    fig4 = figure(4);
    clf(fig4);
    title(['Average Type II PSC Trace (''', dealWithTooShort, ''' mode)']);
    hold on;

    % Plot all PSC traces to be averaged
    plot(tVecAvgTrace, allTypeTwoTraces);

    % Plot average PSC trace in cyan bolded
    plot(tVecAvgTrace, averageTypeTwoTrace, '-c', 'LineWidth', 4);

    % Change x limits
    xlim([0, traceLengthMs]);
 
    % Save figure
    figName4 = fullfile(outputDirectory, ...
                ['Average_Type_II_PSC_Trace_', outputLabel, ...
                 '_mode_', dealWithTooShort]);
    save_all_figtypes(fig4, figName4, figTypes);

    % Close figure
    % close(fig4);
end

% Compute the "average Type III PSC trace"
typeThreeInfo = eventInfo(eventClass == 3, :);
if ~isempty(typeThreeInfo)
    [averageTypeThreeTrace, allTypeThreeTraces] = ...
        compute_average_psc_trace(typeThreeInfo, data, 'SiMs', siMs, ...
                                    'TraceLengthMs', traceLengthMs, ...
                                    'BeforePeakMs', beforePeakMs, ...
                                    'DealWithTooShort', dealWithTooShort, ...
                                    'MessageMode', messageMode);
end

% Plot "average Type III PSC trace"
if isempty(typeThreeInfo) || isempty(allTypeThreeTraces)
    % Add to warning message
    message = {message, 'No Type III PSCs satisfying criteria are found!'};
else
    % Plot averaged Type III PSC trace
    fig5 = figure(5);
    clf(fig5);
    title(['Average Type III PSC Trace (''', dealWithTooShort, ''' mode)']);
    hold on;

    % Plot all PSC traces to be averaged
    plot(tVecAvgTrace, allTypeThreeTraces);

    % Plot average PSC trace in cyan bolded
    plot(tVecAvgTrace, averageTypeThreeTrace, '-c', 'LineWidth', 4);
 
    % Change x limits
    xlim([0, traceLengthMs]);

    % Save figure
    figName5 = fullfile(outputDirectory, ...
                ['Average_Type_III_PSC_Trace_', outputLabel, ...
                 '_mode_', dealWithTooShort]);
    save_all_figtypes(fig5, figName5, figTypes);

    % Close figure
    % close(fig5);
end

% Display warning message if any
if ~isempty(message)
    mTitle = 'Average Type I PSC Trace Warning';
    icon = 'warn';
    print_or_show_message(message, 'MessageMode', messageMode, ...
                            'MTitle', mTitle, 'Icon', icon, 'Verbose', verbose);
end

%% Save psc traces
if strcmp(dealWithTooShort, 'none') || strcmp(dealWithTooShort, 'omit')
    textPath = fullfile(outputDirectory, ...
                    ['All_Type_I_PSC_Traces_', outputLabel, ...
                     '_mode_', dealWithTooShort, '.txt']);

    % Save as an AXON Plain Text File without adding a time column
    dlmwrite(textPath, allTypeOneTraces, 'delimiter', '\t');

    textPath = fullfile(outputDirectory, ...
                    ['All_Type_II_PSC_Traces_', outputLabel, ...
                     '_mode_', dealWithTooShort, '.txt']);

    % Save as an AXON Plain Text File without adding a time column
    dlmwrite(textPath, allTypeTwoTraces, 'delimiter', '\t');

    textPath = fullfile(outputDirectory, ...
                    ['All_Type_III_PSC_Traces_', outputLabel, ...
                     '_mode_', dealWithTooShort, '.txt']);

    % Save as an AXON Plain Text File without adding a time column
    dlmwrite(textPath, allTypeThreeTraces, 'delimiter', '\t');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
