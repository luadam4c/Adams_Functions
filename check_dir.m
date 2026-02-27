function check_dir (directories, varargin)
%% Checks if needed directory(ies) exist and creates them if not
% Usage: check_dir (directories, varargin)
% Arguments:
%       directories - directory(ies) to check
%                   must be a cell array of character arrays
%                       or a scalar text
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'MessageMode' - how message boxes are shown
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'wait'  - stops program and waits for the user
%                                   to close the message box
%                       'show'  - does not stop program but still show the
%                                   message box
%                       'none'  - neither stop program nor show a message box
%                   default == 'wait'
%
% Requires:
%       cd/construct_fullpath.m
%       cd/print_or_show_message.m
%
% Used by:
%       cd/abf2mat.m
%       cd/archive_dependent_scripts.m
%       cd/atf2sheet.m
%       cd/backup_folders.m
%       cd/check_subdir.m
%       cd/compute_and_plot_evoked_LFP.m
%       cd/create_input_file.m
%       cd/create_new_mscript.m
%       cd/crosscorr_profile.m
%       cd/log_arraytext.m
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_initial_slopes.m
%       cd/m3ha_estimate_passive_params.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_optimizer_4compgabab.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_rank_neurons.m
%       cd/parse_all_abfs.m
%       cd/parse_current_family.m
%       cd/parse_multiunit.m
%       cd/plot_all_abfs.m
%       cd/plot_measures.m
%       cd/plot_traces_abf.m
%       cd/plot_traces_EEG.m
%       cd/virt_analyze_sniff_whisk.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting42.m
%       /home/Matlab/function_template.m
%       vIRt-Moore\virt_moore.m
%       vIRt-Moore\virt_run_monte_carlo_simulations.m
%       scAAV\analyze_reachr_motion.m
%       scAAV\analyze_qupath_ddio.m

% File History:
% 2018-06-21 Modified from check_subdir.m
% 2018-09-18 Added input parser and verbose
% 2018-10-03 Now uses print_or_show_message
% 2018-10-03 Now uses isfolder()

%% Hard-coded parameters
validMessageModes = {'wait', 'show', 'none'};
mtitle = 'New Directory Made';              % message title

%% Default values for optional arguments
verboseDefault = true;              % print to standard output by default
messageModeDefault = 'none';        % do not display message box by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'directories', ...
    @(x) iscellstr(directories) || ischar(directories) || ...
        isstring(directories));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MessageMode', messageModeDefault, ...
    @(x) any(validatestring(x, validMessageModes)));

% Read from the Input Parser
parse(iP, directories, varargin{:});
verbose = iP.Results.Verbose;
messageMode = validatestring(iP.Results.MessageMode, validMessageModes);

%% Check directory(ies)
if iscell(directories)
    for k = 1:numel(directories)        
        % Construct the full path to the directory
        directory = construct_fullpath(directories{k});

        % Check the directory and create it if it doesn't already exist
        check_dir_helper(directory, mtitle, messageMode, verbose);
    end
else
    % Construct the full path to the directory
    directory = construct_fullpath(directories);
    
    % Check the directory and create it if it doesn't already exist
    check_dir_helper(directory, mtitle, messageMode, verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_dir_helper(directory, mtitle, messageMode, verbose)

% Check if the directory exists
if ~isfolder(directory)
    % Create the directory
    mkdir(directory);

    % Show message and print to standard output
    msg = sprintf('New directory is made: %s\n\n', directory);
    print_or_show_message(msg, 'MTitle', mtitle, ...
                                'MessageMode', messageMode, ...
                                'Verbose', verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

if verbose
    fprintf('New directory is made: %s\n\n', directory);
end

messageMode= 'show';

if exist(directory, 'dir') ~= 7

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
