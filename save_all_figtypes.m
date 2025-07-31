function save_all_figtypes (fig, filename, varargin)
%% Save figures using all figure types provided
% Usage: save_all_figtypes (fig, filename, varargin)
% Arguments:
%       fig         - figure to save
%                   must be a a figure object or a Simulink block diagram
%       filename    - file name
%                   must be a character vector
%       figtypes    - (opt) figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in 
%                       saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%
% Requires:
%       cd/isfigtype.m
%
% Used by:
%       cd/compare_events_pre_post_stim.m
%       cd/compute_and_plot_all_responses.m
%       cd/compute_and_plot_average_response.m
%       cd/create_waveform_train.m
%       cd/create_pulse_train_series.m
%       cd/crosscorr_profile.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_plot_example_jitter.m.m
%       cd/m3ha_plot_bar3.m
%       cd/m3ha_plot_figure02.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure05.m
%       cd/m3ha_plot_figure08.m
%       cd/m3ha_plot_grouped_scatter.m
%       cd/m3ha_plot_violin.m
%       cd/m3ha_rank_neurons.m
%       cd/m3ha_simulate_population.m
%       cd/parse_current_family.m
%       cd/parse_pleth_trace.m
%       cd/plot_bar.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_chevron.m
%       cd/plot_fitted_traces.m
%       cd/plot_measures.m
%       cd/plot_relative_events.m
%       cd/plot_small_chevrons.m
%       cd/plot_struct.m
%       cd/plot_table_parallel.m
%       cd/plot_traces.m
%       cd/plot_traces_spike2_mat.m
%       cd/plot_tuning_curve.m
%       cd/plot_tuning_map.m
%       cd/save_all_zooms.m
%       cd/virt_golomb_generate_output.m
%       ~/minEASE/minEASE_compute_plot_average_psc.m
%       ~/minEASE/minEASE_detect_gapfree_events.m
%       ~/RTCl/raster_plot.m
%       ~/RTCl/single_neuron.m
%       ~/RTCl/tuning_curves.m
%       ~/RTCl/tuning_maps.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2017-05-09 Created by Adam Lu
% 2017-11-08 Replaced figbase with [figbase, '.', figtypes]
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-09-27 Made sure figure is visible when saved as a .fig file
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'fig');                         % figure to save
    % TODO: validation function %);
addRequired(iP, 'filename', ...                 % file name
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'figtypes', 'png', ...          % figure type(s) for saving
    @(x) min(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, fig, filename, varargin{:});
[~, figtypes] = isfigtype(iP.Results.figtypes, 'ValidateMode', true);

%% Preparation
% Make the figure visible if at least one of the figtypes is 'fig'
if any(strcmpi('fig', figtypes))
    set(fig, 'Visible', 'on');
end

%% Save figure(s)
if ~isempty(figtypes)    % if at least one figtype is provided
    % Break down file name
    [directory, figbase, ~] = fileparts(filename);

    % Save as figtypes
    if iscell(figtypes)        % if many figtypes provided
        % Save figure as each figtype
        nfigtypes = numel(figtypes);    % number of figure types for saving
        for f = 1:nfigtypes
            % Set the figure name
            figName = fullfile(directory, ...
                                [figbase, '.', figtypes{f}]);
            
            % Get the current figure type
            figType = figtypes{f};
            
            % Save figure as the figtype
            saveas(fig, figName, figType);
        end
    elseif ischar(figtypes)        % if only one figtype provided
        % Set the figure name
        figName = fullfile(directory, [figbase, '.', figtypes]);
        
        % Get the figure type
        figType = figtypes;
        
        % Save figure as the figtype
        saveas(fig, figName, figType);
    end
else            % if no figtype provided
    % Set the figure name as the file name
    figName = filename;
    
    % Simply use saveas()
    saveas(fig, figName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
