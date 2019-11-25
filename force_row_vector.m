function vectors = force_row_vector (vectors, varargin)
%% Transform column vector(s) or array(s) to row vector(s)
% Usage: vectors = force_row_vector (vectors, varargin)
% Explanation:
%       Starting with a cell array of mixed vectors, some row and some column,
%           this function makes sure each vector is a row vector.
%       If a single vector is provided, the function makes sure 
%           it's a row vector.
%
% Example(s):
%       vector = force_row_vector(vector);
%       vectors = force_row_vector(vectors);
%       force_row_vector({[3, 4], [5; 6], magic(3)})
%
% Outputs:
%       vectors     - vectors transformed
%
% Arguments:
%       vectors     - original vectors
%       varargin    - see force_column_vector.m
%
% Requires:
%       cd/compute_psth.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%
% Used by:    
%       cd/compute_centers_from_edges.m
%       cd/compute_relative_event_times.m
%       cd/plot_tuning_curve.m
%       cd/m3ha_xolotl_plot.m
%       cd/plot_histogram.m
%       cd/run_neuron.m
%       cd/unique_groups.m

% File History:
% 2018-10-25 Modified from force_column_vector.m
% 2018-10-27 Added 'IgnoreNonVectors' as an optional argument
% 2019-01-13 Now uses force_column_vector.m with 'RowInstead' option true
% 

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
addRequired(iP, 'vectors');

% Read from the Input Parser
parse(iP, vectors, varargin{:});

%% Do the job
vectors = force_column_vector(vectors, 'RowInstead', true, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
