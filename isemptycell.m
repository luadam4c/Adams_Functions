function isEmpty = isemptycell (cellArray)
%% Returns whether each cell of a cell array is empty; if not a cell array, same as isempty()
% Usage: isEmpty = isemptycell (cellArray)
% Explanation:
%       TODO
% Example(s):
%       isemptycell({1:5, [], 2})
%       isemptycell({'', [], 'sdgs'})
% Outputs:
%       isEmpty     - whether each cell or a cell array is empty
%                   specified as a logical array
% Arguments:
%       cellArray   - array or cell array to test
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/all_data_files.m
%       cd/combine_strings.m
%       cd/compute_rms_error.m
%       cd/compute_stats.m
%       cd/create_indices.m
%       cd/distribute_balls_into_boxes.m
%       cd/extract_common_prefix.m
%       cd/extract_distinct_fileparts.m
%       cd/find_in_strings.m
%       cd/istype.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/plot_fitted_traces.m
%       cd/m3ha_select_cells.m
%       cd/plot_traces.m
%       cd/remove_empty.m
%       cd/set_axes_properties.m
%       cd/set_figure_properties.m
%       cd/transform_vectors.m

% File History:
% 2019-01-04 Created by Adam Lu
% 2019-01-10 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
if iscell(cellArray) && ~isempty(cellArray)
    isEmpty = cellfun(@isempty, cellArray);
else
    isEmpty = isempty(cellArray);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
