function ax = remove_non_axes (ax)
%% Removes axes that are titles or labels
% Usage: ax = remove_non_axes (ax)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       ax          - TODO: Description of ax
%                   specified as a TODO
%
% Arguments:
%       ax          - TODO: Description of ax
%                   must be a TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/align_subplots.m
%       cd/update_figure_for_corel.m

% File History:
% 2020-04-22 Moved from update_figure_for_corel.m
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% Extract tags
tags = extract_fields(ax, 'Tag');

% If there are no tags, return
if isempty(tags) 
    return
end

% Hard-coded parameters
tagsToRemove = {'suplabel', 'suptitle'};

% Decide whether to remove axes
toRemove = contains(tags, tagsToRemove);

% Remove the axes containing tags
ax = ax(~toRemove);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%