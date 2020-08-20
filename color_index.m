function index = color_index(value, edges)
%% Find the colormap index for a given value with boundaries set by edges
% Usage: index = color_index(value, edges)
%
% Used by:
%		/media/adamX/Paula_IEIs/paula_iei4.m
%       cd/ZG_fit_IEI_distributions.m
%
% 2017-12-18 Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index = NaN * ones(size(value));
for i = 1:size(value, 1)
    for j = 1:size(value, 2)
        for idx = 1:length(edges)-1
            if value(i, j) >= edges(idx) && value(i, j) < edges(idx + 1)
                index(i, j) = idx;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
