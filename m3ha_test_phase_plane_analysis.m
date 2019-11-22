% m3ha_test_phase_plane_analysis.m
%% Runs phase plane analysis through example m3ha traces
%
% Requires:
%       cd/m3ha_phase_plane_analysis.m

% File History:
% 2016-07-30 Created

true_LTS = {'A092810_0001_12.mat', 'A092810_0002_2.mat', 'A092110_0012_5.mat',  ...
		'A092110_0012_10.mat', 'A092110_0012_25.mat', 'D091710_0007_15.mat', ...
		'A092810_0005_24.mat', 'A092810_0006_8.mat', 'A092810_0007_24.mat', ...
		'A092810_0007_25.mat', 'A092110_0005_17.mat', 'A092910_0002_22.mat'};
false_LTS = {'A092810_0002_7.mat', 'A092110_0005_2.mat', 'A092110_0013_19.mat', ...
		'A092110_0013_21.mat', 'A092810_0004_8.mat', 'A092810_0005_12.mat', ...
		'A092810_0006_3.mat', 'A092910_0001_12.mat', 'A092910_0001_16.mat', ...
		'A092910_0002_1.mat', 'A092110_0012_18.mat', 'A092810_0004_23.mat', ...
		'A092810_0005_8.mat'};
ltswin = [1100 8900];

%% True LTS
parfor f = 1:numel(true_LTS)
	filename = true_LTS{f};
	phase_plane_analysis(filename, ltswin, true);
end

%% False LTS
parfor f = 1:numel(false_LTS)
	filename = false_LTS{f};
	phase_plane_analysis(filename, ltswin, false);
end


