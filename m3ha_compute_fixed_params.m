function [eNa, eK, ehHigh, ehLow] = m3ha_compute_fixed_params(varargin)
%% Compute fixed parameters that are used in the model
% Usage: [eNa, eK, ehHigh, ehLow] = m3ha_compute_fixed_params(varargin)
%
% Requires:    
%       /home/Matlab/Adams_Functions/compute_eRev.m
%
% File History:
% 2017-08-05 - Created
% 2018-01-24 - Added isdeployed

%% Conditions of the experiments
celsius = 33;           % temperature of the experiment [degC]

%% Ion concentrations used in the experiments
naOut = 127.25;         % sodium concentration outside the cell [mM]
naIn = 4.5;             % sodium concentration inside the cell [mM]
kOut = 2.5;             % potassium concentration outside the cell [mM]
kIn = 113;              % potassium concentration inside the cell [mM]

%% H channel
pNaIh = 1;              % permeability of Na set to 1
pKIhLow = 3;            % low end of permeability of K
pKIhHigh = 4;           % high end of permeability of K 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions across servers
if exist('/home/Matlab/', 'dir') == 7
    functionsDirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsDirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsDirectory does not exist!');
end
if ~isdeployed
    addpath_custom(fullfile(functionsDirectory, '/Adams_Functions/')); 
                                        % for compute_eRev.m
end
%% Compute reversal potential for Na
eNa = compute_eRev(naOut, naIn, 1, 1, celsius);

fprintf('eNa = %g mV\n\n', eNa);

%% Compute reversal potential for K
eK = compute_eRev(kOut, kIn, 1, 1, celsius);

fprintf('eK = %g mV\n\n', eK);

%% Compute reversal potential for H channel
ehHigh = compute_eRev([naOut, kOut], [naIn, kIn], ...
                    [pNaIh, pKIhLow], [1, 1], celsius);
ehLow = compute_eRev([naOut, kOut], [naIn, kIn], ...
                    [pNaIh, pKIhHigh], [1, 1], celsius);

fprintf('ehHigh = %g mV\n', ehHigh);
fprintf('ehLow = %g mV\n\n', ehLow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

