function eRev = compute_eRev (cOutArray, cInArray, pArray, zArray, celsius)
% Compute the reversal potential of a channel that passes monovalent ions using the GHK voltage equation
% Usage: eRev = compute_eRev (cOutArray, cInArray, pArray, zArray, celsius)
% Examples: 
%       For a channel selective for a single ion:
%           eNa = compute_eRev (naOut, naIn, 1, 1, celsius);
%       For a channel selective for multiple monovalent ions:
%           eRev = compute_eRev ([naOut, kOut, clOut], [naIn, kIn, clIn], ...
%                           [pNa, pK, pCl], [1, 1, -1], celsius);
%           eh = compute_eRev ([127.25, 2.5], [4.5, 113], [1, 3], [1, 1], 33);
%
% Arguments:
%       TODO
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/compute_fixed_params.m
%
% File History:
% 2017-08-05 Created

% Universal Constants
F = 96485.309;      % Faraday's constant [C/mol], used by NEURON
R = 8.31441;        % Universal gas constant [J/K mol], used by NEURON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 5
    error(['Not enough input arguments, ', ...
            'type ''help compute_Erev'' for usage']);
end

% Make sure all arrays are the same length
if ~isequal(length(cOutArray), length(cInArray), length(pArray), length(zArray))
    error('Input arrays not the same length!');
end

%% Do the job

% Convert temperature from Celsius to Kelvin
T = 273.15 + celsius; % temperature [Kelvin]

% Find out how many ions are passed in
nIons = length(cOutArray);
 
% Determine what concentrations are in the numerator vs in the denominator
cNum = cOutArray;       % initialize all numerator concentrations to Cout
cDen = cInArray;        % initialize all denominator concentrations to Cin
for iIon = 1:nIons
    % If the valence is -1, switch the concentrations
    if zArray(iIon) == -1
        cNum(iIon) = cInArray(iIon);
        cDen(iIon) = cOutArray(iIon);
    end
end

% Compute the numerator
numerator = sum(pArray .* cNum);

% Compute the denominator
denominator = sum(pArray .* cDen);

% Compute reversal potential
eRev = 1000 * (R*T/F) * log(numerator/denominator); % Reversal potential [mV]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
