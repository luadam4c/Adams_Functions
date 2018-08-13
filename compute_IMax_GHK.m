function IMax = compute_IMax_GHK (v, pbar, z, celsius, cOut, cIn)
%% Computes the maximum current [mA/cm^2] using the GHK current equation
% Usage: IMax = compute_IMax_GHK (v, pbar, z, celsius, cOut, cIn)
%
% Arguments:
%       TODO
%
% Used by:    
%       /media/adamX/m3ha/optimizer4gabab/plot_IT.m
%
% File History:
% 2017-08-06 Created

% Conversion Constants
V_PER_mV = 1e-3;
mmol_PER_umol = 1e-3;

% Universal Constants
F = 96485.309;      % Faraday's constant [C/mol], used by NEURON
R = 8.31441;        % Universal gas constant [J/K mol], used by NEURON

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 5
    error(['Not enough input arguments, ', ...
            'type ''help compute_IMax'' for usage']);
end

%% Do the job

% Convert temperature from Celsius to Kelvin
T = 273.15 + celsius; % temperature [Kelvin]

% Compute the exponent zFV/RT, which is dimensionless
exponent = V_PER_mV * z.*F.*v./(R.*T);

% Compute the maximum current [mA/cm^2]
% NOTE: [mM] = [umol/cm3]
%       h(x) = x/(exp(x)-1), which is also dimensionless
IMax = mmol_PER_umol * pbar .* z .* F ...
        .* (cIn .* h(-exponent) - cOut .* h(exponent));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = h(x)
% the function x/(exp(x)-1) with approximation when x is small

val = zeros(size(x));
for i = 1:length(x)
    if abs(x(i)) < 1e-4
        % Apply 1st order Taylor approximation to x/(exp(x)-1)
        val(i) = 1 - x(i)/2;
    else
        % Just use x/(exp(x)-1)
        val(i) = x(i)/(exp(x(i))-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

Without using Taylor approximation:
IMax = mmol_PER_umol * pbar .* (z.*F.*exponent) ...
        .*(cIn-cOut.*exp(-exponent))./(1-exp(-exponent));

From IT.mod: 

FUNCTION ghk(v (mV), ci (mM), co (mM)) (.001 coul/cm3) {
    LOCAL z, eci, eco

    z = (1e-3) * 2 * FARADAY * v / (R * (celsius + 273.15 (degC) ) )
                            : this is ZFV/RT, which is dimensionless 
                            : after applying conversion factor 1e-3 V/mV
    eco = co * efun(z)      : this has units of [mM] = [umol/cm3]
    eci = ci * efun(-z)     : this has units of [mM] = [umol/cm3]
    ghk = (1e-3) * 2 * FARADAY * (eci - eco)
                            : this has units of [mC/cm3]
                            : after applying conversion factor 1e-3 mmol/umol
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1 - z/2      : 1st order Taylor approximation of z/(exp(z) - 1)
    } else {
        efun = z / (exp(z) - 1)
    }
}


%}
