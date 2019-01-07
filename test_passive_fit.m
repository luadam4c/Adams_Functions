%% Tests passive fitting 
%
% Requires:
%       cd/load_examples.m
%       cd/parse_pulse_response.m

% File History:
% 2018-10-11 Created by Adam Lu
% 

% Load examples (see load_examples.m)
load_examples;

% Choose from loaded examples
pulseVectors = myPulse1;
pulseResponse = myPulseResponse1b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse the pulse response to get the pulse width
pulseParams = parse_pulse(pulseVectors);
pulseWidth = pulseParams.pulseWidthSamples * siMs;

% Parse the pulse response
[responseParams, responseData] = ...
    parse_pulse_response(pulseResponse, siMs, 'PulseVectors', pulseVectors);

% Extract the rising and falling phase vectors as a numeric array
tvecsRising = horzcat(responseData.tvecsRising{:});
vvecsRising = horzcat(responseData.vvecsRising{:});
tvecsFalling = horzcat(responseData.tvecsFalling{:});
vvecsFalling = horzcat(responseData.vvecsFalling{:});
tvecsCombined = horzcat(responseData.tvecsCombined{:});
vvecsCombined = horzcat(responseData.vvecsCombined{:});

[fitResultsRising, cfitRising] = ...
    fit_pulse_response(tvecsRising, vvecsRising, pulseWidth);
[fitResultsFalling, cfitFalling] = ...
    fit_pulse_response(tvecsFalling, vvecsFalling, pulseWidth);
[fitResultsCombined, cfitCombined] = ...
    fit_pulse_response(tvecsCombined, vvecsCombined, pulseWidth);

eqRising = fitResultsRising.eqnShortPulseResponse;
eqFalling = fitResultsFalling.eqnShortPulseResponse;
eqCombined = fitResultsCombined.eqnShortPulseResponse;

plot_and_compare(tvecsRising, vvecsRising, cfitRising, eqRising);
plot_and_compare(tvecsFalling, vvecsFalling, cfitFalling, eqFalling);
plot_and_compare(tvecsCombined, vvecsCombined, cfitCombined, eqCombined);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_and_compare (tvecs, vvecs, cfitObject, equationStr)

% Plot the vectors
figure;
hold on
plot(tvecs, vvecs, 'Color', 'k', 'LineWidth', 2);
plot(cfitObject, 'r--');
text(0.1, 0.9, equationStr, 'Units', 'normalized');
hold off

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

[paramsRising, cfitRising] = fit_2exp(tvecsRising, vvecsRising);
[paramsFalling, cfitFalling] = fit_2exp(tvecsFalling, vvecsFalling);
[paramsCombined, cfitCombined] = fit_2exp(tvecsCombined, vvecsCombined);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%