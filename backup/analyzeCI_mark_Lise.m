% function [bigdata] = analyzeCI_mark_Lise(alldata, si, timeC)

%%
ccc
% save('/home/mark/matlab_temp_variables/analyzeCI_mark')
load('/home/mark/matlab_temp_variables/analyzeCI_mark')

%% adam's parameters
nsweeps = 1 ; % adam's code originally parsed out sweeps from abf file.  I'm just sending single sweeps into this function so nsweeps always =1 
tw = 100;       % Time window for computing frequencies in ms 

%% detection code
ntps = size(alldata,1);         % Number of time points
si_insec = si*10^-6;            % Sampling internal in seconds
                
% Find vector for timepoints in msec
tps = ( si/1000 : si/1000 : ntps * si/1000 )';

% Find spikes
is_local_maximum = zeros(size(alldata));
is_local_minimum = zeros(size(alldata));
is_spike = zeros(size(alldata));

for i = 1:nsweeps
% Finds all local maxima and minima
    for j = 2:ntps-1
        is_local_maximum(j,1,i) = alldata(j,1,i) > alldata(j-1,1,i) && alldata(j,1,i) >= alldata(j+1,1,i);
        is_local_minimum(j,1,i) = alldata(j,1,i) < alldata(j-1,1,i) && alldata(j,1,i) <= alldata(j+1,1,i);
    end

% Finds all spike peaks 
% Criteria for a spike: 
%  (1) Must be a local maximum 10 mV higher than the previous local minimum
    for j = 2:ntps-1
        if is_local_maximum(j,1,i)
            plmin = j-1; % possible index of previous local minimum
            while ~ (is_local_minimum(plmin,1,i) || plmin == 1) 
                plmin = plmin - 1;
            end

            if plmin > 1
            % Compare rise in membrane potential to thresholds (10 mV)
                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 10;
            end
        end
    end
                    %  (2) Must be 5 mV higher than the minimum value between the spike and the following spike
    for j = 2:ntps-1
        if is_spike(j,1,i)
            fspike = j+1; % possible index of following spike
            while ~ (is_spike(fspike,1,i) || fspike == ntps) 
                fspike = fspike + 1;
            end
            is_spike(j,1,i) = alldata(j,1,i) - min(alldata(j:fspike,1,i)) > 5;
        end
    end
end

%% CALCULATE SPIKE PROPS

% Variables for all computations
ntw = tw * 1000 / si;       % Number of timepoints within the time window
ntw_h = floor(ntw/2);   % Half of above
tw_insec = tw/1000;             % Time window in seconds
                
% Compute frequencies with a rectangular time window (NOT YET NOTMALIZED AT THE ENDS)
allfreq = zeros(size(alldata));
for j = 1:ntps
    left = max(1, j-ntw_h);
    right = min(ntps, j+ntw-ntw_h-1);
    allfreq(j,1,:) = sum(is_spike(left:right,1,:),1)/tw_insec;
end

% Compute frequencies with an approximate Gaussian time window (NOT YET NOTMALIZED AT THE ENDS)
allfreq_AG = zeros(size(alldata));
AG_window_un = normpdf([1:1:ntw]',ntw_h,floor(ntw/4)); % Approximate Gaussian time window, unnormalized
AG_window = AG_window_un/(sum(AG_window_un)*si_insec); % Approximate Gaussian time window, normalized
for j = 1 : ntw_h
    left = 1;
    right = j + (ntw-ntw_h-1);
    AG_window_shortened = AG_window(ntw_h+2-j:ntw) / (sum(AG_window(ntw_h+2-j:ntw))*si_insec) ;
    allfreq_AG(j,1,:) = AG_window_shortened' * squeeze(is_spike(left:right,1,:));
end

for j = ntw_h + 1 : ntps - (ntw-ntw_h-1)
    left = j - ntw_h;
    right = j + (ntw-ntw_h-1);
    allfreq_AG(j,1,:) = AG_window' * squeeze(is_spike(left:right,1,:));
end

for j = ntps - (ntw-ntw_h-1) + 1 : ntps
    left = j - ntw_h;
    right = ntps;
    AG_window_shortened = AG_window(1:ntps-(j-ntw_h)+1) / (sum(AG_window(1:ntps-(j-ntw_h)+1))*si_insec) ;
    allfreq_AG(j,1,:) = AG_window_shortened' * squeeze(is_spike(left:right,1,:));
end

%% collect data
bigdata.numberSpikes = size(find(is_spike==1),1) ;
bigdata.spikeTimes = timeC(find(is_spike==1)) ; % in sec
bigdata.spikeIDX = find(is_spike==1) ; 
bigdata.detectedSpikes = is_spike ;
bigdata.spikeFreq.reg = allfreq ;
bigdata.spikeFreq.AG = allfreq_AG ;

%% verify detection
plot(timeC, alldata, 'k')
hold on
plot(timeC(bigdata.spikeIDX), alldata(bigdata.spikeIDX), 'ro')
% pause(2)
% close all
