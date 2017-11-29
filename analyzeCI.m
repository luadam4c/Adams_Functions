% function [alldata] = analyzeCI(date)	
date = '20140729';

% Usage: analyzeCI(date)
%        date is in the format 'YYYYMMDD'

% File history
% 20140730 - Adapted from analyzeCI_20140730
% 20160615 - Fixed typo

% Parameters to set
max_celln = 10; % Maximum number of cells in the data
max_protn = 4;  % Maximum number of types of CI protocols in the data
max_recn = 21;  % Maximum number of recordings for each protocol type in the data
tw = 100;       % Time window for computing frequencies in ms 

% Set input file
prot = 'CI';

% Run through all possible files

for celln = 1:1:max_celln
    for protn = 0:1:max_protn
        for recn = 1:1:max_recn
% for celln = 4
%     for protn = 1
%         for recn = 4
            file = sprintf('PClamp_Data/%s%02d_%s%d_%d.abf', date, celln, prot, protn, recn);
            if exist(file, 'file') == 2 % if file exists

                % Read file to alldata
                [alldata, si] = abf2load(file); % si is the sampling interval in usec

                % Find data parameters
                nsweeps = size(alldata, 3);     % Number of sweeps
                ntps = size(alldata,1);         % Number of time points
                si_insec = si*10^-6;            % Sampling interval in seconds
                
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
                %             flmin = j+1; % possible index of following local minimum
                %             while ~ (is_local_minimum(flmin,1,i) || flmin == ntimepoints) 
                %                 flmin = flmin + 1;
                %             end
                %             flmax = j+1; % possible index of following local maximum
                %             while ~ (is_local_maximum(flmax,1,i) || flmax == ntimepoints) 
                %                 flmax = flmax + 1;
                %             end
                %             if plmin > 1 && flmin < ntimepoints
                            if plmin > 1
                                % Compare rise in membrane potential to thresholds (10 mV)
                                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 10;
                %                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 20;
                %                is_spike(j,1,i) = alldata(j,1,i) - alldata(plmin,1,i) > 10 && ...
                %                     (alldata(j,1,i) - alldata(flmin,1,i) > 10 || alldata(flmax,1,i) - alldata(j,1,i) < 10);
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

%                 % Compute frequencies with a Gaussian time window (NOT YET NORMALIZED PROPERLY)
%                 % WARNING: Very slow!
%                 allfreq_G = zeros(size(alldata));
%                 sigma = ntps_tw/4;          % Standard deviation is 1/4 the time window, so that 95.5% of utilized data is within the time window
%                 for j = 1:ntps
%                     G_window_un = normpdf(tps,tps(j),sigma); % Gaussian time window unnormalized
%                     G_window = G_window_un/sum(G_window_un); % Gaussian time window normalized
%                     for i = 1:nsweeps
%                         allfreq_G(j,1,i) = G_window'*is_spike(:,1,i);
%                     end
%                 end                
                
                % Plot raw data with spikes (each sweep individually)
                for i = 1:nsweeps
                    cdata = alldata(:,1,i);
                    sp_ind = find(is_spike(:,1,i));
                    figure(i)
                    plot(tps, cdata, 'k')
                    hold on
                    plot(tps(sp_ind), cdata(sp_ind), 'xr')
                    % axis([0 10000 -160 40])
                    axis([0 4000 -160 40])
                    % xlim([0 4000])
                    title(sprintf('Data for %s%02d_%s%d_%d.abf, Sweep #%d', date, celln, prot, protn, recn, i), 'interpreter', 'none')
                    xlabel('Time (ms)')
                    ylabel('Membrane Potential (mV)')
                    saveas(gcf, sprintf('PClamp_Data/%s%02d_%s%d_%d_sweep%d', date, celln, prot, protn, recn, i), 'png')
                    hold off
                end              
                
% WARNING: The following code crashes Matlab for some reason (but it doesn't in debug mode)
%                 % Plot raw data with spikes & frequency (each sweep individually)
%                 for i = 1:nsweeps
%                     cdata = alldata(:,1,i);
%                     sp_ind = find(is_spike(:,1,i));
%                     padding = zeros(ntps-length(sp_ind),1);
%                     cfreq = allfreq(:,1,i);
%                     figure(i)
%                     [haxes, hline1, hline2] = plotyy([tps, [tps(sp_ind); padding]], [cdata, [cdata(sp_ind); padding]], tps, cfreq);
%                     set(hline1(1,1),'Color','k')
%                     set(hline1(2,1),'LineStyle','none','Marker','x','Color','r')
%                     set(hline2,'Color','b')
%                     xlabel(haxes(1),'Time (ms)','Color','k')
%                     ylabel(haxes(1),'Membrane Potential (mV)','Color','k')
%                     ylabel(haxes(2),'Frequency (Hz)','Color','b')
%                     % axis([0 10000 -160 40])
%                     axis(haxes(1), [0 4000 -160 40])
%                     axis(haxes(2), [0 4000 0 300])
%                     % xlim([0 4000])
%                     title(sprintf('Data for %s%02d_%s%d_%d.abf, Sweep #%d. Time window = %d ms', date, celln, prot, protn, recn, i, tw), 'interpreter', 'none')
%                     saveas(gcf, sprintf('PClamp_Data/%s%02d_%s%d_%d_sweep%d', date, celln, prot, protn, recn, i), 'png')
%                     hold off
%                 end

                % Variables for all plots
                jmap = colormap(jet);

                % Plot raw data (all sweeps together)
                figure(nsweeps + 1)
                for i = 1:nsweeps
                    cdata = alldata(:,1,i);
                    plot(tps, cdata, 'color', jmap((i * floor(size(jmap,1)/10)), :))
                    hold on;
                end
                hold off
                %axis([0 10000 -160 40])
                axis([0 4000 -160 40])
                title(sprintf('Data for %s%02d_%s%d_%d.abf', date, celln, prot, protn, recn), 'interpreter', 'none')
                xlabel('Time (ms)')
                ylabel('Membrane Potential (mV)')
                % HOW TO GENERALIZE THE FOLLOWING?
                legend1 = legend('Sweep #1','Sweep #2','Sweep #3','Sweep #4','Sweep #5','Sweep #6','Sweep #7','Sweep #8','Sweep #9','Sweep #10');
                set(legend1,'Position',[0.132589285714285 0.419345238095237 0.209821428571428 0.501488095238095]);
                saveas(gcf, sprintf('PClamp_Data/%s%02d_%s%d_%d_all_spikes', date, celln, prot, protn, recn), 'png')

                % Plot frequencies (all sweeps together, rectangular time window)
                figure(nsweeps + 2)
                for i = 1:nsweeps
                    cfreq = allfreq(:,1,i);
                    plot(tps, cfreq, 'color', jmap((i * floor(size(jmap,1)/10)), :))
                    hold on;
                end
                hold off
                %axis([0 10000 -160 40])
                axis([0 4000 0 120])
                title(sprintf('Frequency Data for %s%02d_%s%d_%d.abf, Rectangular time window width = %d ms', date, celln, prot, protn, recn, tw), 'interpreter', 'none')
                xlabel('Time (ms)')
                ylabel('Firing rate (Hz)')
                legend1 = legend('Sweep #1','Sweep #2','Sweep #3','Sweep #4','Sweep #5','Sweep #6','Sweep #7','Sweep #8','Sweep #9','Sweep #10');
                set(legend1,'Position',[0.132589285714285 0.419345238095237 0.209821428571428 0.501488095238095]);
                saveas(gcf, sprintf('PClamp_Data/%s%02d_%s%d_%d_all_frequencies', date, celln, prot, protn, recn), 'png')

                % Plot frequencies (all sweeps together, Approximate Gaussian time window)
                figure(nsweeps + 3)
                for i = 1:nsweeps
                    cfreq = allfreq_AG(:,1,i);
                    plot(tps, cfreq, 'color', jmap((i * floor(size(jmap,1)/10)), :))
                    hold on;
                end
                hold off
                %axis([0 10000 -160 40])
                axis([0 4000 0 120])
                title(sprintf('Frequency Data for %s%02d_%s%d_%d.abf, Approximate Gaussian time window sigma = %g ms', date, celln, prot, protn, recn, tw/4), 'interpreter', 'none')
                xlabel('Time (ms)')
                ylabel('Firing rate (Hz)')
                legend1 = legend('Sweep #1','Sweep #2','Sweep #3','Sweep #4','Sweep #5','Sweep #6','Sweep #7','Sweep #8','Sweep #9','Sweep #10');
                set(legend1,'Position',[0.132589285714285 0.419345238095237 0.209821428571428 0.501488095238095]);
                saveas(gcf, sprintf('PClamp_Data/%s%02d_%s%d_%d_all_frequencies_AG', date, celln, prot, protn, recn), 'png')              
                
                % For debug
                c_is_local_maximum = is_local_maximum(:,:,1);
                c_is_local_minimum = is_local_minimum(:,:,1);
                c_is_spike = is_spike(:,:,1);
                c_data = alldata(:,1,1);

            end
        end
    end
end
