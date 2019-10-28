function [stimcellIDs, spacing] = m3ha_network_define_actmode (actmode, actcellID, ncells, RERErad, RETCrad)
%% Determine what cells are stimulated in each activation mode
% 1 - Activate a single RE cell by injecting a train of current pulses
% 2 - Activate every (RERErad + 1)th RE cell by injecting trains of current pulses
% 3 - Activate 3 RE cells by injecting trains of current pulses
% 4 - Activate a single RE cell by changing the membrane potential instantaneously
% 5 - Activate RE cells with a Gaussian likelihood by changing the mp instantaneously
% 6 - Activate every 3rd RE cell by injecting trains of current pulses
% 7 - Activate all RE cells by injecting trains of current pulses
%
% Used by:
%        cd/m3ha_network_launch.m
%        cd/m3ha_network_raster_plot.m

% File History:
% 2017-11-06 Moved from m3ha_launch4.m
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch actmode
case {2, 3}
    spacing = RERErad + 1;
case 8
    spacing = RETCrad;
otherwise
    spacing = 1;
end

switch actmode
case 1
    stimcellIDs = actcellID;
case 2
    ct = 0;
    for c = 0:ncells-1
        if mod(c - actcellID, spacing) == 0
            ct = ct + 1;
            stimcellIDs(ct) = c;
        end
    end
case {3, 8}
    stimcellIDs = [actcellID, actcellID + spacing, actcellID - spacing];
case 6
    ct = 0;
    for c = 0:ncells-1
        if mod(c - actcellID, 3) == 0
            ct = ct + 1;
            stimcellIDs(ct) = c;
        end
    end
case 7
    stimcellIDs = 0:ncells-1;
case 9
    stimcellIDs = actcellID-5:actcellID+4;
case 10
    stimcellIDs = actcellID-10:actcellID+9;
otherwise
    error('actmode undefined');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%