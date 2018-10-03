function [params_final, cpa, pw, dv_rec, rmse_R, rmse_F, params_L_F2, params_S_R2, params_L_F1, params_S_R1] = find_passive_params (tvec0, ivec0s, vvec0s, cpwin, cprwin, cpmid, plotflag, outfolder, filebase, ivec1s, fitmode)
%% Extract passive parameters from both the rising and falling phase of the current pulse response
% Usage: [params_final, cpa, pw, dv_rec, rmse_R, rmse_F, params_L_F2, params_S_R2, params_L_F1, params_S_R1] = find_passive_params (tvec0, ivec0s, vvec0s, cpwin, cprwin, cpmid, plotflag, outfolder, filebase, ivec1s, fitmode)
% Outputs:
%       params        - structure of parameters inferred from fits, 
%                    estimated either from
%                    (1) the extrapolated coefficients of a long pulse response from fits to the falling phase
%                    or (2) the coefficients of the short pulse response from fits to the rising phase
%                    includes:
%                 C0        - Coefficient of first exponential, C0 (mV)
%                 tau0        - Time constant of first exponential, tau0 (ms)
%                 C1        - Coefficient of second exponential, C1 (mV)
%                 tau1        - Time constant of second exponential, tau1 (ms)
%                 cpa_mean    - Average current pulse amplitude (pA)
%                dvss        - Voltage difference at steady state (mV)
%                Rinput        - Input resistance of setup (MOhm)
%                alpha1        - alpha1 = sqrt(tau0/tau1 - 1)
%                L_init        - Initial guess for electrotonic length, L_init
%                L        - Electrotonic length, L
%                rho        - Dendritic to somatic conductance ratio, rho
%                Rmemb        - Input resistance of cell, R_N (Mohm)
%                Rsoma        - Somatic resistance, R_soma (Mohm)
%                Rdend        - Dendritic resistance, R_dend (Mohm)
%                Gmemb        - Input conductance of cell, G_N (uS)
%                Gsoma        - Somatic conductance, G_soma (uS)
%                Gdend        - Dendritic conductance, G_dend (uS)
%                taum        - Membrane time constant, taum (ms)
%                Rm        - Specific membrane resistivity (Ohm-cm^2)
%                rad_soma    - Radius of soma (um)
%                diam_dend    - Diameter of dendrite (um)
%                rad_dend    - Radius of dendrite (um)
%                length_dend    - Length of dendrite (um)
%       cpa         - vector of current pulse amplitudes (pA) for each sweep
%       pw          - vector of pulse width (ms) for each sweep
%       dv_rec      - vector of overall change in membrane potential recorded (mV) for each sweep
%       rmse_R      - vector of root-mean-squared errors (mV) in the rising phase for each sweep
%       rmse_F      - vector of root-mean-squared errors (mV) in the falling phase for each sweep
% Arguments:
%       tvec0       - original time vector, must be a column vector in ms
%       ivec0s      - original current vector(s) in pA
%                   must be a numeric array with columns the same length as tvec0
%       vvec0s      - original voltage vector(s) in mV
%                   must be a numeric array with columns the same length as tvec0
%        cpwin        - (opt) Window in which current pulse would lie, e.g. [95, 115]
%                must be within range of [tbase tfinal]
%                default == [tbase tvec0(end)]
%        cprwin        - (opt) Window in which current pulse response would lie, e.g. [95, 500]
%                must be within range of [tbase tfinal]
%                default == [tbase tvec0(end)]
%        cpmid        - (opt) approximate midpoint of the current pulse (ms)
%                must be within range of cprwin
%                default == 0.5*cpwin(1) + 0.5*cpwin(2)
%        plotflag    - (opt) whether to plot current trace or not
%                must be 0 or 1
%                default == 0
%        outfolder     - (opt) directory to place outputs, e.g. '/media/adamX/m3ha/data_dclamp/take4/test_sweeps/'
%                must be a directory
%                default == pwd
%        filebase     - (opt) base of filename (without extension), e.g. 'A100110_0008_R18'
%                must be a char array
%                default == 'unnamed'
%        ivec1s        - (opt) median-filtered current vector(s) in pA
%                must be a numeric array with columns the same length as tvec0
%                default == medfilt1(ivec0s, round(mfw2/sims))
%        fitmode        - (opt) 0 - all data
%                    1 - all of g incr = 100%, 200%, 400%
%                    2 - all of g incr = 100%, 200%, 400% 
%                        but exclude cell-pharm-g_incr sets containing problematic sweeps
%                must be one of the above values
%                default == does not exist
%
% Requires:    
%        /media/adamX/m3ha/data_dclamp/specs_for_fitmode (if fitmode is used)
%        /home/Matlab/Downloaded_Functions/subplotsqueeze.m
%        /home/Matlab/Adams_Functions/check_subdir.m
%        /home/Matlab/Adams_Functions/print_structure.m
%
% Used by:
%        ~/m3ha/data_dclamp/dclampDataExtractor.m
%        ~/m3ha/data_dclamp/dclampPassiveFitter.m
%
%
% 2016-10-27 - Adapted from find_IPSC_peak.m and dclampDataExtractor.m
% 2016-10-31 - Changed equation form of ft to 'a*(exp(-x/b)-1)+c*(exp(-x/d)-1)' from 'a*exp(-x/b)+c*exp(-x/d)+e'
% 2016-10-31 - Removed '- round(0.5/sims)' from the definition of base_ind
% 2016-10-31 - Made sure tau0 >= tau1
% 2016-11-01 - Added fitmode so that the titles and filenames can change
% 2016-11-01 - Exclude faulty traces, including those with spontaneous spikes, from the fitting
% 2016-11-01 - Added pulse width
% 2016-11-01 - Added long pulse response
% 2016-11-01 - Added L, rho, taum
% 2016-11-01 - Plot ivec1s only (ivec0s is not directly used in the analyses)
% 2016-11-02 - Moved check directories to check_subdir.m
% 2016-11-02 - Added the falling phase of the current pulse response
% 2016-11-12 - Now returns estimates from pooled data if those from averaged data don't exist
% 2016-11-12 - Now outputs all estimates
% 2016-12-04 - Added series resistance Rs and changed the way somatic and dendritic resistances are computed
% 2016-12-04 - Increase the upper bound of tau to 500 ms
% 2016-12-04 - Added rmse_R && rmse_F, the root-mean-squared errors of the rising and falling phase, respectively, for each sweep
% 2016-12-04 - Added the functions fit_setup && measure_error
% 2016-12-04 - Fixed measure_error so that sweeps yielding nonsensical responses have rmse = Inf
% 2016-12-05 - Added tau0_range = [20, 500] and tau1_range = [0, 20]
% 2016-12-05 - Removed typtau and tau_max, changed initial conditions to tau0_range(1) and tau1_range(2)
% 2016-12-05 - Added tau0_range_R and tau1_range_R
% 2016-12-05 - Added typtau0, typtau1, typtau0_R, typtau1_R
% 2017-12-21 - SpecsForFitmode() -> specs_for_fitmode()
% 2018-01-24 - Added isdeployed
% 2018-10-03 - Changed tabs to spaces
%  

%% Flags
printfieldsflag = 0;

%% Parameters used for data analysis
mfw2 = 10;                    % width in ms for the median filter for corrupted data (current traces)
mvw = 0.5;                    % width in ms for calculating mean voltage for input resistance calculations
maxScalingFactor = 10;        % maximum scaling factor from initial voltage amplitude value
typtau0 = 50;        % typical tau0 ~ 50 ms
tau0_range = [20, 500];    % range of tau0
typtau1 = 6;        % typical tau1 ~ 6 ms
tau1_range = [0, 20];    % range of tau1
typtau0_R = 23;        % typical tau0 ~ 50 ms
tau0_range_R = [7, 200];% range of tau0 in rising phase, make different for now
typtau1_R = 0.8;    % typical tau1 ~ 6 ms
tau1_range_R = [0, 7];    % range of tau1 in rising phase, make different for now
sp_thr_init = -45;    % Initial amplitude threshold in mV for detecting a spike 
eqform_R = 'a*(exp(-x/b)-1)+c*(exp(-x/d)-1)';    % double exponential, rising phase
eqform_F = '-a*exp(-x/b)-c*exp(-x/d)';        % double exponential, falling phase

%% Fixed parameters for all cells
Cm = 0.88;        % specific membrane capacitance (uF/cm^2)
Ra = 173;        % axial resistivity (Ohm-cm)
Rs = 10;        % series resistance (MOhm)

%% Directories for placing figures
directories = {'/cprRising/', '/cprFalling/', '/passive/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
    error('Not enough input arguments, type ''help find_passive_params'' for usage');
elseif isempty(tvec0) || isempty(ivec0s) || isempty(vvec0s)
    error('First three inputs cannot be empty!');
elseif ~isnumeric(tvec0) || ~isnumeric(ivec0s) || ~isnumeric(vvec0s)
    error('First three inputs must be numbers or numeric arrays!');
elseif ~isequal(size(tvec0, 1), size(ivec0s, 1))
    error('Time and current vectors do not have the same length!');
elseif ~isequal(size(tvec0, 1), size(vvec0s, 1))
    error('Time and voltage vectors do not have the same length!');
elseif size(tvec0, 2) > 1
    error('tvec0 must be a column vector!');
elseif nargin >= 4 && length(cpwin) < 2
    error('cpwin must have a start and an end!');
elseif nargin >= 5 && length(cprwin) < 2
    error('cprwin must have a start and an end!');
end

%% Extract info from data
sims = tvec0(2) - tvec0(1);                    % sampling interval in ms
ndps = length(tvec0);
nswps = size(ivec0s, 2);
tbase = tvec0(1) - sims;
tfinal = tvec0(end) + sims;                    % in case sims is not divisible by cprwin(2)

%% Check more arguments
if nargin >= 4 && cpwin(1) < tbase
    error('cpwin(1) == %g < tbase == %g is out of range!', cpwin(1), tbase);
elseif nargin >= 4 && cpwin(2) > tfinal
    error('cpwin(2) == %g > tfinal == %g is out of range!', cpwin(2), tfinal);
elseif nargin >= 5 && cprwin(1) < tbase
    error('cprwin(1) == %g < tbase == %g is out of range!', cprwin(1), tbase);
elseif nargin >= 5 && cprwin(2) > tfinal
    error('cprwin(2) == %g > tfinal == %g is out of range!', cprwin(2), tfinal);
elseif nargin >= 6 && (cpmid < cprwin(1) || cpmid > cprwin(2))
    error('cpmid is out of range!');
elseif nargin >= 7 && (plotflag ~= 1 && plotflag ~= 0 && plotflag ~= false && plotflag ~= true)
    error('plotflag is out of range!');
elseif nargin >= 8 && ~ischar(outfolder)
    error('outfolder must be a character array!');
elseif nargin >= 9 && ~ischar(filebase)
    error('filebase must be a character array!');
elseif nargin >= 10 && ~isnumeric(ivec1s)
    error('ivec1s must be a numeric array!');
elseif nargin >= 10 && ~isequal(size(ivec1s, 1), size(tvec0, 1))
    error('ivec1s must have the same column length as tvec0!');
elseif nargin >= 11 && (~isnumeric(fitmode) || ~(fitmode == 0 || fitmode == 1 || fitmode == 2))
    error('fitmode out of range!');
end

%% Set defaults for optional arguments
if nargin < 4
    cpwin = [tbase tvec0(end)];    
end
if nargin < 5
    cprwin = [tbase tvec0(end)];    
end
if nargin < 6
    cpmid = 0.5*cpwin(1) + 0.5*cpwin(2);
end
if nargin < 7
    plotflag = 0;
end
if nargin < 8
    outfolder = pwd;
end
if nargin < 9
    filebase = 'unnamed';
end
if nargin < 10
    ivec1s = zeros(ndps, nswps);        % preallocate
%    for swp = 1:nswps                % FOR each sweep
    parfor swp = 1:nswps                % FOR each sweep
        % Median filter current traces to get rid of corrupted data
        ivec1s(:, swp) = medfilt1(ivec0s(:, swp), round(mfw2/sims));
    end
end

%% Only do the following if fitmode is used
if nargin >= 11
    % Set suffices and title modifications for each fitmode
    [suffix, title_mod] = specs_for_fitmode (fitmode);

    % Change directory names to make it specific to each fitmode
    for k = 1:numel(directories)
        temp = directories{k}(2:end);        % remove the leading '/'
        directories{k} = [ '/', strrep(temp, '/', [suffix, '/']) ];
    end
else
    suffix = '';
    title_mod = '';
end

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                    % for subplotsqueeze.m
end

% Display standard output header
% fprintf('ANALYZING passive parameters for %s ...\n', filebase);
% fprintf('Sampling interval == %g ms\n', sims);

%% Find current pulses for each sweep
first_dip_pt = cell(1, nswps);            % index of current pulse start
cpa = zeros(1, nswps);                % current pulse amplitude (pA)
pw = zeros(1, nswps);                % pulse width (ms)
dv_rec = zeros(1, nswps);            % overall change in membrane potential recorded (mV)
xvecRs = cell(1, nswps);            % adjusted rising phase time vectors (start at 0)
yvecRs = cell(1, nswps);            % adjusted rising phase voltage vectors (start at around 0)
xvecFs = cell(1, nswps);            % adjusted falling phase time vectors (start at 0)
yvecFs = cell(1, nswps);            % adjusted falling phase voltage vectors (asymptote to around 0)
used = zeros(1, nswps);                % whether the sweep is used for fitting
ind_begin = find(tvec0 >= cprwin(1), 1);
ind_mid = find(tvec0 <= cpmid, 1, 'last');
ind_end = find(tvec0 <= cprwin(2), 1, 'last');
ind = ind_begin:ind_end;            % indices of interest
ind1 = ind_begin:ind_mid;            % indices for finding current pulse start
ind2 = ind_mid:ind_end;                % indices for finding current pulse end
mvind = 1:round(mvw/sims);            % indices for taking the mean of voltages
base_ind = zeros(nswps, length(mvind));
last_ind = zeros(nswps, length(mvind));
ivec1s_part1 = ivec1s(ind1, :);            % Use median-filtered current trace
ivec1s_part2 = ivec1s(ind2, :);
ivec1s_part1_begin = ind1(1);
ivec1s_part2_begin = ind2(1);
ivec1s_part3 = ivec1s(ind_begin - 1 + mvind, :);% indices for measuring cp baseline
ivec1s_part4 = ivec1s(ind_mid - 1 + mvind, :);    % indices for measuring cp peak
% for swp = 1:nswps                % FOR each sweep
parfor swp = 1:nswps                % FOR each sweep
    cpa(swp) = mean(ivec1s_part4(:, swp)) - mean(ivec1s_part3(:, swp));    % current pulse amplitude (pA)
    first_dip_pt{swp} = find(ivec1s_part1(:, swp) > cpa(swp)/4, 1, 'last');
    if isempty(first_dip_pt{swp})                % exclude faulty traces
        xvecRs{swp} = [];
        yvecRs{swp} = [];
        xvecFs{swp} = [];
        yvecFs{swp} = [];
        used(swp) = false;
    else
        % Find the voltage trace corresponding to the current pulse response
        cpstart = (ivec1s_part1_begin - 1) + first_dip_pt{swp};    % index of current pulse start
        before_rise_pt = find(ivec1s_part2(:, swp) < cpa(swp) * 3/4, 1, 'last');
        cpend = (ivec1s_part2_begin - 1) + before_rise_pt;    % index of current pulse end
        pw(swp) = (cpend - cpstart) * sims;             % pulse width (ms)
        base_ind(swp, :) = cpstart - fliplr(mvind);        % base indices
        last_ind(swp, :) = cpend - fliplr(mvind);        % last indices of the current pulse
        basev = mean(vvec0s(base_ind(swp, :), swp));
        lastv = mean(vvec0s(last_ind(swp, :), swp));
        dv_rec(swp) = basev - lastv;         % change in membrane potential on voltage trace (mV)
        vmax = max(vvec0s(ind, swp));        % maximum voltage in current pulse response window (mV)
        if dv_rec(swp) <= 0 || cpa(swp) >= 0 || pw(swp) <= 0 ...
            || vmax > sp_thr_init        % exclude faulty traces, including those with spontaneous spikes
            xvecRs{swp} = [];
            yvecRs{swp} = [];
            xvecFs{swp} = [];
            yvecFs{swp} = [];
            used(swp) = false;
        else
            xvecRs{swp} = tvec0(cpstart:cpend) - tvec0(cpstart);
            yvecRs{swp} = vvec0s(cpstart:cpend, swp) - basev;
            xvecFs{swp} = tvec0(cpend:ind_end) - tvec0(cpend);
            yvecFs{swp} = vvec0s(cpend:ind_end, swp) - basev;
            used(swp) = true;
        end
    end
end
if ~isempty(find(dv_rec > 0, 1))
    dv_rec_mean = mean(dv_rec(dv_rec > 0));        % average recorded voltage change (mV), excluding faulty traces
else
    dv_rec_mean = NaN;
end
if ~isempty(find(cpa < 0, 1))
    cpa_mean = mean(cpa(cpa < 0));            % average current pulse amplitude (pA), excluding faulty traces
else
    cpa_mean = NaN;
end
if ~isempty(find(pw > 0, 1))
    pw_mean = mean(pw(pw > 0));            % average pulse width (ms), excluding faulty traces
else
    pw_mean = NaN;
end
fprintf('Average recorded voltage change (mV) is %g\n', dv_rec_mean);
fprintf('Average current pulse amplitude (pA) is %g\n', cpa_mean);
fprintf('Average pulse width (ms) is %g\n', pw_mean);
fprintf('\n');

%% Find a goodness-of-fit measure (root-mean-squared error) for each sweep
rmse_R = measure_error(xvecRs, yvecRs, eqform_R, dv_rec_mean, maxScalingFactor, typtau0_R, typtau1_R, tau0_range_R, tau1_range_R, pw_mean);
rmse_F = measure_error(xvecFs, yvecFs, eqform_F, dv_rec_mean, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range, pw_mean);

%% Put data from all current pulse responses together in two different ways
% Method 1: Put all data points together, excluding faulty traces
xvecR_all = [];    
yvecR_all = [];
xvecF_all = [];    
yvecF_all = [];
for swp = 1:nswps
    xvecR_all = cat(1, xvecR_all, xvecRs{swp});
    yvecR_all = cat(1, yvecR_all, yvecRs{swp});
    xvecF_all = cat(1, xvecF_all, xvecFs{swp});
    yvecF_all = cat(1, yvecF_all, yvecFs{swp});
end

% Method 2: Average all current pulse responses
xvecRs_cl = xvecRs(~cellfun(@isempty, xvecRs));        % Remove all elements that are empty
yvecRs_cl = yvecRs(~cellfun(@isempty, yvecRs));        % ditto
xvecFs_cl = xvecFs(~cellfun(@isempty, xvecFs));        % ditto
yvecFs_cl = yvecFs(~cellfun(@isempty, yvecFs));        % ditto
minlen_R = min(cellfun(@length, xvecRs_cl));        % minimum length of the current pulse responses, excluding faulty traces
minlen_F = min(cellfun(@length, xvecFs_cl));        % ditto for falling phase
if minlen_R > 0
    xvecR_mean = xvecRs_cl{1}(1:minlen_R);            % Restrict adjusted time vector to this minimum length
    for swp = 1:numel(yvecRs_cl)
        yvecRs_cl{swp} = yvecRs_cl{swp}(1:minlen_R);    % Restrict adjusted voltage vectors to this minimum length
    end
    yvecR_mean = mean(cell2mat(yvecRs_cl), 2);        % Calculate the mean of all current pulse responses
else
    xvecR_mean = [];
    yvecR_mean = [];
end
if minlen_F > 0
    xvecF_mean = xvecFs_cl{1}(1:minlen_F);            % ditto for falling phase
    for swp = 1:numel(yvecRs_cl)
        yvecFs_cl{swp} = yvecFs_cl{swp}(1:minlen_F);    % ditto for falling phase
    end
    yvecF_mean = mean(cell2mat(yvecFs_cl), 2);        % ditto for falling phase
else
    xvecF_mean = [];
    yvecF_mean = [];
end

%% Fit voltage trace with two exponentials
if ~isempty(xvecR_all)
    % Fit pooled rising phase data with double exponential
    [cfit_R1, eqnS_R1, eqnL_R1, CS0_R1, tau0_R1, CS1_R1, tau1_R1, ~, ~] = ...
        fit_data ('rising', eqform_R, xvecR_all, yvecR_all, dv_rec_mean, maxScalingFactor, typtau0_R, typtau1_R, tau0_range_R, tau1_range_R, pw_mean);

    % Estimate parameters from short-pulse coefficients
    [params_S_R1] = coeff2params (CS0_R1, tau0_R1, CS1_R1, tau1_R1, cpa_mean, Cm, Ra, Rs);
    if printfieldsflag
        print_structure(params_S_R1);
    end
else
    params_S_R1 = NaN;
end
if ~isempty(xvecF_all)
    % Fit pooled falling phase data with double exponential
    [cfit_F1, eqnS_F1, eqnL_F1, ~, tau0_F1, ~, tau1_F1, CL0_F1, CL1_F1] = ...
        fit_data ('falling', eqform_F, xvecF_all, yvecF_all, dv_rec_mean, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range, pw_mean);

    % Estimate parameters from long-pulse coefficients
    [params_L_F1] = coeff2params (CL0_F1, tau0_F1, CL1_F1, tau1_F1, cpa_mean, Cm, Ra, Rs);
    if printfieldsflag
        print_structure(params_L_F1);
    end
else
    params_L_F1 = NaN;
end
if ~isempty(xvecR_mean)
    % Fit averaged rising phase data with double exponential
    [cfit_R2, eqnS_R2, eqnL_R2, CS0_R2, tau0_R2, CS1_R2, tau1_R2, ~, ~] = ...
        fit_data ('rising', eqform_R, xvecR_mean, yvecR_mean, dv_rec_mean, maxScalingFactor, typtau0_R, typtau1_R, tau0_range_R, tau1_range_R, pw_mean);

    % Estimate parameters from short-pulse coefficients
    [params_S_R2] = coeff2params (CS0_R2, tau0_R2, CS1_R2, tau1_R2, cpa_mean, Cm, Ra, Rs);
    if printfieldsflag
        print_structure(params_S_R2);
    end
else
    params_S_R2 = NaN;
end
if ~isempty(xvecF_all)
    % Fit averaged falling phase data with double exponential
    [cfit_F2, eqnS_F2, eqnL_F2, ~, tau0_F2, ~, tau1_F2, CL0_F2, CL1_F2] = ...
        fit_data ('falling', eqform_F, xvecF_all, yvecF_all, dv_rec_mean, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range, pw_mean);

    % Estimate parameters from long-pulse coefficients
    [params_L_F2] = coeff2params (CL0_F2, tau0_F2, CL1_F2, tau1_F2, cpa_mean, Cm, Ra, Rs);
    if printfieldsflag
        print_structure(params_L_F2);
    end
else
    params_L_F2 = NaN;
end

%% Decide on final values: use falling phase, averaged; if not exist, use rising phase, averaged; 
%    if not exist, use falling phase, pooled; if not exist, use rising phase, pooled
if exist('params_L_F2', 'var') == 1 && isstruct(params_L_F2)
    params_final = params_L_F2;
elseif exist('params_S_R2', 'var') == 1 && isstruct(params_S_R2)
    params_final = params_S_R2;
elseif exist('params_L_F1', 'var') == 1 && isstruct(params_L_F1)
    params_final = params_L_F1;
elseif exist('params_S_R1', 'var') == 1 && isstruct(params_S_R1)
    params_final = params_S_R1;
else
    params_final = NaN;
end    

%% Plot results
if plotflag == 1
    % Check if needed directories exist in outfolder
    check_subdir(outfolder, directories);

    % Plot rising phase of current pulse response
    h = figure(2000);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Rising phase of current pulse response');
    clf(h);
    subplot(2,2,1);
    plot_current_pulse(nswps, used, ind, tvec0, ivec1s, first_dip_pt, base_ind, last_ind, cpwin, cpa_mean);
    subplot(2,2,2);
    plot_voltage_response(nswps, used, ind, tvec0, vvec0s, first_dip_pt, base_ind, last_ind, cpwin)
    subplot(2,2,3);
    if ~isempty(xvecR_all)
        plot_cfit('rising', cpwin, cfit_R1, xvecR_all, yvecR_all, params_S_R1, ...
                eqnS_R1, eqnL_R1, tau0_R1, tau1_R1);
    end
    subplot(2,2,4);
    if ~isempty(xvecR_mean)
        plot_cfit('rising', cpwin, cfit_R2, xvecR_mean, yvecR_mean, params_S_R2, ...
                eqnS_R2, eqnL_R2, tau0_R2, tau1_R2);
    end
    figname = fullfile(outfolder, directories{1}, [filebase, '_cprRising', suffix, '.png']);
    suptitle(strjoin({'Rising phase of current pulse response for', strrep(filebase, '_', '\_'), title_mod}));
    saveas(h, figname);
    % close(h);

    % Plot falling phase of current pulse response
    h = figure(2001);
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Falling phase of current pulse response');
    clf(h);
    subplot(2,2,1);
    plot_current_pulse(nswps, used, ind, tvec0, ivec1s, first_dip_pt, base_ind, last_ind, cprwin, cpa_mean);
    subplot(2,2,2);
    plot_voltage_response(nswps, used, ind, tvec0, vvec0s, first_dip_pt, base_ind, last_ind, cprwin)
    subplot(2,2,3);
    if ~isempty(xvecF_all)
        plot_cfit('falling', cprwin, cfit_F1, xvecF_all, yvecF_all, params_L_F1, ...
                eqnS_F1, eqnL_F1, tau0_F1, tau1_F1, CL0_F1, CL1_F1, pw_mean);
    end
    subplot(2,2,4);
    if ~isempty(xvecF_mean)
        plot_cfit('falling', cprwin, cfit_F2, xvecF_mean, yvecF_mean, params_L_F2, ...
                eqnS_F2, eqnL_F2, tau0_F2, tau1_F2, CL0_F2, CL1_F2, pw_mean);
    end
    figname = fullfile(outfolder, directories{2}, [filebase, '_cprFalling', suffix, '.png']);
    suptitle(strjoin({'Falling phase of current pulse response for', strrep(filebase, '_', '\_'), title_mod}));
    saveas(h, figname);
    % close(h);

    % Plot input resistance analysis with passive parameter fitting
    scrsz = get(groot, 'ScreenSize');        % screen size
    h = figure(2002);
    set(h, 'Position', [1, scrsz(4)/4, scrsz(3)/2, scrsz(4)*3/4]);
    set(h, 'PaperPositionMode', 'auto')
    set(h, 'Visible', 'off');
    set(h, 'Name', 'Passive parameter fitting');
    clf(h);
    subplot(3,2,1);
    if ~isempty(xvecR_all)
        plot_cfit('rising', cpwin, cfit_R1, xvecR_all, yvecR_all, params_S_R1, ...
                eqnS_R1, eqnL_R1, tau0_R1, tau1_R1);
        text(0.1, 1.1, sprintf('Total number of sweeps: %d', nswps), ...
            'Units', 'normalized');
    end
    subplot(3,2,2);
    if ~isempty(xvecR_mean)
        plot_cfit('rising', cpwin, cfit_R2, xvecR_mean, yvecR_mean, params_S_R2, ...
                eqnS_R2, eqnL_R2, tau0_R2, tau1_R2);
    end
    subplot(3,2,3);
    if ~isempty(xvecF_all)
        plot_cfit('falling', cprwin, cfit_F1, xvecF_all, yvecF_all, params_L_F1, ...
                eqnS_F1, eqnL_F1, tau0_F1, tau1_F1, CL0_F1, CL1_F1, pw_mean);
    end
    subplot(3,2,4);
    if ~isempty(xvecF_mean)
        plot_cfit('falling', cprwin, cfit_F2, xvecF_mean, yvecF_mean, params_L_F2, ...
                eqnS_F2, eqnL_F2, tau0_F2, tau1_F2, CL0_F2, CL1_F2, pw_mean);
    end
    subplot(3,2,5);
    if ~isempty(xvecR_all) && ~isempty(xvecF_all)
        plot_geometry_both(params_S_R1, params_L_F1);
    end
    subplot(3,2,6);
    if ~isempty(xvecR_mean) && ~isempty(xvecF_mean)
        plot_geometry_both(params_S_R2, params_L_F2);
    end
    % figname = fullfile(outfolder, directories{3}, [filebase, '_passive', suffix, '.png']);
    figname_nopng = fullfile(outfolder, directories{3}, [strrep(filebase, '.', 'p'), '_passive', suffix]);
    subplotsqueeze(h, 1.1);
    suptitle(sprintf('Passive parameter fitting for %s %s\n', ...
                strrep(filebase, '_', '\_'), title_mod));
    print(h, figname_nopng, '-dpng', '-r0');
    % saveas(h, figname);
    % close(h);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ft, coeff_start, coeff_lower, coeff_upper] = fit_setup(eqform, C_Est, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range)

% Find the order of the coefficients that will be used by fit()
ft = fittype(eqform);                % fit type
coeff_names = coeffnames(ft);            % Most likely: {'a'; 'b'; 'c'; 'd'}
nco = numel(coeff_names);            % number of coefficients
coeff_start = zeros(1, nco);            % start points of coefficients
coeff_lower = zeros(1, nco);            % lower bounds of coefficients
coeff_upper = zeros(1, nco);            % upper bounds of coefficients

% Set up initial conditions and boundary conditions for fit()
for co = 1:nco
    coeff_lower(co) = 0;
    if coeff_names{co} == 'a' || coeff_names{co} == 'c'
        coeff_start(co) = C_Est;
        coeff_upper(co) = maxScalingFactor * C_Est;
    elseif coeff_names{co} == 'b'
        coeff_start(co) = typtau0;
        coeff_lower(co) = tau0_range(1);
        coeff_upper(co) = tau0_range(2);
    elseif coeff_names{co} == 'd'
        coeff_start(co) = typtau1;
        coeff_lower(co) = tau1_range(1);
        coeff_upper(co) = tau1_range(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmse_swp = measure_error(xvecs, yvecs, eqform, C_Est, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range, pulseWidth);
%% Returns root-mean-squared error for each trace fitted individually

% Setup fitting type, initial conditions and boundaries
[ft, coeff_start, coeff_lower, coeff_upper] = fit_setup(eqform, C_Est, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range);

% For each trace, fit to curve and compute sum-of-squares error
nvecs = numel(xvecs);
cfit_swp = cell(1, nvecs);        % curve fits for each individual sweep
gof_swp = cell(1, nvecs);        % goodness-of-fit statistics for each individual sweep
rmse_swp = zeros(1, nvecs);        % root-mean-squared error (mV) for each individual sweep
parfor k = 1:nvecs
    % Fit data to equation
    if ~isempty(xvecs{k}) && ~isempty(yvecs{k})
        [cfit_swp{k}, gof_swp{k}] = fit(xvecs{k}, yvecs{k}, ft, 'StartPoint', coeff_start, ...
                        'Lower', coeff_lower, 'Upper', coeff_upper); 
        rmse_swp(k) = gof_swp{k}.rmse;
    else
        rmse_swp(k) = Inf;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cfit, eqnS, eqnL, CS0, tau0, CS1, tau1, CL0, CL1] = ...
    fit_data (phase, eqform, xvec, yvec, C_Est, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range, pulseWidth)
%% Returns a string for the equation of the fitted curve

% Setup fitting type, initial conditions and boundaries
[ft, coeff_start, coeff_lower, coeff_upper] = fit_setup(eqform, C_Est, maxScalingFactor, typtau0, typtau1, tau0_range, tau1_range);

% Fit data to equation
[cfit, gof] = fit(xvec, yvec, ft, 'StartPoint', coeff_start, ...
            'Lower', coeff_lower, 'Upper', coeff_upper); 

% Rename the coefficients to be consistent with Johnston & Wu, p. 94
if cfit.b >= cfit.d
    CS0 = cfit.a;
    tau0 = cfit.b;
    CS1 = cfit.c;
    tau1 = cfit.d;
else
    CS0 = cfit.c;
    tau0 = cfit.d;
    CS1 = cfit.a;
    tau1 = cfit.b;
end

% Compute corresponding long pulse response coefficients for the falling phase
if strcmp(phase, 'falling')
    [CL0, CL1] = short2long (CS0, tau0, CS1, tau1, pulseWidth);
elseif strcmp(phase, 'rising')
    CL0 = CS0;
    CL1 = CS1;
else
    CL0 = NaN;
    CL1 = NaN;
end

% Generate strings for the fitted equations
eqnS = strrep(eqform, 'a', num2str(CS0, 2));
eqnS = strrep(eqnS, 'b', num2str(tau0, 2));
eqnS = strrep(eqnS, 'c', num2str(CS1, 2));
eqnS = strrep(eqnS, 'd', num2str(tau1, 2));
eqnL = strrep(eqform, 'a', num2str(CL0, 2));
eqnL = strrep(eqnL, 'b', num2str(tau0, 2));
eqnL = strrep(eqnL, 'c', num2str(CL1, 2));
eqnL = strrep(eqnL, 'd', num2str(tau1, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CL0, CL1] = short2long (CS0, tau0, CS1, tau1, pulseWidth)
%% Estimate the long pulse response coefficients from the 
%   short pulse response coefficients
% See equation 4.5.56 on p. 95 of Johnston % Wu

CL0 = CS0 / (1 - exp(-pulseWidth/tau0));
CL1 = CS1 / (1 - exp(-pulseWidth/tau1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params] = coeff2params (C0, tau0, C1, tau1, cpa_mean, Cm, Ra, Rs)
%% Parameters to be estimated from the fitted coefficients

syms x

% Find Rinput
dvss = -(C0 + C1);                % supposed change in membrane potential at steady state (mV)
Rinput = (dvss*10^-3)/(cpa_mean*10^-12)/10^6;    % input resistance (MOhm)

% Find L, rho
alpha1 = sqrt(tau0/tau1 - 1);
if alpha1 < 0
    error('code is wrong!');
end
if alpha1 == 0                    % tau0 == tau1, no dendrite
    L_init = 0;
    L = 0;
    rho = 0;
else
    L_init = pi/alpha1;            % initial guess for L
    L = double(vpasolve(abs(C1/(2*C0*tau1/tau0-C1)) == cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x)), x, L_init));    
                        % electrotonic length
    rho = -alpha1*cot(alpha1*L)/coth(L);    % dendritic to somatic conductance ratio
    while L < 0 || rho < 0                    % Recompute L if L < 0 or rho < 0
        L_init = rand(1);            % Start with a random number between 0 and 1
        L = double(vpasolve(abs(C1/(2*C0*tau1/tau0-C1)) == cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x)), x, L_init));    
                                % electrotonic length
        rho = -alpha1*cot(alpha1*L)/coth(L);        % dendritic to somatic conductance ratio
    end
end

% Find Rsoma, Rdend
Rmemb = Rinput - Rs;                % input resistance of the cell membrane (MOhm)
Gmemb = 1/Rmemb;                % input conductance (uS)
Gsoma = Gmemb/(1 + rho);            % somatic conductance (uS)
Gdend = Gsoma * rho;                % dendritic conductance (uS)
Rsoma = 1/Gsoma;                % somatic resistance (MOhm)
Rdend = 1/Gdend;                % dendritic resistance (MOhm)

% Find taum, Rm
taum = tau0;                    % membrane time constant (ms)
Rm = 10^3 * taum / Cm;                % specific membrane resistivity (Ohm-cm^2)

% Modeling the soma as a sphere, find radius of soma
%     Rsoma = Rm/(4*pi*r^2)
rad_soma = 10^4 * sqrt(Rm/(4*pi*(10^6*Rsoma)));    % radius of soma (um)

% Modeling the dendrite as a cylinder, find radius of dendrite
%     Rdend = (2/pi) * (Rm*Ra)^(1/2) * d^(-3/2) * coth(L)
if isinf(Rdend)
    diam_dend = 0;
else
    diam_dend = 10^4 * ( (2/pi) * (Rm*Ra)^(1/2) * (10^6*Rdend)^(-1) * coth(L) )^(2/3);    % diameter of dendrite (um)
end
rad_dend = diam_dend / 2;            % radius of dendrite (um)

% Find space constant lambda & length of dendrite
lambda = sqrt((rad_dend*Rm)/(2*Ra));        % space constant (um)
length_dend = L * lambda;            % length of dendrite (um)

% Save parameters into params
params.C0 = C0;
params.tau0 = tau0;
params.C1 = C1;
params.tau1 = tau1;
params.cpa_mean = cpa_mean;
params.dvss = dvss;
params.Rinput = Rinput;
params.alpha1 = alpha1;
params.L_init = L_init;
params.L = L;
params.rho = rho;
params.Rmemb = Rmemb;
params.Rsoma = Rsoma;
params.Rdend = Rdend;
params.Gmemb = Gmemb;
params.Gsoma = Gsoma;
params.Gdend = Gdend;
params.taum = taum;
params.Rm = Rm;
params.rad_soma = rad_soma;
params.diam_dend = diam_dend;
params.rad_dend = rad_dend;
params.length_dend = length_dend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_current_pulse (nswps, used, ind, tvec0, ivec1s, first_dip_pt, base_ind, last_ind, win, cpa_mean)
%% Plots current pulses and display total number of sweeps

for swp = 1:nswps            % Plot each current trace
    if used(swp)
        plot(tvec0(ind), ivec1s(ind, swp), '-'); hold on; 
    else                % if not used, plot as dotted line
        plot(tvec0(ind), ivec1s(ind, swp), '--'); hold on; 
    end
    if ~isempty(first_dip_pt{swp})
        plot(tvec0(base_ind(swp, 1)), ivec1s(base_ind(swp, 1), swp), '>');
        plot(tvec0(base_ind(swp, end)), ivec1s(base_ind(swp, end), swp), '<');
        plot(tvec0(last_ind(swp, 1)), ivec1s(last_ind(swp, 1), swp), '>');
        plot(tvec0(last_ind(swp, end)), ivec1s(last_ind(swp, end), swp), '<');
    end
end
xlim(win);
ylim([cpa_mean*6/5 -cpa_mean*1/5]);
ax = gca;
xlimits = get(ax, 'Xlim');
ylimits = get(ax, 'Ylim');
xpos = xlimits(1) + (1/20) * (xlimits(2) - xlimits(1));
ypos = ylimits(1) + (19/20) * (ylimits(2) - ylimits(1));
text(xpos, ypos, ['Total number of sweeps: ', num2str(nswps)]);
xlabel('Time (ms)')
ylabel('Current (pA)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_voltage_response (nswps, used, ind, tvec0, vvec0s, first_dip_pt, base_ind, last_ind, win)
%% Plots voltage responses to current pulses

for swp = 1:nswps            % Plot each voltage trace
    if used(swp)
        plot(tvec0(ind), vvec0s(ind, swp), '-'); hold on; 
    else                % if not used, plot as dotted line
        plot(tvec0(ind), vvec0s(ind, swp), '--'); hold on; 
    end
    if ~isempty(first_dip_pt{swp})
        plot(tvec0(base_ind(swp, 1)), vvec0s(base_ind(swp, 1), swp), '>');
        plot(tvec0(base_ind(swp, end)), vvec0s(base_ind(swp, end), swp), '<');
        plot(tvec0(last_ind(swp, 1)), vvec0s(last_ind(swp, 1), swp), '>');
        plot(tvec0(last_ind(swp, end)), vvec0s(last_ind(swp, end), swp), '<');
    end
end
xlabel('Time (ms)')
ylabel('Voltage (mV)')
xlim(win);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_cfit (phase, win, cfit, xvec, yvec, params, eqnS, eqnL, tau0, tau1, CL0, CL1, pulseWidth)
%% Plots data along with the fitted curve

% Extract from params
dvss = params.dvss;
Rinput = params.Rinput;

% Plot the data with the fitted curve
plot(cfit, xvec, yvec); hold on;
if strcmp(phase, 'falling')
    % Find time before the short pulse response the corresponding long pulse response would begin
    % CSn = CLn*(1 - exp(-w/taun)) = CLn*exp(-t/taun))
    % => t = -log(1 - exp(-w/taun))*taun
    paddedt0 = -log(1 - exp(-pulseWidth/tau0))*tau0;
    % fprintf('paddedt0 == %g\n', paddedt0);
    paddedt1 = -log(1 - exp(-pulseWidth/tau1))*tau1;
    % fprintf('paddedt1 == %g\n', paddedt1);
    paddedt = max(paddedt0, paddedt1);

    % Plot the "corresponding" long pulse response
    xvec2 = floor(-paddedt):1:400;
    yvec2 = -CL0*exp(-(xvec2+paddedt)/tau0) - CL1*exp(-(xvec2+paddedt)/tau1);
    plot(xvec2, yvec2, 'm');
end

% Create new legend
if strcmp(phase, 'rising')
    legend('data', 'fitted');
elseif strcmp(phase, 'falling')
    legend('data', 'fitted', 'LPR');
end

% Adjust the axes
ylim([(dvss - 0.5) 0.5]);        % Add 0.5 mV above and below
if strcmp(phase, 'rising')
    xlim(win - win(1));
elseif strcmp(phase, 'falling')
    xlim([-paddedt win(2)]);
end

% Plot a line for the asymptote
ax = gca;
xlimits = get(ax, 'Xlim');
ylimits = get(ax, 'Ylim');
if strcmp(phase, 'rising')
    color = 'r';
elseif strcmp(phase, 'falling')
    color = 'm';
end
line(xlimits, [dvss dvss], 'Color', color, 'LineStyle', '--', 'LineWidth', 0.5);

% Showing the equation and Rinput
xpos = xlimits(1) + (1/30) * (xlimits(2) - xlimits(1));
ypos = ylimits(1) + (5/30) * (ylimits(2) - ylimits(1));
text(xpos, ypos, eqnS, 'FontSize', 8, 'Color', 'r');
if strcmp(phase, 'falling')
    text(xpos, ypos - (1/15) * (ylimits(2) - ylimits(1)), eqnL, ...
        'FontSize', 8, 'Color', 'm');
end
text(xpos, ypos - (2/15) * (ylimits(2) - ylimits(1)), ['Rin = ', num2str(Rinput), ' MOhm'], ...
    'FontSize', 8, 'Color', color);

% Add axes labels
xlabel('Shifted time (ms)')
ylabel('Shifted voltage (mV)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_geometry_both (params_S, params_L)
%% Plots geometry of model cell from both the short-pulse parameters and the long-pulse parameters

% Plot model cells
plot_geometry(params_S, 'r');
plot_geometry(params_L, 'm');

% Add axes labels
xlabel('Position (um)')
ylabel('Position (um)')

% Adjust axes
xlim([-100 100]);
ylim([-100 100]);

% Show parameter values
ax = gca;
xlimits = get(ax, 'Xlim');
ylimits = get(ax, 'Ylim');
xpos = xlimits(1) + (1/30) * (xlimits(2) - xlimits(1));
ypos = ylimits(1) + (29/30) * (ylimits(2) - ylimits(1));
text(xpos, ypos, ['L = ', num2str(params_S.L)], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (1/15) * (ylimits(2) - ylimits(1)), ['rho_S = ', num2str(params_S.rho)], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (2/15) * (ylimits(2) - ylimits(1)), ['Rsoma_S = ', num2str(params_S.Rsoma, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (3/15) * (ylimits(2) - ylimits(1)), ['Rdend_S = ', num2str(params_S.Rdend, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (4/15) * (ylimits(2) - ylimits(1)), ['taum_S = ', num2str(params_S.taum), ' ms'], ...
    'FontSize', 8, 'Color', 'r');
text(xpos, ypos - (5/15) * (ylimits(2) - ylimits(1)), ['Rm_S = ', num2str(params_S.Rm), ' Ohm-cm^2'], ...
    'FontSize', 8, 'Color', 'r');

text(xpos, ypos - (9/15) * (ylimits(2) - ylimits(1)), ['L_L = ', num2str(params_L.L)], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (10/15) * (ylimits(2) - ylimits(1)), ['rho_L = ', num2str(params_L.rho)], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (11/15) * (ylimits(2) - ylimits(1)), ['Rsoma_L = ', num2str(params_L.Rsoma, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (12/15) * (ylimits(2) - ylimits(1)), ['Rdend_L = ', num2str(params_L.Rdend, '%.0f'), ' MOhm'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (13/15) * (ylimits(2) - ylimits(1)), ['taum_L = ', num2str(params_L.taum), ' ms'], ...
    'FontSize', 8, 'Color', 'm');
text(xpos, ypos - (14/15) * (ylimits(2) - ylimits(1)), ['Rm_L = ', num2str(params_L.Rm), ' Ohm-cm^2'], ...
    'FontSize', 8, 'Color', 'm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_geometry (params, color)
%% Plots geometry of model cell (ball-and-stick) based on fitted parameters

% Extract parameters from params
rad_soma = params.rad_soma;
rad_dend = params.rad_dend;
length_dend = params.length_dend;

% Plot soma
rectangle('Position', [-rad_soma, -rad_soma, 2*rad_soma, 2*rad_soma], ...
        'Curvature', [1 1], 'EdgeColor', color); hold on;

% Plot dendrite
rectangle('Position', [rad_soma, -rad_dend, length_dend, 2*rad_dend], ...
        'Curvature', [0 0], 'EdgeColor', color);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE

ft_R = fittype('a*(exp(-x/b)-1)+c*(exp(-x/d)-1)');    
ft_F = fittype('-a*exp(-x/b)-c*exp(-x/d)');        

dvss_S = params_S.dvss;
Rinput_S = params_S.Rinput;
if strcmp(phase, 'falling')
    dvss_L = params_L.dvss;
    Rinput_L = params_L.Rinput;
end

    dvss_min = min(dvss_L, dvss_S);        % the more negative of the two steady-state voltage changes
    ylim([(dvss_min - 0.5) 0.5]);        % Add 0.5 mV above and below

line(xlimits, [dvss_S dvss_S], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 0.5);
if strcmp(phase, 'falling')
    line(xlimits, [dvss_L dvss_L], 'Color', 'm', 'LineStyle', '--', 'LineWidth', 0.5);
end

if nargin < 10 || (nargin >= 10 && isempty(fitmode))
    fitmode = 0;
end

cfit = cell(1, nswps);                % double exponential fit

        base_ind(swp, :) = cpstart - round(0.5/sims) - fliplr(mvind);    % base indices

ft = fittype('a*exp(-x/b)+c*exp(-x/d)+e');    % double exponential
% coeff = coeffnames(ft);            % {'a'; 'b'; 'c'; 'd'; 'e'}
    cfit_R1 = fit(xvecR_all, yvecR_all, ft, ...
        'StartPoint', [dv_rec_mean, typtau0, typtau1, dv_rec_mean, typtau0, typtau1, -2 * dv_rec_mean], ...
        'Lower', [0, 0, 0, 0, -20 * dv_rec_mean], ...
        'Upper', [10 * dv_rec_mean, 100, 10 * dv_rec_mean, 100, 0]); 

h = figure('Visible', 'off');

        plot(tvec0(ind), ivec0s(ind, swp), 'b'); hold on; 

ypos = dvss + 1;

    if dv_rec(swp) > 0 && cpa(swp) < 0    % exclude faulty traces, but this is really unnecessary
    end

% L = double(solve(abs(C1/(2*C0*tau1/tau0-C1)) == cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x)), x));    

eqnS = [num2str(CS0, 2), '*(exp(-x/', num2str(tau0, 2), ...
    ')-1)+', num2str(CS1, 2), '*(exp(-x/', num2str(tau1, 2), ')-1)'];
eqnL = [num2str(CL0, 2), '*(exp(-x/', num2str(tau0, 2), ...
    ')-1)+', num2str(CL1, 2), '*(exp(-x/', num2str(tau1, 2), ')-1)'];

typtau = 2;        % typical tau ~ 2 ms
tau_max = 500;        % tau cannot exceed 500 ms

        coeff_start(co) = typtau;

% Set up initial conditions and boundary conditions for fit()
for co = 1:nco
    coeff_lower(co) = 0;
    if coeff_names{co} == 'a' || coeff_names{co} == 'c'
        coeff_start(co) = C_Est;
        coeff_upper(co) = maxScalingFactor * C_Est;
    elseif coeff_names{co} == 'b' || coeff_names{co} == 'd'
        coeff_start(co) = mean(tau_min, tau_max);
        coeff_lower(co) = tau_min;
        coeff_upper(co) = tau_max;
    end
end

    % This won't work
        coeff_start(co) = mean(tau0_range(1), tau0_range(2));
        coeff_start(co) = mean(tau1_range(1), tau1_range(2));

        coeff_start(co) = tau0_range(1);
        coeff_start(co) = tau1_range(2);

%} 
