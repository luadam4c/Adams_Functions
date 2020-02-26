% metabolismR01_power_analysis.m
%% Performs power analysis for the metabolism R01 grant renewal 20200225
%
% Requires:
%       cd/compute_sampsizepwr.m
%       cd/create_power_tables.m
%       cd/construct_fullpath.m

% File History:
% 2020-01-30 Modified from plethR01_power_analysis.m

%% Hard-coded parameters
% Flags
computeSwd2DGFlag = false;
computeAmpkByFastFlag = false;
computeSwdADrugFlag = false;
computeS783Flag = false;
computeAmpkByIncFlag = false;
computeGababAmpADrugFlag = false;
computeBurstProbADrugFlag = true;
computeOscDurADrugFlag = false;
computeSwdMetforminFlag = false;
computeBloodLacFlag = false;
computeThalamicLacFlag = false;
computeGababAmpLacFlag = false;
computeOscDurLacFlag = false;

% Paths
grantDir = '/media/shareX/2020marchR01';
powerDir = fullfile(grantDir, 'Power_Analysis');
byEffectSuffix = 'power_table_by_effect_size';

% Strings
nSwdsString = 'SWD_count';
ampkRatioString = 'ampk_ratio';
gababAmpString = 'gababAmp';
oscDurString = 'oscDuration';
lacConcString = 'lactateConcentration';

% Statistical powers to test
statPowerDesired = [0.8; 0.85; 0.95; 0.99];

% Parameters for 2DG SWD counts (Aim 1.1)
swd2DGKeyword = '2DG_infusion_chevron_swds_per_3hours';
swd2DGTestType = 't';
swd2DGVarName1 = 'saline';
swd2DGVarName2 = 'twoDG';
swd2DGToTest = 8:2:14;

% Parameters for AMPK activation by fasting (Aim 1.2)
ampkByFastKeyword = 'ampk_activation_by_fasting_estimate';
ampkByFastTestType = 't2';
ampkByFastVarName1 = 'fed';
ampkByFastVarName2 = 'fasted';
ampkByFastToTest = 1.5:0.5:3.5;

% Parameters for A-Drug infusion SWD counts (Aim 1.3)
swdADrugKeyword = 'adrug_infusion_chevron_swds_per_3hours';
swdADrugTestType = 't';
swdADrugVarName1 = 'saline';
swdADrugVarName2 = 'aDrug';
swdADrugToTest = 8:2:14;

% Parameters for s783 (Aim 1.4)
s783FileBase = 's783';
s783TestType = 't2';
s783MeanNull = 1;
s783Estimates = [1.78, 1.85];   % see s793_estimates.xlsx
s783Stdev = [0.49, 0.23];       % see s793_estimates.xlsx

% Parameters for AMPK phosphorylation by slice incubation (Aim 2.1)
ampkByIncFileBase = 'ampk_activation_by_incubation';
ampkByIncTestType = 't2';
ampkByIncMeanNull = 1;
ampkByIncEstimates = 21.3;      % see ampk_activation_by_incubation_estimate.xlsx
ampkByIncStdev = 1.23;          % see ampk_activation_by_incubation_estimate.xlsx

% Parameters for GABAB IPSC amplitude by A-Drug (Aim 2.2)
% gababAmpADrugKeyword = 'adrug_gabab_ipsc_at_15min_chevron';
gababAmpADrugKeyword = 'adrug_gabab_ipsc_at_6min_chevron';
gababAmpADrugTestType = 't2';
gababAmpADrugVarName1 = 'control';
gababAmpADrugVarName2 = 'aDrug';
gababAmpADrugToTest = 110:10:150;

% Parameters for burst probability by A-Drug (Aim 2.3)
aDrugBurstProbFileBase = 'adrug_burst_probability';
aDrugBurstProbTestTypes = {'p', 't'};

% Parameters for oscillation duration by A-Drug (Aim 2.4)
oscDurADrugKeyword = 'adrug_oscDurationSec_chevron';
oscDurADrugTestType = 't';
oscDurADrugVarName1 = 'baseline';
oscDurADrugVarName2 = 'aDrug';
oscDurADrugToTest = 1:0.5:4;

% Parameters for metformin SWD counts (Aim 3.1a)
swdMetKeyword = 'metformin_injection_chevron_swds_per_hour_timepoint3';
swdMetTestType = 't';
swdMetVarName1 = 'saline';
swdMetVarName2 = 'metformin';
swdMetToTest = 8:2:14;

% Parameters for blood lactate measurements (Aim 3.1b)
bloodLacFileBase = 'blood_lactate';
bloodLacTestType = 't';
bloodLacMeanNull = 0;
bloodLacEstimates = 7.311;      % in mM, see blood_lactate_estimates.xlsx
bloodLacStdev = 2.783;          % in mM, see blood_lactate_estimates.xlsx

% Parameters for thalamic lactate measurements (Aim 3.1c)
thalamicLacKeyword = 'thalamus_lactate_estimates';
thalamicLacTestType = 't';
thalamicLacVarName1 = 'baseline';
thalamicLacVarName2 = 'metformin';
thalamicLacToTest = 100:50:300;

% Parameters for GABAB IPSC amp by lactate agonist (Aim 3.2)
gababAmpLacKeyword = 'lactate_agonist_gabab_ipsc_chevron';
gababAmpLacTestType = 't';
gababAmpLacVarName1 = 'baseline';
gababAmpLacVarName2 = 'agonist';
gababAmpLacToTest = 8:11;

% Parameters for oscillation duration by lactate agonist (Aim 3.3)
oscDurLacKeyword = 'lactate_agonist_oscDurationSec_chevron';
oscDurLacTestType = 't';
oscDurLacVarName1 = 'baseline';
oscDurLacVarName2 = 'agonist';
oscDurLacToTest = 1:0.5:4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Power analysis for 2-DG infusion SWD counts
if computeSwd2DGFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(swd2DGKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, swd2DGToTest, nSwdsString, ...
                        swd2DGVarName1, swd2DGVarName2, ...
                        'TestType', swd2DGTestType);
    useLog2Ratio = true;
    create_power_tables(swd2DGKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, log2(swd2DGToTest), nSwdsString, ...
                        swd2DGVarName1, swd2DGVarName2, ...
                        'TestType', swd2DGTestType);
end

%% Power analysis for AMPK phosphorylation by fasting
if computeAmpkByFastFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(ampkByFastKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, ampkByFastToTest, ampkRatioString, ...
                        ampkByFastVarName1, ampkByFastVarName2, ...
                        'TestType', ampkByFastTestType);
    useLog2Ratio = true;
    create_power_tables(ampkByFastKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, log2(ampkByFastToTest), ...
                        ampkRatioString, ...
                        ampkByFastVarName1, ampkByFastVarName2, ...
                        'TestType', ampkByFastTestType);
end

%% Power analysis for A-Drug infusion SWD counts
if computeSwdADrugFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(swdADrugKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, swdADrugToTest, nSwdsString, ...
                        swdADrugVarName1, swdADrugVarName2, ...
                        'TestType', swdADrugTestType);
    useLog2Ratio = true;
    create_power_tables(swdADrugKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, log2(swdADrugToTest), nSwdsString, ...
                        swdADrugVarName1, swdADrugVarName2, ...
                        'TestType', swdADrugTestType);
end

%% Power analysis for S783 phosphorylation
if computeS783Flag
    % Compute the power table
    s783Table = compute_sampsizepwr('MeanNull', s783MeanNull, ...
                'Stdev', s783Stdev, 'MeanAlt', s783Estimates, ...
                'TestType', s783TestType, 'OutputType', 'scalars');

    % Write the table
    writetable(s783Table, fullfile(powerDir, [s783FileBase, '_', ...
                                    byEffectSuffix, '.csv']));
end

%% Power analysis for AMPK phosphorylation by slice incubation
if computeAmpkByIncFlag
    % Compute the power table
    ampkByIncTable = compute_sampsizepwr('MeanNull', ampkByIncMeanNull, ...
                'Stdev', ampkByIncStdev, 'MeanAlt', ampkByIncEstimates, ...
                'TestType', ampkByIncTestType, 'OutputType', 'scalars');

    % Write the table
    writetable(ampkByIncTable, fullfile(powerDir, [ampkByIncFileBase, '_', ...
                                        byEffectSuffix, '.csv']));
end

%% Power analysis for GABAB IPSC amplitude by A-Drug
if computeGababAmpADrugFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(gababAmpADrugKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, gababAmpADrugToTest, gababAmpString, ...
                        gababAmpADrugVarName1, gababAmpADrugVarName2, ...
                        'TestType', gababAmpADrugTestType);
end

%% Power analysis for burst probability by A-Drug
if computeBurstProbADrugFlag
    [aDrugBurstProbNullValue, aDrugBurstProbAltValue, aDrugBurstProbStdev] = ...
        metabolismR01_extrapolate_burst_probability;

    % Compute the power table
    aDrugBurstProbTable = compute_sampsizepwr('OutputType', 'scalars', ...
                'MeanNull', aDrugBurstProbNullValue, ...
                'MeanAlt', aDrugBurstProbAltValue, ...
                'Stdev', aDrugBurstProbStdev, ...
                'TestType', aDrugBurstProbTestTypes);

    % Write the table
    writetable(aDrugBurstProbTable, ...
                fullfile(powerDir, [aDrugBurstProbFileBase, '_', ...
                                        byEffectSuffix, '.csv']));
end

%% Power analysis for oscillation duration by A-Drug
if computeOscDurADrugFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(oscDurADrugKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, oscDurADrugToTest, oscDurString, ...
                        oscDurADrugVarName1, oscDurADrugVarName2, ...
                        'TestType', oscDurADrugTestType);
end

%% Power analysis for metformin infusion SWD counts
if computeSwdMetforminFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(swdMetKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, swdMetToTest, nSwdsString, ...
                        swdMetVarName1, swdMetVarName2, ...
                        'TestType', swdMetTestType);
    useLog2Ratio = true;
    create_power_tables(swdMetKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, log2(swdMetToTest), nSwdsString, ...
                        swdMetVarName1, swdMetVarName2, ...
                        'TestType', swdMetTestType);
end

%% Power analysis for blood lactate measurements
if computeBloodLacFlag
    % Compute the power table
    bloodLacTable = compute_sampsizepwr('MeanNull', bloodLacMeanNull, ...
                'Stdev', bloodLacStdev, 'MeanAlt', bloodLacEstimates, ...
                'TestType', bloodLacTestType, 'OutputType', 'scalars');

    % Write the table
    writetable(bloodLacTable, fullfile(powerDir, [bloodLacFileBase, '_', ...
                                        byEffectSuffix, '.csv']));
end

%% Power analysis for thalamic lactate measurements
if computeThalamicLacFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(thalamicLacKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, thalamicLacToTest, lacConcString, ...
                        thalamicLacVarName1, thalamicLacVarName2, ...
                        'TestType', thalamicLacTestType);
end

%% Power analysis for GABAB IPSC amplitude by lactate agonist
if computeGababAmpLacFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(gababAmpLacKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, gababAmpLacToTest, gababAmpString, ...
                        gababAmpLacVarName1, gababAmpLacVarName2, ...
                        'TestType', gababAmpLacTestType);
end

%% Power analysis for oscillation duration by lactate agonist
if computeOscDurLacFlag
    % Compute power tables
    useLog2Ratio = false;
    create_power_tables(oscDurLacKeyword, powerDir, statPowerDesired, ...
                        useLog2Ratio, oscDurLacToTest, oscDurString, ...
                        oscDurLacVarName1, oscDurLacVarName2, ...
                        'TestType', oscDurLacTestType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%