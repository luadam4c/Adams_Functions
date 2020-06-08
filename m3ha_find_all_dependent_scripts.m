%% Hard-coded parameters 
scriptDir = fullfile('/home', 'Matlab', 'Adams_Functions');
archiveDir = fullfile('/media', 'adamX', 'm3ha', 'manuscript', 'published_code', 'matlab_code');

%% Find all m3ha master scripts
% Find all paths
[~, m3haScriptPaths] = all_files('Directory', scriptDir, 'Prefix', 'm3ha', ...
                            'Extension', 'm');

% Extract script names
m3haScriptNames = extract_fileparts(m3haScriptPaths, 'base');

% Remove those scripts related to xolotl
m3haScriptNames = m3haScriptNames(~contains(m3haScriptNames, 'xolotl'));

%% Find all dependent scripts
% Extract dependent scripts
m3haDependentScriptPathsAll = ...
    array_fun(@(x) all_dependent_functions(x, 'OriginalOutput', true), ...
            m3haScriptNames, 'UniformOutput', false);

% Find the union
m3haDependentScriptPaths = union_over_cells(m3haDependentScriptPathsAll);

%% Copy paths to code archive directory
copy_into(m3haDependentScriptPaths, archiveDir);