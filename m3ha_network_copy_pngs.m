%% Hard-coded parameters
networdDir = '/media/adamX/m3ha/network_model/';
% iterDirName = '20200312T0130_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_TCepas_varied';
% iterDirName = '20200403T2312_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_TCepas_varied';
% iterDirName = '20200407_using_bestparams_20200203_manual_singleneuronfitting0-102';
% iterDirName = '20200408_using_bestparams_20200203_manual_singleneuronfitting0-102';
% iterDirName = '20200424_using_bestparams_20200203_manual_singleneuronfitting0-102';
% iterDirName = '20200425_using_bestparams_20200203_manual_singleneuronfitting0-102';
% iterDirName = '20200430_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_spikes';
% iterDirName = '20200501_using_bestparams_20200203_manual_singleneuronfitting0-102_2cell_spikes';
% iterDirName = '20200503_using_bestparams_20200203_manual_singleneuronfitting0-102_200cell_spikes';
iterDirName = '20200504_using_bestparams_20200203_manual_singleneuronfitting0-102_hetero_spikes';
iterDir = fullfile(networdDir, iterDirName);
% pngsDir = fullfile(iterDir, 'pngfiles');
pngsDir = fullfile(iterDir, 'rasters');
seedNumbers = transpose(0:79);

pCondStrs = create_labels_from_numbers(1:4, 'Prefix', 'pCond_');
condStrs = strcat(pCondStrs, '_', 'gIncr_16.6667');
condDirs = fullfile(pngsDir, condStrs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create output directories
check_dir(condDirs);

%% Find the rasters
cellfun(@(condDir) copy_pngs_for_cond(iterDir, condDir, seedNumbers), ...
        condDirs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_pngs_for_cond (iterDir, condDir, seedNumbers)

arrayfun(@(seedNum) copy_pngs_for_seedNum(iterDir, condDir, seedNum), ...
            seedNumbers);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_pngs_for_seedNum(iterDir, condDir, seedNum)

seedDir = fullfile(iterDir, ['seedNumber_', num2str(seedNum)]);

condStr = extract_fileparts(condDir, 'name');

keyword = [condStr, '_raster_plot'];
% keyword = condStr;

[~, copyFrom] = all_files('Directory', seedDir, 'Keyword', keyword, ...
                        'Extension', 'png', 'Recursive', true);

if isempty(copyFrom)
    return
end

candLabelRegExp = 'candidateIDs_[0-9,-]*';

pngBase = extract_fileparts(copyFrom, 'base');
cellStr = extract_substrings(copyFrom, 'RegExp', candLabelRegExp);

TCepas = -75 + mod(seedNum, 16);

check_dir(fullfile(condDir, cellStr));

copyTo = fullfile(condDir, cellStr, ...
                strcat('TCepas', num2str(TCepas), ...
                        '_', 'seedNumber', num2str(seedNum), '_', cellStr, ...
                        '_', pngBase, '.png'));

[copyFrom, copyTo] = argfun(@force_column_cell, copyFrom, copyTo);

cellfun(@(a, b) copyfile(a, b), copyFrom, copyTo);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
