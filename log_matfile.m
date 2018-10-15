function log_matfile(fullmatfilepath)
%% Print variables in a MATfile to a comma-separated-value file
% Usage: log_matfile(fullmatfilepath)
%
% TODO: check that matfile contains logheader & logvariables
% TODO: annotate
%
% Requires:
%       /home/Matlab/Adams_Functions/print_next_in_csv.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/remediate_dclampdatalog_take4_20171227.m
%
% 2017-12-27 Created

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: make argument either a full path or a relative path from pwd
% matfilename = 'take4/dclampdatalog_take4.mat';
% fullmatfilepath = fullfile(pwd, matfilename);

% Print status
fprintf('Printing %s as a csv file ... \n', fullmatfilepath);

% Extract directory and filebase of MATfile
[folder, filebase, ~] = fileparts(fullmatfilepath);

% Load variables from MATfile
load(fullmatfilepath);

% Open a comma-separated-value file with same filebase for writing
csvfilepath = fullfile(folder, [filebase, '.csv']);
fid = fopen(csvfilepath, 'w');
fprintf('Creating csv file %s ... \n', csvfilepath);

% Determine the number of variables to print
nHeaders = numel(logheader);
nVars = numel(logvariables);
% TODO: Check if nHeaders == nVars

% TODO: Check if all the variables in logvariables exist in MATfile

% Print header
for iVar = 1:nVars
    fprintf(fid, '%s', logheader{iVar});
    print_next_in_csv(fid, iVar, nVars);
end

% Print names of variables
for iVar = 1:nVars
    fprintf(fid, '%s', logvariables{iVar});
    print_next_in_csv(fid, iVar, nVars);
end

% Place all variables in a cell array
allVar = cell(nVars, 1);
for iVar = 1:nVars
    eval(['allVar{iVar} = ', logvariables{iVar}, ';']);
end

% TODO: Check if all variables have the same number of entries
% Determine the number of entries for each variable
nEntries = numel(allVar{1});

% Print one line for each entry
for iEnt = 1:nEntries
    for iVar = 1:nVars
        % Print the entry based on type
        if iscell(allVar{iVar}) && ischar(allVar{iVar}{iEnt})
            fprintf(fid, '%s', allVar{iVar}{iEnt});
        elseif iscell(allVar{iVar})
            fprintf(fid, '%g', allVar{iVar}{iEnt});
        elseif isnumeric(allVar{iVar})
            fprintf(fid, '%g', allVar{iVar}(iEnt));
        else
            error('Cannot print entry %d, variable %d!!', ...
                    iEnt, iVar);
        end
        print_next_in_csv(fid, iVar, nVars);
    end    
end

% Close the comma-separated-value file
fclose(fid);
fprintf('Closed %s!\n', csvfilepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

eval(['tempVar1 = ', logvariables{1}, ';']);
eval(['fprintf(fid, ''%s, '', ', logvariables{iVar}, '{iSwp});']);
eval(['fprintf(fid, ''%d, '', ', logvariables{iVar}, '(iSwp));']);
eval(['fprintf(fid, ''%g, '', ', logvariables{iVar}, '(iSwp));']);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%