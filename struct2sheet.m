function table = struct2sheet (structArray, sheetName)
%% Converts a structure array into a table and write it to a spreadsheet
% Usage: table = struct2sheet (structArray, sheetName)
% Outputs:
%       table       - converted table
% Arguments:
%       structArray - a structure array with homogeneous fields
%                   must be a structure array
%       sheetName   - the file name to write to
% Used by:

% File History:
% 2018-10-03 Created by Adam Lu
% TODO: Input parser
% TODO: Add 'SheetType' as a parameter, using the provided extension
%           as the default, and if none provided, .xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to a table
table = struct2table(structArray);

% Print the table to a file
writetable(table, sheetName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%