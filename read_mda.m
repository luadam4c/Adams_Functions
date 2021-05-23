function A = read_mda(pathName)
%% Reads the contents of a .mda file. MDA stands for multi-dimensional array
% Usage: A = read_mda(pathName)
%
% See http://magland.github.io//articles/mda-format/
%
% Syntax: A = read_mda(pathName)
%
% Inputs:
%    pathName - path to the .mda file
%
% Outputs:
%    A - the multi-dimensional array
%
% Other m-files required: none
%
% See also: writemda

% Author: Jeremy Magland
% Jan 2015; Last revision: 15-Feb-2106
% Revisions by Adam Lu:
% 2021-05-19 Now loads a column at a time if nDim == 2

% If the file is a csv file, just use textread()
if strcmp(pathName(end-4:end), '.csv') == 1
    A = textread(pathName, '', 'Delimiter', ',');
    return;
end

F = fopen(pathName, 'rb');

try
    code = fread(F, 1, 'int32');
catch
    error('Problem reading file: %s', pathName);
end

if code > 0 
    nDim = code;
    code = -1;
else
    fread(F, 1, 'int32');
    nDim = fread(F, 1, 'int32');    
end

dimTypeStr = 'int32';
if nDim < 0
    nDim = -nDim;
    dimTypeStr = 'int64';
end

S = zeros(1, nDim);
for j = 1:nDim
    S(j) = fread(F, 1, dimTypeStr);
end
N = prod(S);

if nDim == 1
    A = zeros(1, S);
else
    A = zeros(S);
end

switch code
    case -1
        precision = 'float';
    case -2
        precision = 'uchar';
    case -3
        precision = 'float';
    case -4
        precision = 'int16';
    case -5
        precision = 'int32';
    case -6
        precision = 'uint16';
    case -7
        precision = 'double';
    case -8
        precision = 'uint32';
    otherwise
        error('Unsupported data type code: %d', code);
end

if code == -1
    M = zeros(1, N*2);
    M(:) = fread(F, N*2, 'float');
    A(:) = M(1:2:prod(S)*2)+i*M(2:2:prod(S)*2);
else
    if nDim == 2
        for iRow = 1:S(1)
            A(iRow, :) = fread(F, S(2), precision);
        end
    else
        A(:) = fread(F, N, precision);
    end
end

fclose(F);
