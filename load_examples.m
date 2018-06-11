%% load_examples.m
% Load example data structures for testing
% File History:
% 2018-XX-XX Created by Adam LU
% 2018-06-12 Moved xlsData from Test_xlWrite.m

nRows = 15;
nCols = 15;

% Logical arrays
isMarried = [false; false; true; false; false];

% Character arrays
string = 'I love you';

% Numeric arrays
numericScalar = nCols;
numericCol = ones(nRows, 1);
numericRow = ones(1, nCols);
numericArray = magic(nRows);

% Cell arrays
names = {'Adam', 'Ashley', 'Mark', 'Peter', 'Katie'};

% Structures
info = struct;
info.names = names;
info.isMarried = isMarried;

% Sheet Data
xlsData = {'A Number' 'Boolean Data' 'Empty Cells' 'Strings';...
    1 true [] 'String Text';...
    5 false [] 'Another very descriptive text';...
    -6.26 false 'This should have been an empty cell but I made an error' 'This is text';...
    1e8 true [] 'Last cell with text';...
    1e3 false NaN NaN;...
    1e2 true [] 'test'};

% Normally distributed data
normData = randn(100, 1);

% Normally distributed data with outliers
normDataWOutliers = [normData; -100; 200; 1000];

