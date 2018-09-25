%% load_examples.m
% Load example data structures for testing
% File History:
% 2018-XX-XX Created by Adam LU
% 2018-06-12 Moved xlsData from Test_xlWrite.m

nRows = 3;
nCols = 3;

% Logical arrays
myLogicalScalar = true;
myLogicalRow = [true, false, true, true, false];
myLogicalCol = [false; false; true; false; false];
myLogical2D = [true, false; false, true];

% Character arrays
myString = 'I love you';

% Numeric arrays
myNumericScalar = nCols;
myNumericCol = ones(nRows, 1);
myNumericRow = ones(1, nCols);
myNumeric2D = magic(nRows);
myNumeric3D = floor(rand(3, 4, 2) * 100);

% Cell arrays of character arrays
names = {'Adam', 'Ashley', 'Mark', 'Peter', 'Katie'};
myCellStrScalar = {'I love you'};
myCellStrRow = {'I love you', 'You love me', 'Blab hooray!', 'Why?'};
myCellStrCol = {'I love you'; 'You love me'; 'Blab hooray!'; 'Why?'};
myCellStr2D = {'I love you', 'You love me'; 'Blab hooray!', 'Why?'};

% Cell arrays of cell arrays
myCellCell = {myCellStrRow, myCellStrRow, myCellStrRow};

% Structures
blab = struct;
blab.students = names;
blab.isMarried = myLogicalCol;
myStruct = struct;
myStruct.myLogicalScalar = myLogicalScalar;
myStruct.myLogicalCol = myLogicalCol;
myStruct.myLogicalRow = myLogicalRow;
myStruct.myLogical2D = myLogical2D;
myStruct.myNumericScalar = myNumericScalar;
myStruct.myNumericCol = myNumericCol;
myStruct.myNumericRow = myNumericRow;
myStruct.myNumeric2D = myNumeric2D;
myStruct.names = names;
myStruct.myCellStrScalar = myCellStrScalar;
myStruct.myCellStrRow = myCellStrRow;
myStruct.myCellStrCol = myCellStrCol;
myStruct.myCellStr2D = myCellStr2D;
myStruct.blab = blab;

% Structure arrays
myStructArray = [myStruct, myStruct, myStruct];

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

% A random time series
myRandomSignal = rand(100, 1);
