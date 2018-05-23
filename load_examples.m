%% load_examples.m
% Load example data structures for testing

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

