fileName = 'C:\Users\Pinn Analysis\Desktop\Shinnosuke\data\cobalt multi 07192018.adicht';

% Read in the data from the .adicht file
[data, file] = read_adicht(fileName);

% Analysis on data
signalAnalyzer(data);