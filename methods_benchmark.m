clc;
clear;
close all;


GOS = 1; SIRmin = 10;
io = 6; N_Sector = 1;
S = 340; n = 4; Au = 25e-3;

% Calculate cluster size
clusterSize = cellPlanning.getClusterSize(io, SIRmin, n);



% Comparison Between Custom, fzero and fsolve for traffic intensity
AcellCustom = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS,'custom',true);
AcellFsolve = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS,'fsolve',true);
AcellFzero = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS,'fzero',true);
disp(['Traffic Intensity (Custom Function): ', num2str(AcellCustom), ' Erlangs']);
disp(['Traffic Intensity (Fsolve Function): ', num2str(AcellFsolve), ' Erlangs']);
disp(['Traffic Intensity (Fzero  Function): ', num2str(AcellFzero), ' Erlangs']);

% Time the 'custom' method
avgTimeCustom = timeit(@() cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'custom'));

% Time the 'fsolve' method
avgTimeFsolve = timeit(@() cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'fsolve'));

% Time the 'fzero' method
avgTimeFzero = timeit(@() cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'fzero'));

% Time the 'ahmedsaleh' method
avgTimeFsaleh = timeit(@() cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'ahmedsaleh'));


% Convert seconds to milliseconds
avgTimeCustom_ms = avgTimeCustom * 1000;
avgTimeFsolve_ms = avgTimeFsolve * 1000;
avgTimeFzero_ms = avgTimeFzero * 1000;
avgTimeFsaleh_ms = avgTimeFsaleh * 1000;

% Display average execution times in milliseconds
disp(['Average execution time (Custom): ' num2str(avgTimeCustom_ms) ' ms']);
disp(['Average execution time (Fsolve): ' num2str(avgTimeFsolve_ms) ' ms']);
disp(['Average execution time (Fzero): ' num2str(avgTimeFzero_ms) ' ms']);
disp(['Average execution time (ahmedsaleh): ' num2str(avgTimeFsaleh_ms) ' ms']);

