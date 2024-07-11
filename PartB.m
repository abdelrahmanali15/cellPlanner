clc;
clear;
close all;
saveFig = false; %true for saving false to not save

%% PART B Figure 1
%------------------------ PART B Start ------------------------%
% 1) Cluster Size vs SIRmin
% Given Parameters
n = 4;

% Parameters
SIRmin_range = 1:0.01:30; % sweep over SIR in dB
colors = linspecer(3); % Choose colors from colormap
clusterSizes = zeros(length(SIRmin_range), 3);

% Calculate cluster sizes
for i = 1:length(SIRmin_range)
    clusterSizes(i, 1) = cellPlanning.getClusterSize(6, SIRmin_range(i), n);
    clusterSizes(i, 2) = cellPlanning.getClusterSize(2, SIRmin_range(i), n);
    clusterSizes(i, 3) = cellPlanning.getClusterSize(1, SIRmin_range(i), n);
end

% Plot
figure('Position', [100, 100, 700, 500]);
for i = 1:3
    plot(SIRmin_range, clusterSizes(:, i), 'Color', colors(i, :), 'LineWidth', 3);
    hold on;
end

xlabel('SIRmin (dB)');
ylabel('Cluster Size');
legend('Omni-directional', '120° Sectorization', '60° Sectorization');
title('Cluster Size vs SIRmin');
% Insert ticks at 14 and 19 on the x-axis
xticks([0 5 10 14 15 19 20 25 30]);
grid on;
% Plot vertical line at sir = 14
line([14, 14], [0, 40], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');



% Plot vertical line at sir = 14
line([19, 19], [0, 40], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off');
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf,  fullfile('plots','Cluster_Size_vs_SIRmin.png')); % Saves the figure as a PNG file
end



%% PART B Figure 2 and 3
% 2)i)  Plot the number of cells versus GOS -- SIR=19
%   ii) Plot the traffic intensity per cell versus GOS
S = 340; n = 4; Au = 25e-3;

user_density = 1400;
city_area = 100;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
SIRmin = 19;
colors = linspecer(3);
labels = {'Omni-directional', '120° Sectorization', '60° Sectorization'};

GOS = 1:0.01:30; % in percentage

% Specify figure size
figure('Position', [10, 10, 700, 800]);

for i = 1:numel(io_vector)
    % Cluster size calculation
    clusterSize = cellPlanning.getClusterSize(io_vector(i), SIRmin, n);
    
    A_cell = zeros(size(GOS)); % Preallocate the result array       
    % Calculate A_actual for each gosRequired
    for j = 1:numel(GOS)
        A_cell(j) = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(i), GOS(j),'custom');
    end
    
    % Calculate number of cells
    [~, Number_of_cells] = cellPlanning.calculateParameters(A_cell, Au, user_density, city_area);
    
    % Plot number of cells versus GOS
    subplot(2, 1, 1);
    plot(GOS, Number_of_cells,  'Color', colors(i, :), 'LineWidth', 3);
    hold on;
    
    % Plot traffic intensity per cell versus GOS
    subplot(2, 1, 2);
    plot(GOS, A_cell, 'Color', colors(i, :), 'LineWidth', 3);
    hold on;
end

% Configure plots
subplot(2, 1, 1);
xlabel('GOS (%)');
ylabel('Number of cells');
legend(labels);
title('Number of cells vs GOS: SNR = 19 dB');
grid on;

subplot(2, 1, 2);
xlabel('GOS (%)');
ylabel('Traffic Intensity per Cell ');
legend(labels);
title('Traffic Intensity per Cell vs GOS: SNR = 19 dB');
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf, fullfile('plots', 'Plot_Number_of_cells_and_Traffic_Intensity_per_Cell_vs_GOS_19SNR.png')); % Saves the figure as a PNG file
end


%% Part B figure 4 and 5
% 3)i)  Plot the number of cells versus GOS SIR=14
%   ii) Plot the traffic intensity per cell versus GOS

S = 340; n = 4; Au = 25e-3;

user_density = 1400;
city_area = 100;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
SIRmin = 14;
colors = linspecer(3);
labels = {'Omni-directional', '120° Sectorization', '60° Sectorization'};

GOS = 1:0.01:30; % in percentage

% Specify figure size
figure('Position', [10, 10, 700, 800]);

for i = 1:numel(io_vector)
    % Cluster size calculation
    clusterSize = cellPlanning.getClusterSize(io_vector(i), SIRmin, n);

    A_cell = zeros(size(GOS)); % Preallocate the result array       
    % Calculate A_actual for each gosRequired
    for j = 1:numel(GOS)
        A_cell(j) = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(i), GOS(j),'custom');
    end
    
    % Calculate number of cells
    [~, Number_of_cells] = cellPlanning.calculateParameters(A_cell, Au, user_density, city_area);
    
    % Plot number of cells versus GOS
    subplot(2, 1, 1);
    plot(GOS, Number_of_cells,  'Color', colors(i, :), 'LineWidth', 3);
    hold on;
    
    % Plot traffic intensity per cell versus GOS
    subplot(2, 1, 2);
    plot(GOS, A_cell, 'Color', colors(i, :), 'LineWidth', 3);
    hold on;
end

% Configure plots
subplot(2, 1, 1);
xlabel('GOS (%)');
ylabel('Number of cells');
legend(labels);
title('Number of cells vs GOS: SNR = 14 dB');
grid on;

subplot(2, 1, 2);
xlabel('GOS (%)');
ylabel('Traffic Intensity per Cell');
legend(labels);
title('Traffic Intensity per Cell vs GOS: SNR = 14 dB');
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf,  fullfile('plots','Plot_Number_of_cells_and_Traffic_Intensity_per_Cell_vs_GOS_14SNR.png')); % Saves the figure as a PNG file
end



%% Extra: Plot to compare between SNR = 14 and 19
S = 340; n = 4; Au = 25e-3;

% Parameters
user_density = 1400;
city_area = 100;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
colors = linspecer(6);
SIRmin_values = [19, 14]; % SIRmin values to loop over for 2 values to compare
GOS = 1:0.5:30; % in percentage

% Specify figure size
figure('Position', [10, 10, 700, 800]);


% Plotting in 2 subplots all curves for each SIR to compare
for j = 1:numel(SIRmin_values)
    % Set SIRmin
    SIRmin = SIRmin_values(j);
    
    for i = 1:numel(io_vector)
        % Cluster size calculation
        clusterSize = cellPlanning.getClusterSize(io_vector(i), SIRmin, n);
        
        % Calculate traffic intensity in a Cell
        A_cell = zeros(size(GOS)); % Preallocate the result array       
        
        % Calculate A_actual for each gosRequired
        for k = 1:numel(GOS)
            A_cell(k) = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(i), GOS(k),'custom');
        end
        
        % Calculate number of cells
        [~, Number_of_cells] = cellPlanning.calculateParameters(A_cell, Au, user_density, city_area);
        
        % Plot number of cells versus GOS
        % First set of curves use colors different from 2nd set
        subplot(2, 1, 1);
        if (j==1)
            plot(GOS, Number_of_cells,  'Color', colors(i, :) ,'Marker', 'o', 'MarkerSize', 8);
        else
            plot(GOS, Number_of_cells,  'Color', colors(i+3, :), 'LineWidth', 3);
        end
        hold on;
        
        % Plot traffic intensity per cell versus GOS
        % First set of curves use colors different from 2nd set
        subplot(2, 1, 2);
        if (j==1)
            plot(GOS, A_cell, 'Color', colors(i, :), 'Marker', 'o', 'MarkerSize', 8);
        else
            plot(GOS, A_cell, 'Color', colors(i+3, :), 'LineWidth',3);
        end
        
        hold on;
    end
end

% Configure plots
subplot(2, 1, 1);
xlabel('GOS (%)');
ylabel('Number of cells');
legend('Omni-directional (SIRmin = 19)', '120° Sectorization (SIRmin = 19)', '60° Sectorization (SIRmin = 19)', ...
    'Omni-directional (SIRmin = 14)', '120° Sectorization (SIRmin = 14)', '60° Sectorization (SIRmin = 14)');
title('Number of cells vs GOS');
grid on;

subplot(2, 1, 2);
xlabel('GOS (%)');
ylabel('Traffic Intensity per Cell');
legend('Omni-directional (SIRmin = 19)', '120° Sectorization (SIRmin = 19)', '60° Sectorization (SIRmin = 19)', ...
    'Omni-directional (SIRmin = 14)', '120° Sectorization (SIRmin = 14)', '60° Sectorization (SIRmin = 14)');
title('Traffic Intensity per Cell vs GOS');
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf, fullfile('plots', 'Plot_Number_of_cells_and_Traffic_Intensity_per_Cell_vs_GOS_sir_14_19.png')); % Saves the figure as a PNG file
end


%% Part B figure 6 and 7
% 4) (i)  Plot the number of cells versus user density (100 to 2000 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2
% 4) (ii) Plot the cell radius versus user density (100 to 2000 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2
S = 340; n = 4; Au = 25e-3;

% Parameters
SIRmin = 14;
GOS = 2;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
userDensity_range = 100:50:2000; % Users per km^2
colors = linspecer(3);
labels = {'Omni-directional', '120° Sectorization', '60° Sectorization'};

% Initialize arrays to store results
num_cells = zeros(length(userDensity_range), numel(io_vector));
cell_radius = zeros(length(userDensity_range), numel(io_vector));
figure('Position', [10, 10, 700, 800]);
for j = 1:numel(io_vector)
    % Calculate number of cells and cell radius for each user density
    for i = 1:length(userDensity_range)
        % Calculate cluster size
        clusterSize = cellPlanning.getClusterSize(io_vector(j), SIRmin, n);
        
        % Calculate traffic intensity in a Cell
        A_cell = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(j), GOS,'custom');
        
        % Calculate number of cells
        [cell_radius(i, j), num_cells(i, j)] = cellPlanning.calculateParameters(A_cell, Au, userDensity_range(i), city_area);
        
    end
    
    % Plot number of cells versus user density
    subplot(2, 1, 1);
    plot(userDensity_range, num_cells(:, j), 'Color', colors(j, :), 'LineWidth', 3);
    hold on;
    
    % Plot cell radius versus user density
    subplot(2, 1, 2);
    plot(userDensity_range, cell_radius(:, j), 'Color', colors(j, :), 'LineWidth', 3);
    hold on;
end

% Plot configurations
subplot(2, 1, 1);
xlabel('User Density (users/km^2)');
ylabel('Number of Cells');
title('Number of Cells vs User Density: SNR = 14 dB');
legend(labels);
grid on;

subplot(2, 1, 2);
xlabel('User Density (users/km^2)');
ylabel('Cell Radius (km)');
title('Cell Radius vs User Density: SNR = 14 dB');
legend(labels);
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf,  fullfile('plots','Number_of_cells_and_Cell_radius_vs_User_Density_SIR14.png'));
end
%% Part B for snr = 19
% 4) (i)  Plot the number of cells versus user density (100 to 2000 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2
% 4) (ii) Plot the cell radius versus user density (100 to 2000 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2
S = 340; n = 4; Au = 25e-3;

% Parameters
SIRmin = 19;
GOS = 2;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
userDensity_range = 100:50:2000; % Users per km^2
colors = linspecer(3);
labels = {'Omni-directional', '120° Sectorization', '60° Sectorization'};

% Initialize arrays to store results
num_cells = zeros(length(userDensity_range), numel(io_vector));
cell_radius = zeros(length(userDensity_range), numel(io_vector));
figure('Position', [10, 10, 700, 800]);
for j = 1:numel(io_vector)
    % Calculate number of cells and cell radius for each user density
    for i = 1:length(userDensity_range)
        % Calculate cluster size
        clusterSize = cellPlanning.getClusterSize(io_vector(j), SIRmin, n);
        
        % Calculate traffic intensity in a Cell
        A_cell = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(j), GOS,'custom');
        
        % Calculate number of cells
        [cell_radius(i, j), num_cells(i, j)] = cellPlanning.calculateParameters(A_cell, Au, userDensity_range(i), city_area);
        
    end
    
    % Plot number of cells versus user density
    subplot(2, 1, 1);
    plot(userDensity_range, num_cells(:, j), 'Color', colors(j, :), 'LineWidth', 3);
    hold on;
    
    % Plot cell radius versus user density
    subplot(2, 1, 2);
    plot(userDensity_range, cell_radius(:, j), 'Color', colors(j, :), 'LineWidth', 3);
    hold on;
end

% Plot configurations
subplot(2, 1, 1);
xlabel('User Density (users/km^2)');
ylabel('Number of Cells');
title('Number of Cells vs User Density: SNR = 19 dB');
legend(labels);
grid on;

subplot(2, 1, 2);
xlabel('User Density (users/km^2)');
ylabel('Cell Radius (km)');
title('Cell Radius vs User Density: SNR = 19 dB');
legend(labels);
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf,  fullfile('plots','Number_of_cells_and_Cell_radius_vs_User_Density_SIR19.png'));
end

%% Extra: Plot to compare between SNR = 14 and 19
S = 340; n = 4; Au = 25e-3;

SIRmin_values = [14, 19]; % SIRmin values to loop over
GOS = 2;
io_vector = [6, 2, 1];
N_Sectors_vector = [1, 3, 6];
userDensity_range = 100:50:2000; % Users per km^2
colors = linspecer(6);
labels = cell(numel(io_vector)*numel(SIRmin_values), 1); % Initialize cell array for legend labels

% Initialize arrays to store results
num_cells = zeros(length(userDensity_range), numel(io_vector), numel(SIRmin_values));
cell_radius = zeros(length(userDensity_range), numel(io_vector), numel(SIRmin_values));

% Specify figure size
figure('Position', [10, 10, 700, 800]);
for k = 1:numel(SIRmin_values)
    SIRmin = SIRmin_values(k);
    for j = 1:numel(io_vector)
        % Calculate number of cells and cell radius for each user density
        for i = 1:length(userDensity_range)
            % Calculate cluster size
            clusterSize = cellPlanning.getClusterSize(io_vector(j), SIRmin, n);
            
            % Calculate traffic intensity in a Cell
            A_cell = cellPlanning.calculateACell(S, clusterSize, N_Sectors_vector(j), GOS,'custom');
            
            % Calculate number of cells
            [cell_radius(i, j, k), num_cells(i, j, k)] = cellPlanning.calculateParameters(A_cell, Au, userDensity_range(i), city_area);
            
        end
        
        % Plot number of cells versus user density
        subplot(2, 1, 1);
        if (k==1)
            plot(userDensity_range, squeeze(num_cells(:, j, k)), 'Color', colors(j, :), 'Marker', 'o', 'MarkerSize', 8);
        else
            
            plot(userDensity_range, squeeze(num_cells(:, j, k)), 'Color', colors(j+3, :), 'LineWidth',3);
        end
        
        hold on;
        
        % Plot cell radius versus user density
        subplot(2, 1, 2);
        if (k==1)
            plot(userDensity_range, squeeze(cell_radius(:, j, k)), 'Color', colors(j, :), 'Marker', 'o', 'MarkerSize', 8);
        else
            
            plot(userDensity_range, squeeze(cell_radius(:, j, k)), 'Color', colors(j+3, :), 'LineWidth',3);
        end
        hold on;
        
        
    end
end

% Plot configurations
subplot(2, 1, 1);
xlabel('User Density (users/km^2)');
ylabel('Number of Cells');
title('Number of Cells vs User Density');

legend('Omni-directional (SIRmin = 14)', '120° Sectorization (SIRmin = 14)', ...
    '60° Sectorization (SIRmin = 14)', ...
    'Omni-directional (SIRmin = 19)', '120° Sectorization (SIRmin = 19)', ...
    '60° Sectorization (SIRmin = 19)');
grid on;

subplot(2, 1, 2);
xlabel('User Density (users/km^2)');
ylabel('Cell Radius (km)');
title('Cell Radius vs User Density');
legend('Omni-directional (SIRmin = 14)', '120° Sectorization (SIRmin = 14)' ...
    , '60° Sectorization (SIRmin = 14)', ...
    'Omni-directional (SIRmin = 19)', '120° Sectorization (SIRmin = 19)' ...
    , '60° Sectorization (SIRmin = 19)');
grid on;

% Save figure
% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf,  fullfile('plots',['Number_of_cells_and_Cell_radius_vs_' ...
        'User_Density_compare.png']));
end

%------------------------- PART B END -------------------------%



%------------------------- END -------------------------%