%% PART A

clear;
clc;
close all;
saveFig = false; %true for saving false to not save
%% %------------------------- Start -------------------------%

% Given Parameters
S = 340; n = 4; Au = 25e-3;
disp('No. of channels per cluster is 340 path loss exponent is  4 ,  user Traffic intensity is 25e-3' );

%------------------------- Test Cases -------------------------%

%     Test1
%     sectorizationMethod = '120deg';
%     GOS = 1; SIRmin = 19;user_density=1000;city_area=100;

%     Test2
%     sectorizationMethod = '120deg';
%     GOS = 1; SIRmin =14;user_density=1000;city_area=100;

%     Test3
%     sectorizationMethod = '60deg';
%     GOS = 0.2; SIRmin =25;user_density=1400;city_area=100;


%   Prompting user for inputs Using GUI
%------------------------- GUI Input -------------------------%
%
user_inputs = gui_for_parameters();
GOS = user_inputs.GOS;
city_area = user_inputs.city_area;
user_density = user_inputs.user_density;
SIRmin = user_inputs.SIRmin;
sectorizationMethod = user_inputs.sectorizationMethod;

% Sectorization method selection
switch lower(sectorizationMethod)
    case 'omni'
        io = 6;
        N_Sector = 1;
    case '120deg'
        io = 2;
        N_Sector = 3;
    case '60deg'
        io = 1;
        N_Sector = 6;
    otherwise
        error('Invalid sectorization method. Please enter either "omni", "120deg", or "60deg".');
end


%------------------------- PART A Main Logic -------------------------%

% Calculate cluster size
clusterSize = cellPlanning.getClusterSize(io, SIRmin, n);



% Calculate traffic intensity in a Cell using custom function
A_cell = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS,'custom');

% Calculate number of cells
[cell_Radius, Number_of_cells] = cellPlanning.calculateParameters(A_cell, Au, user_density, city_area);


% Hata model of losses
hB = 20; hM = 1.5; P_MS = -95;fc = 900;
% Get Base Station Power using sensitivity
CH = 0.8 + (1.1 * log10(fc) - 0.7) * hM - 1.56 * log10(fc);
Lu = 69.55 + 26.16 * log10(fc) - 13.82 * log10(hB) - CH + (44.9 - 6.55 * log10(hB)) * log10(cell_Radius);
P_BS = P_MS + Lu;

% P_MS with distance relation until the reuse distance (D)
reuse_distance = sqrt(3*clusterSize) * cell_Radius;
reciever_distance = 0:reuse_distance/1000:reuse_distance;
CH = 0.8 + (1.1 * log10(fc) - 0.7) * hM - 1.56 * log10(fc);
Lu = 69.55 + 26.16 * log10(fc) - 13.82 * log10(hB) - CH + (44.9 - 6.55 * log10(hB)) * log10(reciever_distance);
P_MS = P_BS - Lu;


% Plot
figure('Position', [100, 100, 1000, 700]);
plot(reciever_distance, P_MS, 'LineWidth', 2, 'DisplayName', 'Received Power');
title('Received Power vs. Receiver Distance');
xlabel('Receiver Distance (km)');
ylabel('Received Power (dBm)');
xlim([0, reuse_distance]);
grid on;
hold on;

% Place a Marker at P_MS is approximately -95 dBm or at the cell radius
plot(cell_Radius, -95, 'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r','HandleVisibility', 'off');

% Plot horizontal line at y = -95 dBm
line([0, reuse_distance], [-95, -95], 'Color', 'r', 'LineStyle', '--', 'DisplayName', 'SIRmin');

% Plot vertical line at x = cell radius
line([cell_Radius, cell_Radius], [-150, 0], 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'Cell Radius');

% Set y-axis ticks every 5 dBm
yticks(-150:5:max(P_MS(2:length(P_MS)))+10);

% Make grid lines lighter
ax = gca;
ax.GridAlpha = 0.1;

% Add legend
legend('Location', 'southeast');

% Add annotations
text(cell_Radius + 0.1, -100, ['Cell Radius: ', num2str(cell_Radius), ' km'], 'Color', 'k');
text(reuse_distance - 0.5, -95, '-95 dBm', 'Color', 'r');

% Write user inputs under the plot
text(0.02, 0.02, {'Inputs:- ' ...
    ['Sectorization Method: ', sectorizationMethod], ...
    ['GOS: ', num2str(GOS)], ...
    ['City Area: ', num2str(city_area)], ...
    ['User Density: ', num2str(user_density)], ...
    ['SIRmin: ', num2str(SIRmin)]}, ...
    'VerticalAlignment', 'bottom', 'Units', 'normalized');

% Create the folder 'plots' if it doesn't exist
if (saveFig)
    if ~exist('plots', 'dir')
        mkdir('plots');
    end
    saveas(gcf, fullfile('plots', 'Received_Power_vs_Receiver_Distance.png')); % Saves the figure as a PNG file
end

A_sector = A_cell / N_Sector;
% Displaying the outputs
disp('------------------------------------');
disp('Outputs:');
disp(['1) Cluster Size: ', num2str(clusterSize)]);
disp(['2) Number of cells: ', num2str(Number_of_cells), ' cells']);
disp(['3) Cell Radius: ', num2str(cell_Radius), ' kilometer']);
disp(['4) Traffic intensity per cell: ', num2str(A_cell), ' Erlang']);
disp(['5) Traffic intensity per sector: ', num2str(A_sector), ' Erlang']);
disp(['6) Base station transmitted power: ', num2str(P_BS),' dBm']);
disp('------------------------------------');

%------------------------- PART A END -------------------------%

%% GUI INPUT

function user_inputs = gui_for_parameters()
% Initialize flag for valid input
valid_input = false;

while ~valid_input
    % Create a dialog box to capture user inputs
    prompt = {'Grade of Service (in percentage):', ...
        'City area (in square kilometers):', ...
        'User density (users per square kilometer):', ...
        'Minimum Signal-to-Interference Ratio (SIR) (in dB):', ...
        'Sectorization method (omni/120deg/60deg):'};
    dlg_title = 'Input';
    num_lines = 1;
    default_input = {'', '', '', '', ''}; % Default input values
    
    % Display the dialog box
    user_input = inputdlg(prompt, dlg_title, num_lines, default_input);
    
    % Check if any field is empty or contains invalid data
    if isempty(user_input{1}) || isempty(user_input{2}) || isempty(user_input{3}) || isempty(user_input{4}) || isempty(user_input{5})
        % Display message and continue the loop
        h = msgbox('All fields are required. Please enter valid data.', 'Error', 'error');
        % Wait for the user to acknowledge the message box
        uiwait(h);
        
    else
        % Convert user inputs to appropriate data types
        GOS = str2double(user_input{1});
        city_area = str2double(user_input{2});
        user_density = str2double(user_input{3});
        SIRmin = str2double(user_input{4});
        sectorizationMethod = user_input{5}; % No conversion needed for string
        
        % Check if any conversion resulted in NaN (invalid data)
        if any(isnan([GOS, city_area, user_density, SIRmin]))
            % Display message and continue the loop
            h = msgbox(['Invalid data. Please enter numeric values for Grade of Service,' ...
                ' City area, User density, and Minimum SIR.'], 'Error', 'error');
            uiwait(h);
            
        else
            % If all inputs are valid, set flag to true to exit the loop
            valid_input = true;
        end
    end
end

% Store user inputs in a structure
user_inputs.GOS = GOS;
user_inputs.city_area = city_area;
user_inputs.user_density = user_density;
user_inputs.SIRmin = SIRmin;
user_inputs.sectorizationMethod = sectorizationMethod;

% Display the collected inputs
disp('User Inputs:');
disp(['Grade of Service: ', num2str(GOS)]);
disp(['City area: ', num2str(city_area)]);
disp(['User density: ', num2str(user_density)]);
disp(['Minimum SIR: ', num2str(SIRmin)]);
disp(['Sectorization method: ', sectorizationMethod]);
end



