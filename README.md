# MATLAB Cell Planning Scripts

## Table of Contents

1. [Overview](#overview)
2. [Files](#files)
   - [methods_benchmark.m](#1-methods_benchmarkm)
   - [PartA.m](#2-parta)
   - [PartB.m](#3-partb)
   - [cellPlanning.m](#4-cellplanningm)
3. [CellPlanning Class](#cellplanning-class)
4. [Usage of the Class](#usage-of-the-class)
5. [Notes](#notes)
6. [How to Use This Repo](#how-to-use-this-repo)
7. [Requirements](#requirements)
8. [Lessons Learned](#lessons-learned)
9. [License](#license)

## Overview

This repository contains MATLAB scripts for cell planning simulations. The main components include methods for calculating cluster sizes, traffic intensity, and various parameters related to cell planning using different sectorization methods. The scripts also benchmark different methods for calculating traffic intensity.
This is done as a part of the course work for Wireless Networks course at Ain Shams University under supervision of Dr. Michael Gad

## Files

### 1. `methods_benchmark.m`

This script benchmarks various methods for calculating traffic intensity in a cell. It includes custom calculations as well as built-in MATLAB functions `fsolve` and `fzero`.

#### Main Functions:

- `cellPlanning.getClusterSize(io, SIRmin, n)`: Calculates the cluster size.
- `cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, method)`: Calculates the traffic intensity in a cell using various methods such as ('custom', 'fsolve', 'fzero', 'ahmedsaleh').

#### Example Execution:

```matlab
clc;
clear;
close all;

GOS = 1; SIRmin = 10;
io = 6; N_Sector = 1;
S = 340; n = 4; Au = 25e-3;

clusterSize = cellPlanning.getClusterSize(io, SIRmin, n);
AcellCustom = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'custom', true);
disp(['Traffic Intensity (Custom Function): ', num2str(AcellCustom), ' Erlangs']);
```

### 2. `PartA.m`

This script handles the first part of the cell planning process, including user input via GUI, calculation of cluster size, cell radius, number of cells, and received power vs. receiver distance plotting.

#### Main Functions:

- `gui_for_parameters()`: Prompts the user for inputs using a GUI.
- `cellPlanning.getClusterSize(io, SIRmin, n)`: Calculates the cluster size.
- `cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, method)`: Calculates traffic intensity in a cell using the custom method.
- `cellPlanning.calculateParameters(A_cell, Au, user_density, city_area)`: Calculates the cell radius and number of cells.

#### Example Execution:

```matlab
clear;
clc;
close all;
saveFig = false;

user_inputs = gui_for_parameters();
GOS = user_inputs.GOS;
city_area = user_inputs.city_area;
user_density = user_inputs.user_density;
SIRmin = user_inputs.SIRmin;
sectorizationMethod = user_inputs.sectorizationMethod;

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
        error('Invalid sectorization method.');
end

clusterSize = cellPlanning.getClusterSize(io, SIRmin, n);
A_cell = cellPlanning.calculateACell(S, clusterSize, N_Sector, GOS, 'custom');
[cell_Radius, Number_of_cells] = cellPlanning.calculateParameters(A_cell, Au, user_density, city_area);
```

### 3. `PartB.m`

This script handles the second part of the cell planning process, focusing on plotting cluster size vs. SIRmin for different sectorization methods and other related figures.

#### Main Functions:

- `cellPlanning.getClusterSize(io, SIRmin, n)`: Calculates the cluster size for different SIRmin values.

#### Example Execution:

```matlab
clc;
clear;
close all;
saveFig = false;

n = 4;
SIRmin_range = 1:0.01:30;
colors = linspecer(3);
clusterSizes = zeros(length(SIRmin_range), 3);

for i = 1:length(SIRmin_range)
    clusterSizes(i, 1) = cellPlanning.getClusterSize(6, SIRmin_range(i), n);
    clusterSizes(i, 2) = cellPlanning.getClusterSize(2, SIRmin_range(i), n);
    clusterSizes(i, 3) = cellPlanning.getClusterSize(1, SIRmin_range(i), n);
end

figure;
for i = 1:3
    plot(SIRmin_range, clusterSizes(:, i), 'Color', colors(i, :), 'LineWidth', 3);
    hold on;
end
xlabel('SIRmin (dB)');
ylabel('Cluster Size');
legend('Omni-directional', '120° Sectorization', '60° Sectorization');
title('Cluster Size vs SIRmin');
grid on;
```

### 4. `cellPlanning.m`

This Contains a class responsible for the whole cell planning and optimization process.

#### Main Functions:

- `getClusterSize(io, SIRmin, n)`: Calculates the cluster size.
- `calculateACell(S, clusterSize, N_Sector, GOS, method)`: Calculates the traffic intensity in a cell.
- `calculateParameters(A_cell, Au, user_density, city_area)`: Calculates the cell radius and number of cells.

# CellPlanning Class

The `CellPlanning` class provides methods for calculating various cell parameters used in cellular network simulations. This includes calculating cluster sizes, traffic intensity per cell, cell radius, and the number of cells needed to cover a specific area.

## Class Definition

### Class: `cellPlanning`

#### Methods

1. **`getClusterSize(io, minSIR, pathLossExponent)`**

   - Calculates the cluster size based on the input interference, minimum signal-to-interference ratio, and path loss exponent.
   - **Parameters:**
     - `io`: Input interference.
     - `minSIR`: Minimum signal-to-interference ratio.
     - `pathLossExponent`: Path loss exponent.
   - **Returns:**
     - `clusterSize`: The calculated cluster size.

2. **`calculateACell(S, clusterSize, numberSectors, GOS, MethodSelect, debugMode)`**

   - Calculates the traffic intensity per cell.
   - **Parameters:**
     - `S`: Total number of channels.
     - `clusterSize`: Cluster size.
     - `numberSectors`: Number of sectors.
     - `GOS`: Grade of Service.
     - `MethodSelect`: Method for calculation (`'custom'`, `'fzero'`, `'fsolve'`, `'ahmedsaleh'`).
     - `debugMode` (optional): Enable debug mode for detailed output.
   - **Returns:**
     - `A_cell`: The calculated traffic intensity per cell.

3. **`calculateParameters(A_cell, Au, userDensity, cityArea)`**

   - Calculates cell radius and the number of cells needed.
   - **Parameters:**
     - `A_cell`: Traffic intensity per cell.
     - `Au`: Traffic per user.
     - `userDensity`: User density.
     - `cityArea`: Total city area.
   - **Returns:**
     - `cellRadius`: The calculated cell radius.
     - `numberOfCells`: The estimated number of cells.

4. **`calculateASectorFsolve(gosRequired, C, debugMode)`**

   - Calculates the actual offered traffic using `fsolve` and inverse Erlang B method.
   - **Parameters:**
     - `gosRequired`: Grade of Service required.
     - `C`: Number of channels per sector.
     - `debugMode` (optional): Enable debug mode for detailed output.
   - **Returns:**
     - `A_sector`: The calculated traffic intensity per sector.

5. **`calculateASectorFzero(gosRequired, C, debugMode)`**

   - Calculates the actual offered traffic using `fzero` function.
   - **Parameters:**
     - `gosRequired`: Grade of Service required.
     - `C`: Number of channels per sector.
     - `debugMode` (optional): Enable debug mode for detailed output.
   - **Returns:**
     - `A_sector`: The calculated traffic intensity per sector.

6. **`calculateASectorCustom(gosRequired, C, debugMode)`**

   - Calculates the actual offered traffic using a custom Newton's optimization method.
   - **Parameters:**
     - `gosRequired`: Grade of Service required.
     - `C`: Number of channels per sector.
     - `debugMode` (optional): Enable debug mode for detailed output.
   - **Returns:**
     - `A_actual`: The calculated traffic intensity per sector.

7. **`ahmedsaleh(C, GOS, user_num)`**

   - Calculates the traffic intensity per sector using Ahmed Saleh's method.
   - **Parameters:**
     - `C`: Number of channels per sector.
     - `GOS`: Grade of Service required.
     - `user_num`: Number of users.
   - **Returns:**
     - `A_sec`: The calculated traffic intensity per sector.

8. **`erlangB(A, C)`**
   - Calculates blocking probability and its gradient for Erlang B formula.
   - **Parameters:**
     - `A`: Total offered traffic.
     - `C`: Number of channels.
   - **Returns:**
     - `GOS`: Grade of Service.
     - `gosDiff`: Gradient of Grade of Service.

## Usage of the Class

1. **Calculate Cluster Size:**

   ```matlab
   clusterSize = cellPlanning.getClusterSize(io, minSIR, pathLossExponent);
   ```

2. **Calculate Traffic Intensity Per Cell:**

   ```matlab
   A_cell = cellPlanning.calculateACell(S, clusterSize, numberSectors, GOS, MethodSelect, debugMode);
   ```

3. **Calculate Cell Parameters:**

   ```matlab
   [cellRadius, numberOfCells] = cellPlanning.calculateParameters(A_cell, Au, userDensity, cityArea);
   ```

4. **Calculate Traffic Intensity Per Sector:**

   ```matlab
   A_sector = cellPlanning.calculateASectorFsolve(gosRequired, C, debugMode);
   A_sector = cellPlanning.calculateASectorFzero(gosRequired, C, debugMode);
   A_sector = cellPlanning.calculateASectorCustom(gosRequired, C, debugMode);
   A_sec = cellPlanning.ahmedsaleh(C, GOS, user_num);
   ```

5. **Calculate Grade of Service:**
   ```matlab
   [GOS, gosDiff] = cellPlanning.erlangB(A, C);
   ```

## Notes

- Ensure that the `GOS` input in `calculateACell` and `calculateASector*` methods is between 0 and 100.
- Debug mode can be enabled in methods to print detailed computation steps and number of iterations.

## How to use this Repo

1. Clone the repository:

   ```bash
   git clone https://github.com/abdelrahmanali15/CircuitCruncher.git
   ```

2. Open MATLAB and navigate to the cloned repository.

3. Run the desired script:

4. Follow prompts and view results.

## Requirements

- MATLAB R2021a or later.
- Signal Processing Toolbox.

## Lessons Learned

All the Lessons learned for this project can be found at the pdf report file in the Repo. [Wireless Project Report](Wireless_Poject_Report_V04.pdf)

## License

This project is licensed under the MIT License.

---
