classdef cellPlanning
    % CellPlanning: A class for calculating cell parameters.
    %   This class provides methods to calculate various cell parameters
    %   used in cellular network simulations.
    
    methods (Static)
        
        function clusterSize = getClusterSize(io, minSIR, pathLossExponent)
            % Calculate cluster size.
            SIR_ratio = 10 .^ (minSIR / 10);
            originalClusterSize = (1/3) * (((io * SIR_ratio).^(1/pathLossExponent))+1).^2;
            
            % Compute cluster size values
            maxIndex = ceil(sqrt(originalClusterSize));
            [I, K] = meshgrid(0:maxIndex);
            newClusterSize = I.^2 + K.^2 + I.*K;
            
            % Find the minimum cluster size larger than the original cluster size
            clusterSize = min(newClusterSize(newClusterSize >= originalClusterSize));
        end
        
        function A_cell = calculateACell(S, clusterSize, numberSectors, GOS, MethodSelect,debugMode)
            % Calculate the traffic intensity per cell.
            
            if any(GOS < 0) || any(GOS > 100)
                error('Invalid input. Grade of Service should be between 0 and 100.');
            end
            
            % Set default value for debugMode if not provided
            if nargin < 6
                debugMode = false;
            end
            
            % Channels per cell
            k = S / clusterSize;
            
            % Channels per sector
            c = floor(k / numberSectors);
            
            if any(c < 1)
                error('Invalid SIR: Not enough channels for each sector.');
            end
            
            % Find traffic intensity per sector
            if strcmp(MethodSelect,'custom')
                % Calculate using our own implemented function
                A_sector = cellPlanning.calculateASectorCustom(GOS, c ,debugMode);
            elseif strcmp(MethodSelect,'fzero')
                % Calculate using fzero() function
                A_sector = cellPlanning.calculateASectorFzero(GOS, c ,debugMode);
            elseif strcmp(MethodSelect,'fsolve')
                % Calculate using fsolve() function
                A_sector = cellPlanning.calculateASectorFsolve(GOS, c ,debugMode);
            elseif strcmp(MethodSelect,'ahmedsaleh')
                A_sector = cellPlanning.ahmedsaleh(c, GOS ,500);
            else
                error('Invalid MethodSelect: Please choose "custom" or "fzero" or "fsolve".');
            end
            
            % Find traffic intensity per cell
            A_cell = A_sector .* numberSectors;
        end
        
        function [cellRadius, numberOfCells] = calculateParameters(A_cell, Au, userDensity, cityArea)
            % Calculate cell radius and number of cells.
            
            % Calculate the number of users per cell
            usersCell = floor(A_cell / Au);
            
            % Calculate the area of each cell
            areaCell = usersCell ./ userDensity;
            
            % Calculate the radius of each cell
            cellRadius = sqrt(2 * areaCell / (sqrt(3) * 3));
            
            % Estimate the number of cells needed to cover the entire city area
            numberOfCells = ceil(cityArea ./ areaCell);
        end
        
        function A_sector = calculateASectorFsolve(gosRequired, C,debugMode)
            % Calculate the actual offered traffic using fsolve and inverse ErlangB method function.
            
            gosRequired = gosRequired / 100;
            
            objectiveFcn = @(A) cellPlanning.erlangB(A, C) - gosRequired;
            initialGuess = C; % or any other suitable initial guess
            % options = optimoptions('fsolve', 'Display', 'final','TolFun',1e-9,'TolX',1e-9);
            options = optimoptions('fsolve','TolFun',1e-14,'TolX',1e-14,'Display', 'off');
            [A_sector,~,~,out] = fsolve(objectiveFcn, initialGuess, options);
             
            % Print number of iterations for debug mode
            iterations = out.iterations;
            if nargin > 2 && debugMode
               disp(['Number of iterations (Fsolve Function): ', num2str(iterations)]);
            end
        end
        
        function A_sector = calculateASectorFzero(gosRequired, C, debugMode)
            % Calculate the actual offered traffic using fzero function and normal erlangB.

            gosRequired = gosRequired / 100;
        
            objectiveFcn = @(A) (((A^C) / factorial(C)) - gosRequired * (sum(A.^(0:C) ./ factorial(0:C))));
            initialGuess = C; % or any other suitable initial guess
            options = optimset(@fzero); % default options
            % options = optimset('Display', 'final');
            [A_sector,~,~,out] = fzero(objectiveFcn, initialGuess, options);
            
            % Print number of iterations for debug mode
            iterations = out.iterations;
            if nargin > 2 && debugMode
                disp(['Number of iterations (Fzero Function): ', num2str(iterations)]);
            end
        end

        
        
        % This is our own implemented method to calculate fzero
        function A_actual = calculateASectorCustom(gosRequired, C, debugMode)
            % Calculate the actual offered traffic using custom implemented newton's optimization method and inverse Erlang B computaion.
            
            
            gosRequired = gosRequired / 100; % convert percentage to number
            A_guess = C; % Initial guess
            maxIterations = 200; % Maximum iterations
            tolerance = 1e-6; % Convergence criteria
            
            for i = 1:maxIterations
                [gosCalculated, gosPrime] = cellPlanning.erlangB(A_guess, C);
                gosError = gosCalculated - gosRequired;
                
                % Check convergence criteria
                if abs(gosError) < tolerance
                    break
                end
                
                % Update guess
                A_next = A_guess - (gosError / gosPrime);
                
                % Handle negative guesses
                if A_next < 0
                    A_next = C * rand();
                end
                
                A_guess = A_next;
            end
            
            if i == maxIterations
                warning('calculateASector:MaxIterationsReached', ...
                    'Maximum iterations reached. Solution may not be accurate.');
                % Try solving using matlab fzero function
                A_actual = cellPlanning.calculateAFsolve(gosRequired, A_guess);
            else
                % Return actual offered traffic calculated by this function
                A_actual = A_guess;

                 % Print number of iterations for debug mode
                if nargin > 2 && debugMode
                    disp(['Number of iterations (Custom Function): ', num2str(i)]);
                end
                
            end
        end
        function A_sec = ahmedsaleh(C, GOS, user_num)
            % Find traffic intensity per sector
            %max_intensity = user_num^(1/(c+1));
            max_intensity = user_num;
            % Calculate maximum offered traffic using Erlang B formula
            %max_intensity = (GOS * C) / (1 - GOS);
            A = 0:0.01:max_intensity;
            num = A .^ C / factorial(C);
            den = zeros(1, length(A));
            for k = 0:C
                den = den + A .^ k / factorial(k);
            end
            p = num ./ den;
            A_sec = zeros(1, length(GOS));
            for gosCounter = 1:length(GOS)
                A_sec(gosCounter) = mean(A(abs(p - GOS(gosCounter)) / GOS(gosCounter) <= 0.01));
            end
        end
        function [GOS, gosDiff] = erlangB(A, C)
            % Calculate blocking probability and its gradient for Erlang B formula.
            % This uses Inverse Erlang B which is better in computation and support large
            % factorial (number of channels 100 and more)
            
            % Input validation
            if A <= 0 || C < 0 || mod(C, 1) ~= 0
                error('Invalid input. Total offered traffic should be positive and channel number should be a non-negative integer.');
            end
            
            inv_b = 0;
            inv_b_prime = 0;
            
            % Calculate inverse erlang_b function and its gradient
            for i = 0:C
                iter = 0 : C-i-1;
                product = prod((C - iter)./A, 2);
                inv_b = inv_b + product;
                inv_b_prime = inv_b_prime + product * (i - C) / A;
            end
            
            % Calculate GOS and its gradient
            if inv_b ~= 0
                GOS = 1 / inv_b;
                gosDiff = -inv_b_prime / (inv_b^2);
            else
                % Handle division by zero
                GOS = Inf;
                gosDiff = NaN;
            end
        end
    end
end
