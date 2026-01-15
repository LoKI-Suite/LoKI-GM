% LoKI-B solves a time and space independent form of the two-term 
% electron Boltzmann equation (EBE), for non-magnetised non-equilibrium 
% low-temperature plasmas excited by DC/HF electric fields from 
% different gases or gas mixtures.
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

classdef Grid < handle
  %Grid Class that defines a grid and its interactions with other objects
  %   Class that defines a grid and its interactions with other objects. A
  %   grid object stores three properties: the values of the grid at node 
  %   positions, the values of the grid at cell positions (middle point) and 
  %   the step of the grid. 
  %
  %   Node values -> |     |     |     |     |     |     |     |
  %   Grid        -> o-----o-----o-----o-----o-----o-----o-----o
  %   Cell values ->    |     |     |     |     |     |     |
  
  properties

    node = [];              % values of the grid at node positions
    cell = [];              % values of the grid at cell positions (i.e. between two consecutive nodes)
    step = [];              % difference between the values at two consecutive nodes (scalar for uniform, first element for variable)
    energyStep = [];        % step values between consecutive nodes (array for variable grid, constant array for uniform grid)
    cellNumber = [];        % number of cells in the energy grid
    maxEnergy = [];         % maximum energy of the grid (in eV)
    variableGrid = false;   % whether to use variable grid (default: false)
    ratio = 1.1;            % ratio 'a' for geometric progression (used only for variable grid)
    firstEnergyStep = [];   % first energy step (required if variableGrid is true)
    isSmart = false;        % smart properties of the energy grid (deactivated by default)
    minEedfDecay = [];      % minimum number of decades of decay for the EEDF
    maxEedfDecay = [];      % maximum number of decades of decay for the EEDF
    updateFactor = [];      % percentage factor to update the maximum value of the energy grid
    
  end
  
  events
    updatedMaxEnergy1
    updatedMaxEnergy2
  end

  methods

    function grid = Grid(gridProperties)
      grid.cellNumber = gridProperties.cellNumber;
      grid.maxEnergy = gridProperties.maxEnergy;
      
      % Check if variableGrid is provided and set it
      if isfield(gridProperties, 'variableGrid')
        grid.variableGrid = gridProperties.variableGrid;
      end
      
      % If variableGrid is true, firstEnergyStep must be provided
      if grid.variableGrid
        if ~isfield(gridProperties, 'firstEnergyStep')
          error('firstEnergyStep must be provided when variableGrid is true.');
        end
        grid.firstEnergyStep = gridProperties.firstEnergyStep;
        
        % Perform self-diagnostic test
        grid.selfDiagnosticTest();
        % Compute the progression ratio 'a'
        grid.computeProgressionRatio();
        % Compute the nodes using geometric progression
        grid.computeGeometricGrid();
      else
        % Compute the nodes using constant energy step (uniform grid)
        grid.computeUniformGrid();
      end
      
      % Set smart grid properties if provided
      if isfield(gridProperties, 'smartGrid')
        grid.isSmart = true;
        grid.minEedfDecay = gridProperties.smartGrid.minEedfDecay;
        grid.maxEedfDecay = gridProperties.smartGrid.maxEedfDecay;
        grid.updateFactor = gridProperties.smartGrid.updateFactor;
      end
    end
    
    function computeProgressionRatio(grid)
      % Define the function F(a) for solving progression ratio 'a'
      F = @(a) exp(log((grid.maxEnergy / grid.firstEnergyStep) * (a - 1) + 1) / (grid.cellNumber)) - a;
      
      % Use fzero to solve for progression ratio 'a'
      options = optimset('Display', 'iter'); % show iterations
      solution = fzero(F, [1+1E-15, 3], options); %default: 1+1E-15
      
      grid.ratio = solution;
      
      % Check if the solution is valid
      if grid.ratio <= 1
        error('Invalid progression ratio "a". Ensure that maxEnergy and firstEnergyStep are consistent.');
      end
    end
    
    function selfDiagnosticTest(grid)
      % Self-diagnostic test to ensure the progression ratio 'a' is valid
      if grid.variableGrid
        % Check if the total sum of steps would exceed maxEnergy
        if grid.firstEnergyStep * grid.cellNumber >= grid.maxEnergy
          error(['Invalid progression ratio "a" (>= 1). Ensure that maxEnergy and firstEnergyStep are consistent. ' ...
                 'Consider using a smaller value for firstEnergyStep.']); 
        end
      end
    end

    function computeGeometricGrid(grid)
      % Compute node and cell values for a non-uniform geometric progression grid
      grid.node = zeros(1, grid.cellNumber + 1); % Initialize node array
      grid.energyStep = zeros(1, grid.cellNumber); % Initialize step sizes
      
      % Define the first step size (u1)
      u1 = grid.firstEnergyStep;
      
      % Compute nodes and steps using geometric progression
      grid.node(1) = 0; % Starting energy
      for n = 1:grid.cellNumber
        grid.energyStep(n) = u1 * grid.ratio^(n - 1);
        grid.node(n + 1) = grid.node(n) + grid.energyStep(n);
      end
      
      % % new
      % % Adjust last node to exactly match maxEnergy (for consistency)
      % if abs(grid.node(end) - grid.maxEnergy) > 1e-10
      %   grid.node(end) = grid.maxEnergy;
      %   grid.energyStep(end) = grid.maxEnergy - grid.node(end-1);
      % end

      % Compute cell values as midpoints between nodes
      grid.cell = 0.5 * (grid.node(1:end-1) + grid.node(2:end));
      
      % new
      % Set step property to first energy step for backward compatibility
      grid.step = grid.energyStep(1);
    end

    function computeUniformGrid(grid)
      
      % % new  
      % % Compute node and cell values for a uniform energy grid
      % grid.step = grid.maxEnergy / grid.cellNumber; % constant step size
      % grid.node = (0:grid.cellNumber) * grid.step; % Uniformly spaced nodes
      % grid.energyStep = ones(1, grid.cellNumber) * grid.step; % Constant step sizes (as array)
      % 
      % % Compute cell values as midpoints between nodes
      % grid.cell = ((1:grid.cellNumber) - 0.5) * grid.step;
      % Compute node and cell values for a uniform energy grid

      grid.node = linspace(0, grid.maxEnergy, grid.cellNumber + 1); % Uniformly spaced nodes
      grid.energyStep = diff(grid.node); % Constant step sizes
            
      % Compute cell values as midpoints between nodes
      grid.cell = 0.5 * (grid.node(1:end-1) + grid.node(2:end));

      % new
      % Set step property to first energy step for backward compatibility
      grid.step = grid.energyStep(1);
      
    end
    
    function updateMaxValue(grid, maxValue)
      % resize the grid with a new maximum value
      grid.maxEnergy = maxValue;
      if grid.variableGrid
        grid.computeProgressionRatio();
        grid.computeGeometricGrid();
      else
        grid.computeUniformGrid();
      end
      
      % broadcast change in the grid
      notify(grid, 'updatedMaxEnergy1'); % update interpolations of collision cross-sections
      notify(grid, 'updatedMaxEnergy2'); % update operators of the Boltzmann equation
    end

    function cellNumber = findCellNumber(grid, energy)
      % Calculate the cell number for a given energy value
      if energy < grid.node(1) || energy > grid.node(end)
        error('Energy value is outside the grid range.');
      end
      
      if grid.variableGrid
        % For variable grid, use the geometric progression formula
        u1 = grid.energyStep(1);
        a = grid.ratio;
        cellNumber = floor(log(energy / u1 * (a - 1) + 1) / log(a) + 1);
        
        % % new
        % % Ensure cellNumber is within bounds
        % cellNumber = min(cellNumber, grid.cellNumber);

      else
        % For uniform grid, use linear interpolation
        cellNumber = floor(energy / (grid.maxEnergy / grid.cellNumber));
        
        % % new
        % % Ensure cellNumber is within bounds
        % cellNumber = min(cellNumber, grid.cellNumber);

      end
    end

  end

end
