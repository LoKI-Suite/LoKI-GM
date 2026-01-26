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
    step = [];              % difference between the values at two consecutive nodes
    cellNumber = [];        % number of cells in the energy grid
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
      grid.step = gridProperties.maxEnergy/grid.cellNumber;
      grid.node = (0:grid.cellNumber)*grid.step;
      grid.cell = ((1:grid.cellNumber)-0.5)*grid.step; 
      if isfield(gridProperties, 'smartGrid')
        grid.isSmart = true;
        grid.minEedfDecay = gridProperties.smartGrid.minEedfDecay;
        grid.maxEedfDecay = gridProperties.smartGrid.maxEedfDecay;
        grid.updateFactor = gridProperties.smartGrid.updateFactor;
      end
    end
    
    function updateMaxValue(grid, maxValue)
      % resize the grid with a new maximum value
      grid.step = maxValue/grid.cellNumber;
      grid.node = (0:grid.cellNumber)*grid.step;
      grid.cell = ((1:grid.cellNumber)-0.5)*grid.step; 
      
      % broadcast change in the grid
      notify(grid, 'updatedMaxEnergy1'); % update interpolations of collision cross-sections
      notify(grid, 'updatedMaxEnergy2'); % update operators of the Boltzmann equation
    end

  end

end
