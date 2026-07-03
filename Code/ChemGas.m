% LoKI-GM comprises two modules, that can run self-consistently coupled 
% or as standalone tools.
% 
% LoKI-B, which solves the space independent form of the two-term 
% electron Boltzmann equation (EBE) to calculate the isotropic and the 
% anisotropic parts of the electron distribution function, 
% and the associated electron macroscopic parameters. 
% LoKI-B applies to non-magnetised non-equilibrium LTPs, excited by 
% DC/HF electric fields or time-dependent (non-oscillatory) electric fields 
% from different gases or gas mixtures. 
% The tool uses a stationary description for DC fields, 
% a Fourier time-expansion description for HF fields, 
% and a time-dependent description for time-varying fields.
% 
% LoKI-C, which solves the system of zero-dimensional (volume average) 
% rate balance equations for the most relevant 
% charged and neutral species in the plasma. 
% LoKI-C receives as input data the kinetic schemes for the gas/plasma/
% surface system under study, via an intuitive csv-like input file, 
% and gives as output the particle densities of the different gas/plasma/
% surface species, the corresponding creation/destruction reaction rates, 
% and the reduced electric field (and any related quantity, such as 
% the discharge current or the discharge power-density).
% The tool uses several modules to describe the mechanisms 
% (collisional, radiative and transport) controlling the
% creation/destruction of species, namely various transport models 
% for the charged particles and for the neutral particles. 
% LoKI-C includes also a gas/plasma thermal model, for the self-consistent 
% calculation of the gas temperature, and supports multicomponent 
% mean-field microkinetic mesoscopic models to handle surface kinetics 
% in a fully coupled way with volume kinetics.
%
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% Copyright (C) 2026 L. L. Alves, A. Tejero-del-Caz, T. C. Dias, 
% A. Gonçalves, L. Marques, P. Pereira, N. Pinhão, C. D. Pintassilgo, 
% T. Silva, P. Viegas and V. Guerra
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

classdef ChemGas < Gas
  
  properties
    
    isVolumeSpecies;                  % boolean variable indicating if the gas (and its states) is a volume specie
    isSurfaceSpecies;                 % boolean variable indicating if the gas (and its states) is a surface specie
    
    reactionArray = Reaction.empty;   % array of handles to the reactions in wich this gas is involved
    
    eedfEquivalent = EedfGas.empty;   % handle to the equivalent EedfGas (in case it exists) 
    
  end
  
  events
    
  end
  
  methods
    
    function gas = ChemGas(gasName)
      
      persistent lastID;
      if isempty(lastID)
        lastID = 0;
      end
      lastID = lastID + 1;
      gas.ID = lastID;
      gas.name = gasName;
      gas.stateArray = ChemState.empty;
      
      % find if the species is in the volume or at the surface (keyword 'wall_' case insensitive)
      if length(gasName)>5 && strcmpi(gasName(1:5), 'wall_')
        gas.isVolumeSpecies = false;
        gas.isSurfaceSpecies = true;
      else
        gas.isVolumeSpecies = true;
        gas.isSurfaceSpecies = false;
      end
      
    end
    
    function linkWithElectronKineticsGas(chemGas, eedfGas)
      
      % save handle to the equivalent electron kinetic gas
      chemGas.eedfEquivalent = eedfGas;
      eedfGas.chemEquivalent = chemGas;
      
      % copy properties of the equivalent electron kinetic gas (avoiding 'ID' and 'stateArray' properties)
      eedfGasProperties = fields(eedfGas)';
      chemGasProperties = fields(chemGas)';
      for property = eedfGasProperties
        if ~strcmp(property{1},'ID') && ~strcmp(property{1}, 'stateArray') && any(strcmp(property, chemGasProperties))
          chemGas.(property{1}) = eedfGas.(property{1});
        end
      end
      
    end
    
    function updateElectronKineticsEquivalentPopulations(gas)
    % updateElectronKineticsEquivalentPopulations is a function that updates the populations of the equivalent
    % gases that are used to solve the electron kinetics. The function avoids the electron kinetics gases for which
    % the user did not specify a cross section set (dummy gases). Because of the same reason the function also
    % avoids to update electronic or ionic states for which the user did not specify a non zero population.
      
      % avoid gases that does not have an eedf equivalent or gases whos states are not target of e-impact collisions
      if ~isempty(gas.eedfEquivalent) && ~isempty(gas.eedfEquivalent.collisionArray)
        % loop over all the states of a gas
        for state = gas.stateArray
          % avoid states without an eedf equivalent state and states who are not target of any e-impact collision
          if ~isempty(state.eedfEquivalent) && ~isempty(state.eedfEquivalent.collisionArray)
            state.eedfEquivalent.population = state.population;
          end
        end
      end
      
    end
    
    function checkThermalModelData(gas)
    % checkThermalModelData checks for the data needed for the thermal module to be activated to be present in the gas
    % properties (heat capacity and heat conductivity)
      
      % check for the definition of the gas mass
      if isempty(gas.mass)
        error(['A value for the mass is not found for gas %s.\n' ...
          'Thermal model can not be activated without it.\n' ...
          'Please check your setup file.'], gas.name);
      end
      
      % check for the definition of the heat conductivity
      if isempty(gas.thermalConductivity)
        error(['A value for the thermal conductivity is not found for gas %s.\n' ...
          'Thermal model can not be activated without it.\n' ...
          'Please check your setup file.'], gas.name);
      end
      
      % check for the definition of the heat capacity
      if isempty(gas.heatCapacity)
        error(['A value for the heat capacity is not found for gas %s.\n' ...
          'Thermal model can not be activated without it.\n' ...
          'Please check your setup file.'], gas.name);
      end
      
    end

    function checkInFlowPopulationNorms(gas)
    % checkInFlowPopulationNorms checks for the inflow population of the different states of the gas to be properly normalised,
    % i. e. the inflow populations of all sibling states should add to one.

      % avoid gases not present in the mixture
      if gas.inFlowFraction == 0
        return;
      end

      % check norm of electronic/ionic states
      gasNorm = 0;
      electronicStatesToBeChecked = true;
      ionicStatesToBeChecked = true;
      for state = gas.stateArray
        % norm of electronic states
        if strcmp(state.type, 'ele') && electronicStatesToBeChecked
          for eleState = [state state.siblingArray]
            if eleState.inFlowPopulation ~= 0
              gasNorm = gasNorm + eleState.inFlowPopulation;
              % check norm of vibrational states (if they exist)
              if ~isempty(eleState.childArray)
                vibNorm = 0;
                for vibState = eleState.childArray
                  if vibState.inFlowPopulation ~= 0
                    vibNorm = vibNorm + vibState.inFlowPopulation;
                    % check norm of rotational states (if they exist)
                    if ~isempty(vibState.childArray)
                      rotNorm = 0;
                      for rotState = vibState.childArray
                        if rotState.inFlowPopulation ~= 0
                          rotNorm = rotNorm + rotState.inFlowPopulation;
                        end
                      end
                      if abs(rotNorm-1) > 10*eps(1)
                        stateDistName = vibState.name;
                        stateDistName = [stateDistName(1:end-1) ',J=*)'];
                        error('Rotational inFlow distribution %s is not properly normalised. (Error = %e)\n', ...
                          stateDistName, rotNorm-1);
                      end
                    end
                  end
                end
                if abs(vibNorm-1) > 10*eps(1)
                  stateDistName = eleState.name;
                  stateDistName = [stateDistName(1:end-1) ',v=*)'];
                  error('Vibrational inFlow distribution %s is not properly normalised. (Error = %e)\n', ...
                    stateDistName, vibNorm-1);
                end
              end
            end
          end
          electronicStatesToBeChecked = false;
        end
        if strcmp(state.type, 'ion') && ionicStatesToBeChecked
          for ionState = [state state.siblingArray]
            if ionState.inFlowPopulation ~= 0
              gasNorm = gasNorm + ionState.inFlowPopulation;
            end
          end
          ionicStatesToBeChecked = false;
        end
        if ~electronicStatesToBeChecked && ~ionicStatesToBeChecked
          break;
        end
      end
      if abs(gasNorm-1) > 10*eps(1)
        stateDistName = [gas.name '(*)'];
        error('Electronic/ionic inFlow distribution %s is not properly normalised. (Error = %e)\n', stateDistName, gasNorm-1);
      end

    end      

    function evaluateInFlowDensities(gas)
    % evaluateInFlowDensities is a function that evaluates the inflow relative densities of all the states of the gas from their
    % populations and the gas fraction

      % loop over all the states of the gas
      for state = gas.stateArray
        state.evaluateInFlowRelDensity();
      end

    end      
        
  end
  
  methods (Static)
    
    function updateElectronKineticsFractions(gasArray)
    % updateElectronKineticsFractions is a function that updates (properly reescaling) the fractions of the
    % equivalent electron kinetics gases. The function avoids gases that do not have an electron kinetics
    % equivalent with a cross section set (dummy gases).
      
      % find gases to update and evaluate the fractions norm
      gasesToUpdate = [];
      fractionNorm = 0;
      for i = 1:length(gasArray)
        if ~isempty(gasArray(i).eedfEquivalent) && ~isempty(gasArray(i).eedfEquivalent.collisionArray)
          gasesToUpdate(end+1) = i;
          fractionNorm = fractionNorm + gasArray(i).fraction;
        end
      end
      
      % update gases
      for i = gasesToUpdate
        gasArray(i).eedfEquivalent.fraction = gasArray(i).fraction/fractionNorm;
      end
      
    end

    function checkInFlowFractionNorm(gasArray)
    % checkInFlowFractionNorm checks for the inflow fractions of the different gases in gasArray to be properly normalized. 
    % Moreover, it assures that only volume species have an inflow different than zero

      norm = 0;
      for gas = gasArray
        if gas.isSurfaceSpecies && gas.inFlowFraction ~= 0
          error('The gas "%s" cannot have an inFlowFraction different than zero, as it is a surface species!',gas.name);
        end  
        norm = norm + gas.inFlowFraction;
      end
      if abs(norm-1) > 10*eps(1)
        error('inFlowFractions are not properly normalised (Error = %e).\nPlease, check input file.', norm-1);
      end

    end     
    
  end
  
end