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

classdef Gas < handle
  %Gas Class that stores the information of a certain gas
  %   Class that stores the information of a certain gas present in the
  %   mixture of gases included in the discharge. The electron collisions
  %   and the populations of the different states are to be used in a
  %   Boltzmann solver to obtain the EEDF.

  properties

    ID = -1;                             % ID that identifies the gas in the gas array
    name = '';                           % name of the gas
                                         
    mass = [];                           % mass of the gas
    massFunc = [];                       % handle to function that evaluates the mass of the gas
    massParams = {};                     % cell array of parameters needed by massFunc
                                         
    harmonicFrequency = [];              % harmonic oscillator frequency (molecules)
    harmonicFrequencyFunc = [];          % handle to function that evaluates harmonic oscillator frequency of the gas
    harmonicFrequencyParams = {};        % cell array of parameters needed by harmonicFrequencyFunc 
                                         
    anharmonicFrequency = [];            % anharmonic oscillator frequency (molecules)
    anharmonicFrequencyFunc = [];        % handle to function that evaluates anharmonic oscillator frequency of the gas
    anharmonicFrequencyParams = {};      % cell array of parameters needed by anharmonicFrequencyFunc 
                                         
    rotationalConstant = [];             % rotational constant (molecules)
    rotationalConstantFunc = [];         % handle to function that evaluates rotational constant of the gas
    rotationalConstantParams = {};       % cell array of parameters needed by rotationalConstantFunc 
                                         
    lennardJonesDistance = [];           % sigma parameter of Lennard-Jones potential
    lennardJonesDistanceFunc = [];       % handle to function that evaluates Lennard-Jones distance of the gas
    lennardJonesDistanceParams = {};     % cell array of parameters needed by lennardJonesDistanceFunc
                                         
    lennardJonesDepth = [];              % epsilon parameter of Lennard-Jones potential
    lennardJonesDepthFunc = [];          % handle to function that evaluates Lennard-Jones depth of the gas
    lennardJonesDepthParams = {};        % cell array of parameters needed by lennardJonesDepthFunc
                                         
    electricDipolarMoment = [];          % electric dipolar moment (molecules)
    electricDipolarMomentFunc = [];      % handle to function that evaluates dipolar moment of the gas
    electricDipolarMomentParams = {};    % cell array of parameters needed by electricDipolarMomentFunc
                                         
    electricQuadrupoleMoment = [];       % electric cuadrupole moment (molecules)
    electricQuadrupoleMomentFunc = [];   % handle to function that evaluates quadrupolar moment of the gas
    electricQuadrupoleMomentParams = {}; % cell array of parameters needed by electricQuadrupoleMomentFunc
    
    polarizability = [];                 % polarizability (molecules)
    polarizabilityFunc = [];             % handle to function that evaluates polarizability of the gas
    polarizabilityParams = {};           % cell array of parameters needed by polarizabilityFunc
                                         
    fraction = 0;                        % fraction in the gas mixture
    fractionFunc = [];                   % handle to function that evaluates fraction of the gas
    fractionParams = {};                 % cell array of parameters needed by fractionFunc
    
    heatCapacity = [];                   % constant pressure heat capacity 
    heatCapacityFunc = [];               % handle to function that evaluates the constant pressure heat capacity 
    heatCapacityParams = {};             % cell array of parameters needed by heatCapacity
                                         
    thermalConductivity = [];            % thermal conductivity
    thermalConductivityFunc = [];        % handle to function that evaluates the thermal conductivity
    thermalConductivityParams = {};      % cell array of parameters needed by thermalConductivityFunc

    inFlowFraction = 0;                  % inflow fraction in the gas mixture
    inFlowFractionFunc = [];             % handle to function that evaluates inflow fraction of the gas
    inFlowFractionParams = {};           % cell array of parameters needed by inFlowFractionFunc
                                         
    stateArray;                          % State array with all the states of the gas (initialized in the subclass)

  end
  
  events
    
  end

  methods (Access = public)

    function disp(gas)
      for field = fieldnames(gas)'
        if ~isempty(gas.(field{1}))
          if isnumeric(gas.(field{1}))
            fprintf('%s: ', field{1});
            fprintf('%g ', gas.(field{1}));
            fprintf('\n');
          elseif ischar(gas.(field{1}))
            fprintf('%s: %s\n', field{1}, gas.(field{1}));
          end
        end
      end
    end

    function atoms = getAtoms(gas)
    % getAtoms returns the atomic composition of the gas
      
      % avoid "gases without atoms" (empty sites for surface model)
      if strcmp(gas.name, "wall_EmptySite")
        atoms = [];
        return
      end

      % initialize structure with atom names (strings) and count (integers)
      atoms = struct('name', [], 'count', []);

      % parse gas name to obtain atoms and their count
      atomsInfo = regexp(gas.name, '([A-Z][a-z]*)(\d*)', 'tokens');
      for i = 1:length(atomsInfo)
        atoms(end+1).name = atomsInfo{i}{1};
        if isempty(atomsInfo{i}{2})
          atoms(end).count = 1;
        else
          atoms(end).count = str2double(atomsInfo{i}{2});
        end
      end

      % remove first element (empty) and return
      atoms = atoms(2:end);
    end

    function dispStates(gas)
      fprintf('%s states:\n', gas.name);
      for i = 1:length(gas.stateArray)
        fprintf('\t %d.- %s\n', i, gas.stateArray(i).name);
      end
    end

    function dispFamilyTree(gas)
      fprintf('%s states: (family tree)\n', gas.name);
      fprintf('%s\n', gas.name);
      for eleState = gas.stateArray
        if strcmp(eleState.type, 'ele')
          fprintf('\t|%s\n', eleState.name);
          for vibState = eleState.childArray
            fprintf('\t|\t|%s\n', vibState.name);
            for rotState = vibState.childArray
              fprintf('\t|\t|\t|%s\n', rotState.name);
            end
          end
        end
      end
      for ionState = gas.stateArray
        if strcmp(ionState.type, 'ion')
          fprintf('%s\n', ionState.name);
        end
      end
    end

    function checkPopulationNorms(gas)
    % checkPopulationNorms checks for the population of the different states of the gas to be properly normalised,
    % i. e. the populations of all sibling states should add to one.

      % avoid gases not present in the mixture
      if gas.fraction == 0
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
            if eleState.population ~= 0
              gasNorm = gasNorm + eleState.population;
              % check norm of vibrational states (if they exist)
              if ~isempty(eleState.childArray)
                vibNorm = 0;
                for vibState = eleState.childArray
                  if vibState.population ~= 0
                    vibNorm = vibNorm + vibState.population;
                    % check norm of rotational states (if they exist)
                    if ~isempty(vibState.childArray)
                      rotNorm = 0;
                      for rotState = vibState.childArray
                        if rotState.population ~= 0
                          rotNorm = rotNorm + rotState.population;
                        end
                      end
                      if abs(rotNorm-1) > 10*eps(1)
                        stateDistName = vibState.name;
                        stateDistName = [stateDistName(1:end-1) ',J=*)'];
                        error('Rotational distribution %s is not properly normalised. (Error = %e)\n', ...
                          stateDistName, rotNorm-1);
                      end
                    end
                  end
                end
                if abs(vibNorm-1) > 10*eps(1)
                  stateDistName = eleState.name;
                  stateDistName = [stateDistName(1:end-1) ',v=*)'];
                  error('Vibrational distribution %s is not properly normalised. (Error = %e)\n', ...
                    stateDistName, vibNorm-1);
                end
              end
            end
          end
          electronicStatesToBeChecked = false;
        end
        if strcmp(state.type, 'ion') && ionicStatesToBeChecked
          for ionState = [state state.siblingArray]
            if ionState.population ~= 0
              gasNorm = gasNorm + ionState.population;
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
        error('Electronic/ionic distribution %s is not properly normalised. (Error = %e)\n', stateDistName, gasNorm-1);
      end

    end  

    function renormalizeWithDensities(gas)
    % renormalizeWithDensities is a function that obtain new "populations" of the states of a gas once their
    % "densities" have been updated. Note that the "populations" of the states are just the "densities" normalized
    % to the sibling states and that the "densities" are the absolute densites (in m^-3) normalized to the total
    % gas density.

      % loop over electronic/ionic states to obtain the new gas fraction (sum of electronic/ionic densities)
      gasFraction = 0;
      eleStatesNeedToBeChecked = true;
      ionStatesNeedToBeChecked = true;
      for state = gas.stateArray
        % check electronic states
        if strcmp(state.type, 'ele') && eleStatesNeedToBeChecked
          for eleState = [state.siblingArray state]
            gasFraction = gasFraction + eleState.density;
          end
          eleStatesNeedToBeChecked = false;
        end
        % check ionic states
        if strcmp(state.type, 'ion') && ionStatesNeedToBeChecked
          for ionState = [state.siblingArray state]
            gasFraction = gasFraction + ionState.density;
          end
          ionStatesNeedToBeChecked = false;
        end
        % if everything is checked exit the loop and save the new gas fraction
        if ~eleStatesNeedToBeChecked && ~ionStatesNeedToBeChecked
          gas.fraction = gasFraction;
          break;
        end
      end

      % avoid NaNs for "empty" gases
      if gasFraction == 0
        return
      end

      % loop over electronic and ionic states to obtain the new populations (for states with a chemistry equivalent)
      eleStatesNeedToBeChecked = true;
      ionStatesNeedToBeChecked = true;
      for state = gas.stateArray
        % check electronic states
        if strcmp(state.type, 'ele') && eleStatesNeedToBeChecked
          for eleState = [state.siblingArray state]
            if isempty(eleState.chemEquivalent)
              continue
            end
            eleState.population = eleState.density/gasFraction;
            % obtain the new populations of the vibrational distribution (in case it exists)
            if ~isempty(eleState.childArray)
              for vibState = eleState.childArray
                if isempty(vibState.chemEquivalent)
                  continue
                end
                vibState.population = vibState.density/eleState.density;
                % obtain the new populations of the rotational distribution (in case it exists)
                if ~isempty(vibState.childArray)
                  for rotState = vibState.childArray
                    if isempty(rotState.chemEquivalent)
                      continue
                    end
                    rotState.population = rotState.density/vibState.density;
                  end
                end
              end
            end
          end
          eleStatesNeedToBeChecked = false;
        end
        % check ionic states
        if strcmp(state.type, 'ion') && ionStatesNeedToBeChecked
          for ionState = [state.siblingArray state]
            if isempty(ionState.chemEquivalent)
              continue
            end
            ionState.population = ionState.density/gasFraction;
          end
          ionStatesNeedToBeChecked = false;
        end
        % exit the loop if everything is checked
        if ~eleStatesNeedToBeChecked && ~ionStatesNeedToBeChecked
          break;
        end
      end

      % evaluate density of states with the updated populations
      gas.evaluateDensities();

    end

    function evaluateDensities(gas)
    % evaluateDensities is a function that evaluates the densities of all the states of the gas from their
    % populations and the gas fraction. Note that the densities of the states are normalized to the total gas
    % density.

      % loop over all the states of the gas
      for state = gas.stateArray
        state.evaluateDensity();
      end

    end  
    
    function massLocal = evaluateMass(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.massFunc)
        % return fixed parameter
        massLocal = gas.mass;
      else
        % call function to evaluate the value of the property
        massLocal = gas.massFunc(gas, gas.massParams, workCond);
        % save local value in object properties
        gas.mass = massLocal;
      end
      
    end
    
    function harmonicFrequencyLocal = evaluateHarmonicFrequency(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.harmonicFrequencyFunc)
        % return fixed parameter
        harmonicFrequencyLocal = gas.harmonicFrequency;
      else
        % call function to evaluate the value of the property
        harmonicFrequencyLocal = gas.harmonicFrequencyFunc(gas, gas.harmonicFrequencyParams, workCond);
        % save local value in object properties
        gas.harmonicFrequency = harmonicFrequencyLocal;
      end
      
    end
    
    function anharmonicFrequencyLocal = evaluateAnharmonicFrequency(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.anharmonicFrequencyFunc)
        % return fixed parameter
        anharmonicFrequencyLocal = gas.anharmonicFrequency;
      else
        % call function to evaluate the value of the property
        anharmonicFrequencyLocal = gas.anharmonicFrequencyFunc(gas, gas.anharmonicFrequencyParams, workCond);
        % save local value in object properties
        gas.anharmonicFrequency = anharmonicFrequencyLocal;
      end
      
    end
    
    function rotationalConstantLocal = evaluateRotationalConstant(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.rotationalConstantFunc)
        % return fixed parameter
        rotationalConstantLocal = gas.rotationalConstant;
      else
        % call function to evaluate the value of the property
        rotationalConstantLocal = gas.rotationalConstantFunc(gas, gas.rotationalConstantParams, workCond);
        % save local value in object properties
        gas.rotationalConstant = rotationalConstantLocal;
      end
      
    end
    
    function lennardJonesDistanceLocal = evaluateLennardJonesDistance(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.lennardJonesDistanceFunc)
        % return fixed parameter
        lennardJonesDistanceLocal = gas.lennardJonesDistance;
      else
        % call function to evaluate the value of the property
        lennardJonesDistanceLocal = gas.lennardJonesDistanceFunc(gas, gas.lennardJonesDistanceParams, workCond);
        % save local value in object properties
        gas.lennardJonesDistance = lennardJonesDistanceLocal;
      end
      
    end
    
    function lennardJonesDepthLocal = evaluateLennardJonesDepth(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.lennardJonesDepthFunc)
        % return fixed parameter
        lennardJonesDepthLocal = gas.lennardJonesDepth;
      else
        % call function to evaluate the value of the property
        lennardJonesDepthLocal = gas.lennardJonesDepthFunc(gas, gas.lennardJonesDepthParams, workCond);
        % save local value in object properties
        gas.lennardJonesDepth = lennardJonesDepthLocal;
      end
      
    end
    
    function electricDipolarMomentLocal = evaluateElectricDipolarMoment(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.electricDipolarMomentFunc)
        % return fixed parameter
        electricDipolarMomentLocal = gas.electricDipolarMoment;
      else
        % call function to evaluate the value of the property
        electricDipolarMomentLocal = gas.electricDipolarMomentFunc(gas, gas.electricDipolarMomentParams, workCond);
        % save local value in object properties
        gas.electricDipolarMoment = electricDipolarMomentLocal;
      end
      
    end
    
    function electricQuadrupoleMomentLocal = evaluateElectricQuadrupoleMoment(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.electricQuadrupoleMomentFunc)
        % return fixed parameter
        electricQuadrupoleMomentLocal = gas.electricQuadrupoleMoment;
      else
        % call function to evaluate the value of the property
        electricQuadrupoleMomentLocal = ...
          gas.electricQuadrupoleMomentFunc(gas, gas.electricQuadrupoleMomentParams, workCond);
        % save local value in object properties
        gas.electricQuadrupoleMoment = electricQuadrupoleMomentLocal;
      end
      
    end
    
    function polarizabilityLocal = evaluatePolarizability(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.polarizabilityFunc)
        % return fixed parameter
        polarizabilityLocal = gas.polarizability;
      else
        % call function to evaluate the value of the property
        polarizabilityLocal = gas.polarizabilityFunc(gas, gas.polarizabilityParams, workCond);
        % save local value in object properties
        gas.polarizability = polarizabilityLocal;
      end
      
    end
    
    function fractionLocal = evaluateFraction(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.fractionFunc)
        % return fixed parameter
        fractionLocal = gas.fraction;
      else
        % call function to evaluate the value of the property
        fractionLocal = gas.fractionFunc(gas, gas.fractionParams, workCond);
        % save local value in object properties
        gas.fraction = fractionLocal;
      end
      
    end
    
    function heatCapacityLocal = evaluateHeatCapacity(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.heatCapacityFunc)
        % return fixed parameter
        heatCapacityLocal = gas.heatCapacity;
      else
        % call function to evaluate the value of the property
        heatCapacityLocal = gas.heatCapacityFunc(gas, gas.heatCapacityParams, workCond);
        % save local value in object properties
        gas.heatCapacity = heatCapacityLocal;
      end
      
    end
    
    function thermalConductivityLocal = evaluateThermalConductivity(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.thermalConductivityFunc)
        % return fixed parameter
        thermalConductivityLocal = gas.thermalConductivity;
      else
        % call function to evaluate the value of the property
        thermalConductivityLocal = gas.thermalConductivityFunc(gas, gas.thermalConductivityParams, workCond);
        % save local value in object properties
        gas.thermalConductivity = thermalConductivityLocal;
      end
      
    end

    function inFlowFractionLocal = evaluateInFlowFraction(gas, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(gas.inFlowFractionFunc)
        % return fixed parameter
        inFlowFractionLocal = gas.inFlowFraction;
      else
        % call function to evaluate the value of the property
        inFlowFractionLocal = gas.inFlowFractionFunc(gas, gas.inFlowFractionParams, workCond);
        % save local value in object properties
        gas.inFlowFraction = inFlowFractionLocal;
      end
      
    end    
        
  end

  methods (Static)

    function [gasArray, gasID] = add(gasName, gasArray)
      % addGas looks for gasName in Gas array gas. If a Gas with the same
      % name is found the function returns its ID. If no Gas is found with
      % the same name, a new Gas is added to the array and its ID is returned.

      gasID = Gas.find(gasName, gasArray);
      if gasID ~= -1
        return
      end
      switch class(gasArray)
        case 'EedfGas'
          gasArray(end+1) = EedfGas(gasName);
        case 'ChemGas'
          gasArray(end+1) = ChemGas(gasName);
        otherwise
          error('Gas class ''%s'' not recognised when trying to create gas %s\n.', class(gas), gas.name);
      end
      gasID = gasArray(end).ID;

    end

    function gasID = find(gasName, gasArray)

      for gas = gasArray
        if strcmp(gasName, gas.name)
          gasID = gas.ID;
          return
        end
      end
      gasID = -1;

    end

    function checkFractionNorm(gasArray)
    % checkFractionsNorms checks for the fractions of the different gases in gasArray to be properly normalized. In case
    % of gasArray being EedfGas it checks that the sum of the fractions of non-dummy gases (gases with non empty
    % collisionArry) add to 1. In case of gasArray being ChemGas it checks that both, volume/surface phase, gas
    % fractions add to 1. In case of non normalized gas fractions, an error message is thrown.

      volumeNorm = 0;
      surfaceNorm = 0;
      volumeNormNeedsToBeChecked = false;
      surfaceNormNeedsToBeChecked = false;
      for gas = gasArray
        if (isa(gas, 'EedfGas') && ~isempty(gas.collisionArray)) || (isa(gas, 'ChemGas') && gas.isVolumeSpecies)
          volumeNorm = volumeNorm + gas.fraction;
          volumeNormNeedsToBeChecked = true;
        elseif isa(gas, 'ChemGas') && gas.isSurfaceSpecies
          surfaceNorm = surfaceNorm + gas.fraction;
          surfaceNormNeedsToBeChecked = true;
        end
      end
      if volumeNormNeedsToBeChecked
        if abs(volumeNorm-1) > 10*eps(1)
          switch class(gasArray)
            case 'EedfGas'
              auxStr = 'Electron kinetics gas';
            case 'ChemGas'
              auxStr = 'Chemistry volume gas';
          end
          error('%s fractions are not properly normalised (Error = %e).\nPlease, check input file.', auxStr, ...
            volumeNorm-1);
        end
      end
      if surfaceNormNeedsToBeChecked
        if abs(surfaceNorm-1) > 10*eps(1)
          error(['Chemistry surface gas fractions are not properly normalised (Error = %e).\nPlease, check input ' ...
            'file.'], surfaceNorm-1);
        end
      end

    end   

  end

end
