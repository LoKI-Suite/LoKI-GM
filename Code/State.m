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

classdef State < handle
  %State Class that stores the information of a certain state of a gas
  %   Class that stores the information of a certain state of a gas. The
  %   information stored here (in particular the population of the state)
  %   is to be used in a Boltzmann solver to obtain the EEDF.
  
  properties

    ID = -1;                        % ID that identifies the state in the state array of a gas
                                    
    type = '';                      % type of state, options are: 'ele', 'vib', 'rot' or 'ion'
    ionCharg = [];                  % ionic level of the state
    eleLevel = [];                  % electronic level of the state
    vibLevel = [];                  % vibrational level of the state
    rotLevel = [];                  % rotational level of the state
    name = [];                      % name of the state gas(ionCharg,eleLevel,vibLevel,rotLevel)
                                    
    gas;                            % handle to the gas that the state belongs to (initialized in the subclass)
    parent;                         % handle to the parent state (e.g. the parent of N2(X,v=0) is -> N2(X)), electronic 
                                    %  states have no parent) (initicialized in the subclass)
    siblingArray;                   % handle array to the states with the same parent (e.g. siblings of N2(X) are N2(A),
                                    %  N2(B), ...) (initicialized in the subclass)
    childArray;                     % handle array to the child states (e.g. children of N2(X) are N2(X,v=0), N2(X,v=1), 
                                    %  N2(X,v=2),...) (initicialized in the subclass)
                                    
    energy = [];                    % energy of the state
    energyFunc = [];                % handle to function that evaluates the energy of the state
    energyParams = {};              % cell array of parameters needed by energyFunc
                                    
    statisticalWeight = [];         % statistical weight of the state
    statisticalWeightFunc = [];     % handle to function that evaluates the statistical weight of the state
    statisticalWeightParams = {};   % cell array of parameters needed by statisticalWeightFunc
                                    
    population = 0;                 % population of the state relative to its siblings
    populationFunc = [];            % handle to function that evaluates the population 
    populationParams = {};          % cell array of parameters needed by populationFunc
                                    
    density = [];                   % absolute density (relative to the gas density)
                                    
    reducedDiffCoeff = [];          % reduced free diffusion coefficient
    reducedDiffCoeffFunc = [];      % handle to function that evaluates the reduced free diffusion coefficient
    reducedDiffCoeffParams = {};    % cell array of parameters needed by reducedDiffCoeffFunc
                                    
    reducedMobility = [];           % reduced free mobility
    reducedMobilityFunc = [];       % handle to function that evaluates the reduced free mobility
    reducedMobilityParams = {};     % cell array of parameters needed by reducedMobillityFunc
    
    flowBarrier = [];               % flow (factor) barrier 
    
    inFlowPopulation = 0;           % inflow population of the state relative to its siblings
    inFlowPopulationFunc = [];      % handle to function that evaluates the inflow population 
    inFlowPopulationParams = {};    % cell array of parameters needed by inFlowPopulationFunc    
    
    gasTemperatureListener = [];    % handle to the listener of changes in gas temperature (working conditions property)
    
  end
  
  events
    
  end

  methods (Access = public)
    
    function addFamily(state)
      % addFamily find all the relatives of state and stores the
      % corresponding information in their properties parent, sibling and
      % child.
      
      switch state.type
        
        case 'rot'
          for auxState = state.gas.stateArray
            if ( strcmp(auxState.type, 'vib') && ...
                 strcmp(auxState.eleLevel, state.eleLevel) && ...
                 strcmp(auxState.vibLevel, state.vibLevel) )
              state.parent = auxState;
              auxState.childArray(end+1) = state;
            elseif ( strcmp(auxState.type, 'rot') && ...
                     strcmp(auxState.eleLevel, state.eleLevel) && ...
                     strcmp(auxState.vibLevel, state.vibLevel) )
              state.siblingArray(end+1) = auxState;
              auxState.siblingArray(end+1) = state;
            end
          end
          
        case 'vib'
          for auxState = state.gas.stateArray
            if ( strcmp(auxState.type, 'ele') && ...
                 strcmp(auxState.eleLevel, state.eleLevel) )
              state.parent = auxState;
              auxState.childArray(end+1) = state;
            elseif ( strcmp(auxState.type, 'vib') && ...
                     strcmp(auxState.eleLevel, state.eleLevel) )
              state.siblingArray(end+1) = auxState;
              auxState.siblingArray(end+1) = state;
            elseif ( strcmp(auxState.type, 'rot') && ...
                     strcmp(auxState.eleLevel, state.eleLevel) && ...
                     strcmp(auxState.vibLevel, state.vibLevel) )
              state.childArray(end+1) = auxState;
              auxState.parent = state;
            end
          end
          
        case 'ele'
          for auxState = state.gas.stateArray
            if strcmp(auxState.type, 'ele')
              state.siblingArray(end+1) = auxState;
              auxState.siblingArray(end+1) = state;
            elseif ( strcmp(auxState.type, 'vib') && ...
                     strcmp(auxState.eleLevel, state.eleLevel) )
              state.childArray(end+1) = auxState;
              auxState.parent = state;
            end
          end
          
        case 'ion'
          for auxState = state.gas.stateArray
            if strcmp(auxState.type, 'ion')
              state.siblingArray(end+1) = auxState;
              auxState.siblingArray(end+1) = state;
            end
          end
      end
      
    end
    
    function evaluateDensity(state)
    % evaluateDensity evaluates the absolute density of a certain state.
    % For a rotational state this means its population, multiplied by the
    % population of its vibrational parent, multiplied by the population of
    % the corresponding electronic parent, multiplied by the gas fraction.
    % For vibrational or electronic states the evaluation is analogue. 
    % 
    % NOTE! this density is still relative to the total gas density.
    
      switch state.type
        case 'rot'
        state.density = state.population*state.parent.population*...
          state.parent.parent.population*state.gas.fraction;
        case 'vib'
        state.density = state.population*state.parent.population*...
          state.gas.fraction;
        case 'ele'
        state.density = state.population*state.gas.fraction;
        case 'ion'
        state.density = state.population*state.gas.fraction;
      end
        
    end
    
    function disp(state)
      fprintf('ID: %d\n', state.ID);
      fprintf('name: %s\n', state.name);
      if ~isempty(state.energy)
        fprintf('energy: %g\n', state.energy);
      end
      if ~isempty(state.statisticalWeight)
        fprintf('statisticalWeight: %g\n', state.statisticalWeight);
      end
      if ~isempty(state.population)
        fprintf('population: %g\n', state.population);
      end
      if ~isempty(state.density)
        fprintf('density: %g\n', state.density);
      end
      if ~isempty(state.parent)
        fprintf('parent: %s\n', state.parent.name);
      end
      if ~isempty(state.siblingArray)
        fprintf('siblings:\n');
        for sibling = state.siblingArray
          fprintf('\t %s\n', sibling.name);
        end
      end
      if ~isempty(state.childArray)
        fprintf('childrends:\n')
        for child = state.childArray
          fprintf('\t %s\n', child.name);
        end
      end
    end
  
    function stateName = evaluateName(state)
      stateName = [state.gas.name '('];
      if ~isempty(state.ionCharg)
        stateName = [stateName state.ionCharg ','];
      end
      stateName = [stateName state.eleLevel];
      if ~isempty(state.vibLevel)
        stateName = [stateName ',v=' state.vibLevel];
        if ~isempty(state.rotLevel)
          stateName = [stateName ',J=' state.rotLevel];
        end
      end
      stateName = [stateName ')'];
      state.name = stateName;
    end
    
    function gasMass = mass(state)
      % mass is an alias function to obtain the mass of a certain state which is the mass of the parent gas
      
      gasMass = state.gas.mass;
      
    end
    
    function energyLocal = evaluateEnergy(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.energyFunc)
        % return fixed parameter
        energyLocal = state.energy;
      else
        % call function to evaluate the value of the property
        energyLocal = state.energyFunc(state, state.energyParams, workCond);
        % save local value in object properties
        state.energy = energyLocal;
      end
      
    end
    
    function statisticalWeightLocal = evaluateStatisticalWeight(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.statisticalWeightFunc)
        % return fixed parameter
        statisticalWeightLocal = state.statisticalWeight;
      else
        % call function to evaluate the value of the property
        statisticalWeightLocal = state.statisticalWeightFunc(state, state.statisticalWeightParams, workCond);
        % save local value in object properties
        state.statisticalWeight = statisticalWeightLocal;
      end
      
    end
    
    function populationLocal = evaluatePopulation(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.populationFunc)
        % return fixed parameter
        populationLocal = state.population;
      else
        % call function to evaluate the value of the property
        populationLocal = state.populationFunc(state, state.populationParams, workCond);
        % save local value in object properties
        state.population = populationLocal;
      end
      
    end
    
    function reducedDiffCoeffLocal = evaluateReducedDiffCoeff(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.reducedDiffCoeffFunc)
        % return fixed parameter
        reducedDiffCoeffLocal = state.reducedDiffCoeff;
      else
        % call function to evaluate the value of the property
        reducedDiffCoeffLocal = state.reducedDiffCoeffFunc(state, state.reducedDiffCoeffParams, workCond);
        % save local value in object properties
        state.reducedDiffCoeff = reducedDiffCoeffLocal;
      end
      
    end
    
    function reducedMobilityLocal = evaluateReducedMobility(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.reducedMobilityFunc)
        % return fixed parameter
        reducedMobilityLocal = state.reducedMobility;
      else
        % call function to evaluate the value of the property
        reducedMobilityLocal = state.reducedMobilityFunc(state, state.reducedMobilityParams, workCond);
        % save local value in object properties
        state.reducedMobility = reducedMobilityLocal;
      end
      
    end

    function inFlowPopulationLocal = evaluateInFlowPopulation(state, workCond)
      
      % checks if the property is to be evaluated with a function or a fixed parameter
      if isempty(state.inFlowPopulationFunc)
        % return fixed parameter
        inFlowPopulationLocal = state.inFlowPopulation;
      else
        % call function to evaluate the value of the property
        inFlowPopulationLocal = state.inFlowPopulationFunc(state, state.inFlowPopulationParams, workCond);
        % save local value in object properties
        state.inFlowPopulation = inFlowPopulationLocal;
      end
      
    end    
        
  end
  
  methods (Static)
    
    function [stateArray, stateID] = add(gas, ionCharg, eleLevel, ...
      vibLevel, rotLevel, stateArray)
      % addState analyse a state string read from an LXCat file and finds if
      % the state is already defined in the corresponding gas. If the state
      % (or even the gas) is not defined the function creates it. The
      % function returns the Gas array gas and the ID of the state.
      
      stateID = State.find(gas.name, ionCharg, eleLevel, vibLevel, ...
        rotLevel, stateArray);
      if stateID ~= -1
        return;
      end
      switch class(gas)
        case 'EedfGas'
          stateArray(end+1) = EedfState(gas, ionCharg, eleLevel, vibLevel, rotLevel);
        case 'ChemGas'
          stateArray(end+1) = ChemState(gas, ionCharg, eleLevel, vibLevel, rotLevel);
        otherwise
          error('Gas class ''%s'' not recognised when trying to create an state of %s\n.', class(gas), gas.name);
      end
      stateID = stateArray(end).ID;
    
    end
    
    function stateID = find(gasName, ionCharg, eleLevel, vibLevel, ...
      rotLevel, stateArray)
      
      stateID = [];
      if strcmp(eleLevel, Parse.wildCardChar)
        for i = 1:length(stateArray)
          state = stateArray(i);
          if ( strcmp(gasName, state.gas.name) && ...
               strcmp('ele', state.type) )
            stateID = [state.ID state.siblingArray.ID];
            break;
          end
        end
      elseif strcmp(vibLevel, Parse.wildCardChar)
        for i = 1:length(stateArray)
          state = stateArray(i);
          if ( strcmp(gasName, state.gas.name) && ...
               strcmp(eleLevel, state.eleLevel) && ...
               strcmp('ele', state.type) )
            stateID = [state.childArray.ID];
            break;
          end
        end
      elseif strcmp(rotLevel, Parse.wildCardChar)
        for i = 1:length(stateArray)
          state = stateArray(i);
          if ( strcmp(gasName, state.gas.name) && ...
               strcmp(eleLevel, state.eleLevel) && ...
               strcmp(vibLevel, state.vibLevel) && ...
               strcmp('vib', state.type) )
            stateID = [state.childArray.ID];
            break;
          end
        end
      else
        for i = 1:length(stateArray)
          state = stateArray(i);
          if ( strcmp(gasName, state.gas.name) && ...
               strcmp(ionCharg, state.ionCharg) && ...
               strcmp(eleLevel, state.eleLevel) && ...
               strcmp(vibLevel, state.vibLevel) && ...
               strcmp(rotLevel, state.rotLevel) )
            stateID = state.ID;
            break;
          end
        end
      end
      if isempty(stateID)
        stateID = -1;
      end
      
    end
    
    function stateArray = fixOrphanStates(stateArray)
      % fixOrphanStates looks for orphan states and create the needed
      % parent states
      
      for state = stateArray
        if strcmp(state.type, 'rot') && isempty(state.parent)
          stateArray = State.add(state.gas, state.ionCharg, ...
            state.eleLevel, state.vibLevel, '', stateArray);
        end
      end
      for state = stateArray
        if strcmp(state.type, 'vib') && isempty(state.parent)
          stateArray = State.add(state.gas, state.ionCharg, ...
            state.eleLevel, '', '', stateArray);
        end
      end
      
    end
    
    function equalStates = compareArrays(stateArray1, stateArray2)
      
      equalStates = false;
      numberOfStates = length(stateArray1);
      if numberOfStates == length(stateArray2)
        if numberOfStates == 0
          equalStates = true;
        elseif numberOfStates == 1
          if stateArray1.ID == stateArray2.ID
            equalStates = true;
          end
        else
          for i = 1:numberOfStates
            for j = 1:numberOfStates
              if stateArray1(i).ID == stateArray2(j).ID
                break;
              elseif j == length(stateArray2)
                return;
              end
            end
          end
          equalStates = true;
        end
      end
    end
    
  end

end

