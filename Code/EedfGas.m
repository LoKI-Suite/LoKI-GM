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

classdef EedfGas < Gas
  
  properties
    
    collisionArray = Collision.empty; % electron collisions array with all the collisions of the gas
    collisionArrayExtra = Collision.empty;
    effectivePopulations = [];        % populations to evaluate an elastic cross section from an effective one
    
    OPBParameter = [];                % 

    chemEquivalent = ChemGas.empty;
    
  end
  
  events
    
  end
  
  methods
    
    function gas = EedfGas(gasName)
      persistent lastID;
      if isempty(lastID)
        lastID = 0;
      end
      lastID = lastID + 1;
      gas.ID = lastID;
      gas.name = gasName;
      gas.stateArray = EedfState.empty;
    end
    
    function dispCollisions(gas)
      fprintf('%s collisions:\n', gas.name);
      for i = 1:length(gas.collisionArray)
        fprintf('\t %d.- %s\n', i, gas.collisionArray(i).description);
      end
    end
    
    function checkPopulationNorms(gas)
      % checkPopulationNorms checks for the population of the different states of the gas to be properly normalised,
      % i. e. the populations of all sibling states should add to one. This overloaded function avoids "dummy" states,
      % i. e. states that appear only on as product of a collision or as targets of extra collisions.
      
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
            if eleState.isTarget && eleState.population ~= 0
              gasNorm = gasNorm + eleState.population;
              % check norm of vibrational states (if they exist)
              if ~isempty(eleState.childArray)
                vibNorm = 0;
                for vibState = eleState.childArray
                  if vibState.isTarget && vibState.population ~= 0
                    vibNorm = vibNorm + vibState.population;
                    % check norm of rotational states (if they exist)
                    if ~isempty(vibState.childArray)
                      rotNorm = 0;
                      for rotState = vibState.childArray
                        if rotState.isTarget && rotState.population ~= 0
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
    
    function checkCARconditions(gas)
    % checkCARconditions checks that the particular gas fulfils the conditions to use the continuous approximation for 
    % rotation: not to have electron collisions defined that causes a rotational transition and to have defined a 
    % value for the electri quadrupole moment and the rotational constant. This function is called whenever the 
    % continuous approximation for rotations is considered for this particular gas.
      
      % check for the definition of the electric quadrupole moment
      if isempty(gas.electricQuadrupoleMoment)
        error(['A value for the electric quadrupole moment is not found for gas %s.\n' ...
          'CAR can not be taken into account without it.\n' ...
          'Please check your setup file.'], gas.name);
      end
      
      % check for the definition of the rotational constant
      if isempty(gas.rotationalConstant)
        error(['A value for the rotational constant is not found for gas %s.\n' ...
          'CAR can not be taken into account without it.\n' ...
          'Please check your setup file.'], gas.name);
      end
      
      % check for the absence of 'Rotational' collisions
      for collision = gas.collisionArray
        if strcmp(collision.type, 'Rotational')
          error(['Collision ''%s'' was found while CAR is activated for this gas.\n' ...
            'Please check your setup file'], collision.description, gas.name);
        end
      end
      
    end
    
    function collisionArray = checkElasticCollisions(gas, collisionArray)
    % checkElasticCollisions checks for each target state of the gas to have an elastic collision defined. In case the
    % Elastic collision is not found, the functions tries to obtain it from an Effective one. The function returns the
    % collision array of the gas with the Elastic collision added (if needed).
    %
    % Note: the function avoid dummy gases (gases without collisions, created for the sake of a pretty output)
      
      % avoid dummy gases
      if isempty(gas.collisionArray)
        return;
      end
      
      % loop over all the states of the gas to find an electronic state
      for state = gas.stateArray
        if strcmp(state.type, 'ele')

          % check for effective cross sections among siblings of electronic target states
          [~] = state.checkSiblingsEffectiveCrossSections();

          % check elastic collisions for electronic target states (recursively checked for children states)
          [allElasticCollisionsFound, ~] = state.checkSiblingsElasticCollisions();
          break;
        end
      end
      % create Elastic collisions in case they are not found
      if ~allElasticCollisionsFound
        rawElasticCrossSection = gas.elasticFromEffectiveCrossSection;
        % loop over all the states of the gas to find an electronic state
        for state = gas.stateArray
          if strcmp(state.type, 'ele')
            % loop over every electronic state of the gas
            for eleState = [state state.siblingArray]
              % create Elastic collision with every electronic target state
              if eleState.isTarget
                collisionArray = Collision.add('Elastic', false, eleState, eleState, 1, false, 0.0, ...
                  rawElasticCrossSection, collisionArray, 0);
              end
            end
            break;
          end
        end
      end      
    end
    
    function rawElasticCrossSection = elasticFromEffectiveCrossSection(gas)
    % elasticFrommEffectiveCrossSection evaluates an elastic cross section from an Effective one. In order to do that, 
    % it subtracts from the Effective cross section all the Excitation, Vibrational, Rotational, Ionization and 
    % Attachment cross sections of the gas weighted by the populations of the corresponding targets. 
    %
    % Note: in case that the user do not specify some particular populations for effective cross section in the setup
    % file (exceptional case), they are assumed to be those of equilibrium at 300k (regular case), which is our best 
    % guess for the conditions of measurement of an 'Effective' cross section. That is:
    %     - electronic distribution collapsed to the ground state
    %     - vibrational distribution assumed to be boltzmann at 300k
    %     - rotational distribution assumed to be boltzmann at 300k
      
      % look for Effective collision
      effectiveCollision = Collision.empty;
      for collision = gas.collisionArray
        if strcmp(collision.type, 'Effective')
          effectiveCollision = collision;
          break;
        end
      end
      if isempty(effectiveCollision)
        error(['Gas ''%s'' do not have an ''Effective'' collision defined so ''Elastic'' collisions can not be ' ...
          'evaluated from it.\nPlease, check the corresponding LXCat file.'], gas.name);
      end
      
      % initialize Elastic cross section as the Effective one
      rawElasticCrossSection = effectiveCollision.rawCrossSection;
      
      % find maximum ID of the states of the gas
      maxID = gas.stateArray(1).ID;
      for state = gas.stateArray
        if maxID < state.ID
          maxID = state.ID;
        end
      end
      
      % calculate populations for the evaluation of the elastic cross section (in case it is needed)
      if isempty(gas.effectivePopulations)
        
        % initialize vector of populations
        gas.effectivePopulations(maxID) = 0;
      
        % look for electronic ground state (target of the 'Effective' collision) and asing it a population of 1
        electronicGroundState = effectiveCollision.target;
        gas.effectivePopulations(electronicGroundState.ID) = 1;
        
        % look for vibrational states of the electronic ground state
        if ~isempty(electronicGroundState.childArray)
          % evaluate vibrational distribution (boltzmann at 300k)
          norm = 0;
          vibrationalGroundState = electronicGroundState.childArray(1);
          for state = electronicGroundState.childArray
            if isempty(state.energy)
              error(['Unable to find %s energy for the evaluation of ''Elastic'' cross section of %s.\n'...
                'Check input file'], state.name, state.gas.name);
            elseif isempty(state.statisticalWeight)
              error(['Unable to find %s statistical weight for the evaluation of ''Elastic'' cross section of %s.\n'...
                'Check input file'], state.name, state.gas.name);
            elseif state.energy < vibrationalGroundState.energy
              vibrationalGroundState = state;
            end
            gas.effectivePopulations(state.ID) = state.statisticalWeight*...
              exp(-state.energy/(Constant.boltzmannInEV*300));
            norm = norm + gas.effectivePopulations(state.ID);
          end
          for state = electronicGroundState.childArray
            gas.effectivePopulations(state.ID) = gas.effectivePopulations(state.ID)/norm;
          end
          
          % look for rotational states of the vibrational ground state of the electronic ground state
          if ~isempty(vibrationalGroundState.childArray)
            % evaluate rotational distribution (boltzmann at 300k)
            norm = 0;
            for state = vibrationalGroundState.childArray
              if isempty(state.energy)
                error(['Unable to find %s energy for the evaluation of ''Elastic'' cross section of %s.\n'...
                  'Check input file'], state.name, state.gas.name);
              elseif isempty(state.statisticalWeight)
                error(['Unable to find %s statistical weight for the evaluation of ''Elastic'' cross section of %s.'...
                  '\nCheck input file'], state.name, state.gas.name);
              end
              gas.effectivePopulations(state.ID) = state.statisticalWeight*...
                exp(-state.energy/(Constant.boltzmannInEV*300));
              norm = norm + gas.effectivePopulations(state.ID);
            end
            for state = vibrationalGroundState.childArray
              gas.effectivePopulations(state.ID) = gas.effectivePopulations(vibrationalGroundState.ID)*...
                gas.effectivePopulations(state.ID)/norm;
            end
          end
          
        end
        
      elseif length(gas.effectivePopulations)<maxID
        gas.effectivePopulations(maxID)=0;
      end
      
      % remove contributions to the effective due to the different collisional mechanisms
      for collision = gas.collisionArray
       if strcmp(collision.type, 'Effective') || strcmp(collision.type, 'Elastic')
          continue
       end
       [integralCS, momentumTransferCS] = collision.interpolatedCrossSection(rawElasticCrossSection(1,:));
       if isempty(momentumTransferCS)
         rawElasticCrossSection(2,:) = rawElasticCrossSection(2,:) - gas.effectivePopulations(collision.target.ID)*...
           integralCS;
       else
         rawElasticCrossSection(2,:) = rawElasticCrossSection(2,:) - gas.effectivePopulations(collision.target.ID)*...
           momentumTransferCS;
       end
       % remove contributions to the effective due to the super-elastic mechanism (in case it is defined)
       if collision.isReverse
         [integralCS, momentumTransferCS] = collision.superElasticCrossSection(rawElasticCrossSection(1,:));
         if isempty(momentumTransferCS)
           rawElasticCrossSection(2,:) = rawElasticCrossSection(2,:) - ...
             gas.effectivePopulations(collision.productArray.ID)*integralCS;
         else
           rawElasticCrossSection(2,:) = rawElasticCrossSection(2,:) - ...
             gas.effectivePopulations(collision.productArray.ID)*momentumTransferCS;
         end
       end
      end
      
      % check for negative values in the Elastic cross section
      if any(rawElasticCrossSection(2,:)<0)
        rawElasticCrossSection(2,rawElasticCrossSection(2,:)<0) = 0;
        warning(['Negative values obtained when evaluating an Elastic cross section from an Effective one (%s).\n'...
          'Negative values have been clipped to 0 and unreliable results may be obtained.\n'...
          'Please, carefully check inputs and outputs of your simulation'], gas.name);
      end
      
    end
    
  end
  
  methods (Static)
    
  end
  
end