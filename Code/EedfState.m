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

classdef EedfState < State
  
  properties
    
    isTarget = false;                 % true if the state is the target of a collision (not extra), false otherwise
    collisionArray = Collision.empty; % handle array to the collisions of which the target is state
    collisionArrayExtra = Collision.empty;

    chemEquivalent = ChemState.empty;

  end
  
  events
    
  end
  
  methods
    
    function state = EedfState(gas, ionCharg, eleLevel, vibLevel, rotLevel)
      persistent lastID;
      if isempty(lastID)
        lastID = 0;
      end
      lastID = lastID + 1;
      state.ID = lastID;
      state.gas = gas;
      state.ionCharg = ionCharg;
      state.eleLevel = eleLevel;
      state.vibLevel = vibLevel;
      state.rotLevel = rotLevel;
      if isempty(ionCharg)
        if isempty(rotLevel)
          if isempty(vibLevel)
            state.type = 'ele';
          else
            state.type = 'vib';
          end
        else
          state.type = 'rot';
        end
      else
        state.type = 'ion';
      end
      state.parent = EedfState.empty;
      state.siblingArray = EedfState.empty;
      state.childArray = EedfState.empty;
      state.addFamily;
      gas.stateArray(end+1) = state;
      state.evaluateName;
    end

    function [doubleEffectiveCrossSectionFound] = checkSiblingsEffectiveCrossSections(state)
      % Check if there are two target siblings of the state (including itself) with effective cross sections defined
      firstSiblingWithEffectiveCrossSection = false;
      doubleEffectiveCrossSectionFound = false;
      for sibling = [state, state.siblingArray]
        % Avoid dummy states
        if ~sibling.isTarget
          continue;
        end
        % Check if current sibling has an effective cross section defined
        for collision = sibling.collisionArray
          if strcmp(collision.type, 'Effective')
            if ~firstSiblingWithEffectiveCrossSection
                firstSiblingWithEffectiveCrossSection = true;
            elseif ~doubleEffectiveCrossSectionFound
                doubleEffectiveCrossSectionFound = true;
                error(['Effective cross section duplicated among %s and its siblings.\n' ...
                    'Please, check the corresponding LXCat file(s).'], sibling.name);
            end    
          end
        end
      end 
    end 

    function [allElasticCollisionsFound, anyElasticCollisionFound] = checkSiblingsElasticCollisions(state)
      % Check if all target siblings (including itself) of the state have elastic collisions defined
      anyElasticCollisionFound = false;
      allElasticCollisionsFound = true;
      for sibling = [state, state.siblingArray]
        % Avoid dummy states
        if ~sibling.isTarget
          continue;
        end
        % Check if current sibling has an elastic collision defined explicitly
        currentSiblingDoesNotHaveElasticCollision = true;
        for collision = sibling.collisionArray
          if strcmp(collision.type, 'Elastic')
            currentSiblingDoesNotHaveElasticCollision = false;
            anyElasticCollisionFound = true;
            break;
          end
        end
        % Check elastic collisions among its chldren (if any)
        if ~isempty(sibling.childArray)
          [allElasticCollisionFoundAmongChildren, anyElasticCollisionFoundAmongChildren] = ...
            sibling.childArray(1).checkSiblingsElasticCollisions();
          % check if information is provided/duplicated with elastic collisions among its children
          if currentSiblingDoesNotHaveElasticCollision
            currentSiblingDoesNotHaveElasticCollision = ~allElasticCollisionFoundAmongChildren;
          else
            if anyElasticCollisionFoundAmongChildren
              error(['Elastic collision information is duplicated among %s and its children.\n'...
                'Please, check the corresponding LXCat file(s).'], sibling.name);
            end
          end
        end
        % If current sibling finally does not have elastic collision, set global flag to false and return
        if currentSiblingDoesNotHaveElasticCollision
          allElasticCollisionsFound = false;
        end
      end
      % If any, but not all, elastic collisions are found among siblings, raise error
      if anyElasticCollisionFound && ~allElasticCollisionsFound
        error(['Some elastic collision information is missing among %s siblings.\n'...
          'Please, check the corresponding LXCat file(s).'], state.name);
      end
    end
    
  end
  
end