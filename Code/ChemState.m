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

classdef ChemState < State
  
  properties
    
    reactionsCreation = Reaction.empty;     % handle to the reactions where this state is created (forward reaction)
    reactionsDestruction = Reaction.empty;  % handle to the reactions where this state is destroyed (forward reaction)
    
    eedfEquivalent = EedfState.empty;       % handle to the corresponding equivalent eedf state (in case it exists)
    
    inFlowRelDensity;                       % handle to the inflow population of a state relative to its siblings

  end
  
  events
    
  end
  
  methods (Access = public)
    
    function state = ChemState(gas, ionCharg, eleLevel, vibLevel, rotLevel)
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
      state.parent = ChemState.empty;
      state.siblingArray = ChemState.empty;
      state.childArray = ChemState.empty;
      state.addFamily;
      gas.stateArray(end+1) = state;
      state.evaluateName;
    end
    
    function isVolumeSpecies = isVolumeSpecies(state)
    % isVolumeSpecies return a boolean value that is true if the parent gas of the state is defined as a volume species
    % and false otherwise
      
      isVolumeSpecies = state.gas.isVolumeSpecies;

    end

    function isSurfaceSpecies = isSurfaceSpecies(state)
    % isSurfaceSpecies return a boolean value that is true if the parent gas of the state is defined as a surface species
    % and false otherwise
      
      isSurfaceSpecies = state.gas.isSurfaceSpecies;

    end

    function linkWithElectronKineticsState(chemState, eedfState)
      
      % save handle to the equivalent electron kinetic state
      chemState.eedfEquivalent = eedfState;
      eedfState.chemEquivalent = chemState;
      
      % copy properties of the equivalent electron kinetic state (avoiding 'ID', 'gas', 'parent', 'siblingArray' and 
      % 'childArray'  properties)
      eedfStateProperties = fields(eedfState)';
      chemStateProperties = fields(chemState)';
      for property = eedfStateProperties
        if ~strcmp(property{1},'ID') && ~strcmp(property{1}, 'gas') && ~strcmp(property{1}, 'parent') && ...
            ~strcmp(property{1}, 'siblingArray') && ~strcmp(property{1}, 'childArray') && ...
            any(strcmp(property, chemStateProperties))
          chemState.(property{1}) = eedfState.(property{1});
        end
      end
      
    end
    
    function evaluateInFlowRelDensity(state)
    % evaluateInFlowRelDensity evaluates the inflow "relative density" of a certain state.
    % For a rotational state this means its inFlowPopulation, 
    % multiplied by the inFlowPopulation of its vibrational parent, 
    % multiplied by the inFlowPopulation of the corresponding electronic parent, 
    % multiplied by the gas inFlowFraction.
    % For vibrational or electronic states the evaluation is analogue. 
    % 
    % NOTE! this value is relative
    
      switch state.type
        case 'rot'
        state.inFlowRelDensity = state.inFlowPopulation*state.parent.inFlowPopulation*...
          state.parent.parent.inFlowPopulation*state.gas.inFlowFraction;
        case 'vib'
        state.inFlowRelDensity = state.inFlowPopulation*state.parent.inFlowPopulation*...
          state.gas.inFlowFraction;
        case 'ele'
        state.inFlowRelDensity = state.inFlowPopulation*state.gas.inFlowFraction;
        case 'ion'
        state.inFlowRelDensity = state.inFlowPopulation*state.gas.inFlowFraction;
      end
        
    end    
    
  end
  
end
