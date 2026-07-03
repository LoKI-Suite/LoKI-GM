function [rateCoeff, dependent] = wallReactionSurface(~, ~, ~, reaction, rateCoeffParams, chemistry)
% wallReactionSurface evaluates the transport rate coefficient of a particular species interacting with the wall, 
%  assuming a flux reaching the wall (with a certain wall recombination/deactivation probability), 
%  and neglecting the diffusion transport 
% The rate coefficient is modified to consider the collision of the volume species with ONE surface species 
%  (dividing by surfaceSiteDensity*AreaOverVolume).
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Chantry P J 1987 Journal of Applied Physics, 62(4) 1141-1148. https://doi.org/10.1063/1.339662
%
% For more info check documentation.

  persistent gammaSumArray
  persistent dependentInfo;

  % local save of the ID of the volume reactant species
  for i = 1:length(reaction.reactantArray)
    if reaction.reactantArray(i).isVolumeSpecies
      volumeReactantID = reaction.reactantArray(i).ID;
    end
  end

  % --- performance sensitive calculations ---

  % calculations performed once per simulation
  if isempty(dependentInfo)
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % initialize array of gammaSum (total loss wall probability) and thermal velocities
    gammaSumArray = zeros(size(chemistry.stateArray));
  end
  % evaluate total wall reaction probability (total loss) and thermal velocities
  if gammaSumArray(volumeReactantID) == 0
    % evaluate total wall reaction coefficient
    for rxn = chemistry.stateArray(volumeReactantID).reactionsDestruction
      if rxn.isTransport
        if strcmp(rxn.type, reaction.type)
          if ~isempty(rxn.rateCoeffParams) && ...
              isnumeric(rxn.rateCoeffParams{1}) && rxn.rateCoeffParams{1}<=1 && rxn.rateCoeffParams{1}>0
            gammaSumArray(volumeReactantID) = gammaSumArray(volumeReactantID) + rxn.rateCoeffParams{1};
          else % error checking
            error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''wallCoefficient'' ' ...
              'not probided (or wrong value) in the corresponding ''.chem'' file.\nPlease, fix the problem and ' ...
              'run the code again'], reaction.type, rxn.description);
          end
        else % error checking
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            'It has been found another transport reaction of a different type for the same species:\n%s.\nPlease, ' ...
            'fix the problem and run the code again'], reaction.type, ...
            reaction.description, rxn.description);
        end
      end
    end
    % error checking for limiting values of gammaSumArray
    if gammaSumArray(volumeReactantID) == 0
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability must be different than zero.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    elseif gammaSumArray(volumeReactantID) > 1
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability can not be larger than 1.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    end
    % error checking for species mass
    if isempty(reaction.reactantArray(1).gas.mass)
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not ' ...
        'found.\nPlease, fix the problem and run the code again'], reaction.type, ...
        reaction.description, reaction.reactantArray(1).name);
    end
  end

  % --- regular calculations ---

  % evaluate wall reaction coefficient
  gamma = rateCoeffParams{1};

  % evaluate thermal velocity
  temperatureStr = rateCoeffParams{2};
  if ~any(strcmp(temperatureStr, {'gasTemperature' 'nearWallTemperature' 'wallTemperature'}))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' is not a valid temperature ' ...
      'to evaluate the thermal velocity with.\nChoose one of the following: ''gasTemperature'', ''nearWallTemperature'' ' ...
      'or ''wallTemperature''.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{2});
  elseif isempty(chemistry.workCond.(temperatureStr))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' not found in the working ' ...
      'conditions object.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{2});
  else
    thermalVelocity = sqrt(8*Constant.boltzmann*chemistry.workCond.(temperatureStr)/(pi*chemistry.stateArray(volumeReactantID).gas.mass));
  end

  % evaluate rate coefficient
  rateCoeff = chemistry.workCond.areaOverVolume*thermalVelocity/4*gamma/(1-gammaSumArray(volumeReactantID)/2);
  rateCoeff = rateCoeff/(chemistry.workCond.surfaceSiteDensity*chemistry.workCond.areaOverVolume);

  % set function dependencies
  dependent = dependentInfo;

end
