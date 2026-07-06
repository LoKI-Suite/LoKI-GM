function [rateCoeff, dependent] = wallReaction(~, ~, ~, reaction, rateCoeffParams, chemistry)
% wallReaction evaluates the transport rate coefficient of a particular species interacting with the wall, 
%  assuming a flux reaching the wall (with a certain wall recombination/deactivation probability), 
%  and neglecting the diffusion transport 
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Chantry P J 1987 Journal of Applied Physics, 62(4) 1141-1148. https://doi.org/10.1063/1.339662

% For more info check documentation.
  
  persistent gammaSumArray
  persistent dependentInfo;

  % local save of the ID of the reactant species
  reactantID = reaction.reactantArray(1).ID;
  
  % --- performance sensitive calculations ---
  
  % calculations performed once per simulation 
  if isempty(dependentInfo)
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % initialize array of gammaSum (total loss wall probability) and thermal velocities
    gammaSumArray = zeros(size(chemistry.stateArray));
  end
  % evaluate total wall reaction probability (total loss) and thermal velocities
  if gammaSumArray(reactantID) == 0
    % evaluate total wall reaction coefficient
    for reactionDestruction = chemistry.stateArray(reactantID).reactionsDestruction
      if reactionDestruction.isTransport
        if strcmp(reactionDestruction.type, reaction.type)
          if ~isempty(reactionDestruction.rateCoeffParams) && isnumeric(reactionDestruction.rateCoeffParams{1}) && ...
            reactionDestruction.rateCoeffParams{1}<=1 && reactionDestruction.rateCoeffParams{1}>0
            gammaSumArray(reactantID) = gammaSumArray(reactantID) + reactionDestruction.rateCoeffParams{1};
          else % error checking
            error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''wallCoefficient'' ' ...
              'not probided (or wrong value) in the corresponding ''.chem'' file.\nPlease, fix the problem and ' ...
              'run the code again'], reaction.type, reaction.description);
          end
        else % error checking
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            'It has been found another transport reaction of a different type for the same species:\n%s.\nPlease, ' ...
            'fix the problem and run the code again'], reaction.type, reaction.description, ...
            reactionDestruction.description);
        end
      end
    end
    % error checking for limiting values of gammaSumArray
    if gammaSumArray(reactantID) == 0
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability must be different than zero.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    elseif gammaSumArray(reactantID) > 1
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability can not be larger than 1.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    end
    % error checking for species mass
    if isempty(reaction.reactantArray(1).gas.mass)
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not ' ...
        'found.\nPlease, fix the problem and run the code again'], reaction.type, reaction.description, ...
        reaction.reactantArray(1).name);
    end
  end
    
  % --- regular calculations ---
  
  % local copy of working conditions
  workCond = chemistry.workCond;
  % evaluate wall reaction coefficient
  gamma = rateCoeffParams{1};

  % evaluate thermal velocity
  temperatureStr = rateCoeffParams{2};
  if ~any(strcmp(temperatureStr, {'gasTemperature' 'nearWallTemperature' 'wallTemperature'}))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' is not a valid temperature ' ...
      'to evaluate the thermal velocity with.\nChoose one of the following: ''gasTemperature'', ''nearWallTemperature'' ' ...
      'or ''wallTemperature''.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{2});
  elseif isempty(workCond.(temperatureStr))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' not found in the working ' ...
      'conditions object.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{2});
  else
    thermalVelocity = sqrt(8*Constant.boltzmann*workCond.(temperatureStr)/(pi*chemistry.stateArray(reactantID).gas.mass));
  end

  % evaluate rate coefficient
  rateCoeff = workCond.areaOverVolume*thermalVelocity/4*gamma/(1-gammaSumArray(reactantID)/2);
  
  % set function dependencies
  dependent = dependentInfo;

end