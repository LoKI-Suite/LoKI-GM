function [rateCoeff, dependent] = binaryTransportChantry(~, ~, totalGasDensity, reaction, rateCoeffParams, chemistry)
% binaryTransportChantry evaluates the diffusion rate coefficient of a particular neutral species interacting with the wall, 
%  adopting a binary diffusion model where the species diffuse in the gas (neglecting the effects of multicomponent transport),
%  taking into account both the binary diffusion time and the wall reaction time assuming a flux reaching the wall 
%  (with a certain wall recombination/deactivation probability). 
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003

% NOTE: this implementation of the binary transport uses the "heuristic" formula provided by Chantry 
%  for the effective diffusion length.
% Chantry P J 1987 Journal of Applied Physics, 62(4) 1141-1148. https://doi.org/10.1063/1.339662

% For more info check the documentation.

  persistent charDiffLengthSquared;
  persistent gammaSumArray;
  persistent dependentInfo;

  % local save of the ID of the reactant species
  reactantID = reaction.reactantArray(1).ID;

  % --- performance sensitive calculations ---

  % calculations performed once per simulation
  if isempty(dependentInfo)
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % evaluate geometrical parameters for a cylinder (including infinitely long cylinder and slab limiting cases)
    if chemistry.workCond.chamberLength == 0        % infinitely long cylinder
      charDiffLengthSquared = (chemistry.workCond.chamberRadius/firstBesselZero)^2;
    elseif chemistry.workCond.chamberRadius == 0    % infinitely wide cylinder (slab)
      charDiffLengthSquared = (chemistry.workCond.chamberLength/pi)^2;
    else                                  % finite cylinder
      charDiffLengthSquared = 1/((firstBesselZero/chemistry.workCond.chamberRadius)^2+(pi/chemistry.workCond.chamberLength)^2);
    end
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % initialize array of gammaSum (total loss wall probability)
    gammaSumArray = zeros(size(chemistry.stateArray));
  end

  % evaluate total wall reaction probability for the reactant specie (only done once per lost species)
  if gammaSumArray(reactantID) == 0
    % evaluate total wall reaction coefficient
    for rxn = chemistry.stateArray(reactantID).reactionsDestruction
      if rxn.isTransport
        if strcmp(rxn.type, reaction.type)
          if ~isempty(rxn.rateCoeffParams) && ...
              isnumeric(rxn.rateCoeffParams{1}) && rxn.rateCoeffParams{1}<=1 && rxn.rateCoeffParams{1}>0
            gammaSumArray(reactantID) = gammaSumArray(reactantID) + rxn.rateCoeffParams{1};
          else % error checking
            error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''wallCoefficient'' ' ...
              'not probided (or wrong value) in the corresponding ''.chem'' file.\nPlease, fix the problem and ' ...
              'run the code again'], rxn.type, rxn.description);
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
    if gammaSumArray(reactantID) == 0
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability must be different than zero.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    elseif gammaSumArray(reactantID) > 1
      error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
        'Total wall reaction probability can not be larger than 1.\nPlease, fix the problem and run the code ' ...
        'again'], reaction.type, reaction.description);
    end
  end

  % --- regular calculations ---

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
    thermalVelocity = sqrt(8*Constant.boltzmann*chemistry.workCond.(temperatureStr)/(pi*chemistry.stateArray(reactantID).gas.mass));
  end

  % evaluate diffusion rate coefficient
  redDiffCoeff = reaction.reactantArray(1).evaluateReducedDiffCoeff(chemistry.workCond);

  % evaluate wall reaction coefficient
  gamma = rateCoeffParams{1};
  % evaluate rate coefficient
  rateCoeff = (gamma/gammaSumArray(reactantID)) / ( charDiffLengthSquared*totalGasDensity/redDiffCoeff + ...
    chemistry.workCond.volumeOverArea*4*(1-gammaSumArray(reactantID)/2)/(thermalVelocity*gammaSumArray(reactantID)) );

  % set function dependencies
  dependent = dependentInfo;
  if ~dependent.onGasTemperature && ~isempty(reaction.reactantArray(1).reducedDiffCoeffFunc)
    for param = reaction.reactantArray(1).reducedDiffCoeffParams
      if strcmp(param{1}, 'gasTemperature')
        dependent.onGasTemperature = true;
      end
    end
  end

end
