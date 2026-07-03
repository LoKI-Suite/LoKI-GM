function [rateCoeff, dependent] = binaryDiff(time, ~, totalGasDensity, reaction, rateCoeffParams, chemistry)
% binaryDiff evaluates the diffusion rate coefficient of a particular neutral species interacting with the wall, 
%  adopting a binary diffusion model where the species diffuse in the gas (neglecting the effects of multicomponent transport),
%  assuming a density profile that is zero at the walls with wall recombination at a certain probability
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003

% For more info check documentation.
  
  persistent firstBesselZero;
  persistent initialTime;
  persistent characteristicLengthSquared;
  persistent dependentInfo;
  
  % --- performance sensitive calculations ---
  
  % initialize variables the first time the binaryDiff function is called
  if isempty(firstBesselZero)
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    initialTime = time;
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', false, 'onElectronKinetics', false);
  end
  
  % --- time independent calculations (also performance sensitive) ---
  
  % evaluate the squared characteristic length (every initial time, performance sensitive)
  workCond = chemistry.workCond;
  if time == initialTime
    if workCond.chamberLength == 0
      characteristicLengthSquared = (workCond.chamberRadius/firstBesselZero)^2;
    elseif workCond.chamberRadius == 0
      characteristicLengthSquared = (workCond.chamberLength/pi)^2;
    else
      characteristicLengthSquared = 1/((workCond.chamberRadius/firstBesselZero)^-2+(workCond.chamberLength/pi)^-2);
    end
  end
  
  % --- regular calculations ---
  
  % evaluate diffusion rate coefficient
  wallProbability = rateCoeffParams{1};
  redDiffCoeff = reaction.reactantArray.evaluateReducedDiffCoeff(workCond);
  
  if isempty(redDiffCoeff)
    error(['Error found when evaluating the rate coefficient of the reaction:\n%s\n'...
      '''reducedDiffCoeff'' property of the state ''%s'' not found.\n'...
      'Please check your setup file'], reaction.descriptionExtended, reaction.reactantArray.name);
  end
  rateCoeff = wallProbability*redDiffCoeff/(characteristicLengthSquared*totalGasDensity);
  
  % set function dependencies
  dependent = dependentInfo;
  if ~dependent.onGasTemperature && ~isempty(reaction.reactantArray.reducedDiffCoeffFunc)
    for param = reaction.reactantArray.reducedDiffCoeffParams
      if strcmp(param{1}, 'gasTemperature')
        dependent.onGasTemperature = true;
      end
    end
  end   

end