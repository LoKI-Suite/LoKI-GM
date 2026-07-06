function [rateCoeff, dependent] = multicomponentTransportSurface(~, densitiesAll, ~, reaction, rateCoeffParams, chemistry)
% multicomponentTransportSurface evaluates the transport rate coefficient of a particular species interacting with the wall, 
%  taking into account both the binary diffusion time and the wall reaction time assuming a flux reaching the wall 
%  (with a certain wall recombination/deactivation probability). 
% The rate coefficient is modified to consider the collision of the volume species with ONE surface species 
%  (dividing by surfaceSiteDensity*AreaOverVolume).
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Chantry P J 1987 Journal of Applied Physics, 62(4) 1141-1148. https://doi.org/10.1063/1.339662
%
% For more info check documentation.

  persistent firstBesselZero;
  persistent neutralSpeciesIDs;
  persistent reducedMass;
  persistent sigma;
  persistent epsilon;
  persistent gammaSumArray;
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
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % find IDs of volume neutral species that are going to be taken into account for the multicomponent transport
    % (childless neutral volume species)
    neutralSpeciesIDs = [];
    for i = 1:length(chemistry.stateArray)
      if isempty(chemistry.stateArray(i).ionCharg) && isempty(chemistry.stateArray(i).childArray) && chemistry.stateArray(i).isVolumeSpecies
        neutralSpeciesIDs(end+1) = i;
      end
    end
    % evaluate parameters needed for the calculation of the inverse of the binary diffusion coefficients
    reducedMass = zeros(length(chemistry.stateArray));
    sigma = zeros(length(chemistry.stateArray));
    epsilon = zeros(length(chemistry.stateArray));
    for idx = 1:length(neutralSpeciesIDs)
      for jdx = idx+1:length(neutralSpeciesIDs)
        i = neutralSpeciesIDs(idx);
        j = neutralSpeciesIDs(jdx);
        % error checking
        if isempty(chemistry.stateArray(i).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(j).gas.name);
        elseif isempty(chemistry.stateArray(i).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(j).gas.name);
        elseif isempty(chemistry.stateArray(i).gas.lennardJonesDepth)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDepth'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.lennardJonesDepth)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDepth'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(j).gas.name);
        end
        % evaluation of different auxiliary variables needed
        reducedMass(i,j) = chemistry.stateArray(i).gas.mass*chemistry.stateArray(j).gas.mass/(chemistry.stateArray(i).gas.mass+chemistry.stateArray(j).gas.mass);
        sigma(i,j) = 0.5*(chemistry.stateArray(i).gas.lennardJonesDistance+chemistry.stateArray(j).gas.lennardJonesDistance);
        epsilon(i,j) = sqrt(chemistry.stateArray(i).gas.lennardJonesDepth*chemistry.stateArray(j).gas.lennardJonesDepth);
        reducedMass(j,i) = reducedMass(i,j);
        sigma(j,i) = sigma(i,j);
        epsilon(j,i) = epsilon(i,j);
      end
    end
    % initialize array of gammaSum (total loss wall probability)
    gammaSumArray = zeros(size(chemistry.stateArray));
  end

  % evaluate total wall reaction probability for the reactant specie (only done once per lost species)
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
  end

  % --- regular calculations ---

  % evaluate total neutral density (neutralSpeciesIDs only include childless neutral species)
  neutralDens = sum(densitiesAll(neutralSpeciesIDs));

  % evaluation of the inverse of the reduced binary diffusion coeficients
  inverseReducedBinaryDiffCoeff = zeros(size(densitiesAll));
  for i = neutralSpeciesIDs
    if i == volumeReactantID
      continue;
    end
    Taux = chemistry.workCond.gasTemperature*Constant.boltzmann/epsilon(i,volumeReactantID);
    % fit to first order collision integral from (R.J. Kee, M.E. Coltrin, P. Glarborg "Chemically Reacting Flow:
    % Theory and Practice" John Wiley (2003), Pag. 492)
    collIntegral = 1.0548*Taux^-0.15504+(Taux+0.55909)^-2.1705;
    % actual evaluation of the inverse of the reduced binary diffusion coeficient between species i and j
    inverseReducedBinaryDiffCoeff(i) = 16/3*sqrt(pi*reducedMass(i,volumeReactantID)/(2*Taux*epsilon(i,volumeReactantID)))*...
      sigma(i,volumeReactantID)^2*collIntegral;
  end

  % evaluate effective diffusion coefficient (using Wilke's formula)
  effDiffCoeff = (1-densitiesAll(volumeReactantID)/neutralDens)/sum(densitiesAll.*inverseReducedBinaryDiffCoeff);

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

  % evaluate aux variable
  aux = thermalVelocity*gammaSumArray(volumeReactantID)/(effDiffCoeff*4*(1-gammaSumArray(volumeReactantID)/2));

  % evaluate effective diffusion length for a cylinder (including infinitely long cylinder and slab limiting cases)
  L = chemistry.workCond.chamberLength;
  R = chemistry.workCond.chamberRadius;
  if L == 0        % infinitely long cylinder
    squaredEffectiveDiffusionLength = ( R / ...
      fzero(@(x) aux-besselj(1,x)*x/(besselj(0,x)*R), [0 firstBesselZero-eps(firstBesselZero)]) )^2;
  elseif R == 0    % infinitely wide cylinder (slab)
    squaredEffectiveDiffusionLength = ( L / (2*...
      fzero(@(x) aux-tan(x)*x*2/L, [0 pi/2-eps(pi/2)])) )^2;
  else             % finite cylinder
    squaredEffectiveDiffusionLength = 1 / ( ...
      (fzero(@(x) aux-besselj(1,x)*x/(besselj(0,x)*R), [0 firstBesselZero-eps(firstBesselZero)])/R)^2 + ...
      (fzero(@(x) aux-tan(x)*x*2/L, [0 pi/2-eps(pi/2)])*2/L)^2 );
  end

  % evaluate wall reaction coefficient
  gamma = rateCoeffParams{1};

  % evaluate rate coefficient
  rateCoeff = (gamma/gammaSumArray(volumeReactantID))*(effDiffCoeff/squaredEffectiveDiffusionLength);
  rateCoeff = rateCoeff/(chemistry.workCond.surfaceSiteDensity*chemistry.workCond.areaOverVolume);

  % set function dependencies
  dependent = dependentInfo;

end
