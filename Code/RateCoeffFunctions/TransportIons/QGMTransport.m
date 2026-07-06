function [rateCoeff, dependent] = QGMTransport(~, densitiesAll, ~, reaction, rateCoeffParams, chemistry)
% QGMTransport evaluates the transport rate coefficient of a particular positive ion
%  for a plasma with multiple positive and negative ions (assumed singly-charged)
% - assuming that the electrons and the negative ions are in Boltzmann equilibrium with the space-charge potential,
%  thus considering (grad n_neg) = alpha * gamma_neg (grad ne)
% - moderate electronegativities, satisfying alpha(mu_pos,neg / mu_e) << 1
% - low ion mobilities, satisfying mu_pos ~ mu_neg << mu_e
% - with similar temperatures, T_pos ~ T_neg
% using the transport model described in
%  E. Stoffels, et al, Contributions to Plasma Physics, 35 (1995) 331-357
% and adopted by the Quantemol Global Model (QGM)
%  J. Tennyson et al. 2022 Plasma Sources Sci. Technol. 31 095020

% The full model (case "transport") considers not only the effects of the ambipolar transport but also the thermal losses at the wall,
%  described by characteristic frequencies estimated from the volume average of the corresponding particle balance equation
%  In this case, the free diffusion coefficient of the positive ion species is weighted by the population of these species
% The model can also run in a simplified version (case "diffusion"), in which case it considers only the effects of
%  the ambipolar transport, using the binary diffusion coefficient as free diffusion coefficient of the positive ion species
% The rateCoeffParams are
%  - (1) the wall recombination probability of the positive ion
%  - (2) a boolean operator (true/false) to include/exclude the effect in the diffusion coefficient
%        of the positive ion - (positive/negative)ion collisions
%  - (3) a string ("transport"/"diffusion") to select the full/simplified version of the model
%
% For more info check
%  LL Alves and A Tejero-del-Caz, Plasma Sources Sci. Technol., 32 (2023) 054003

  persistent neutralSpeciesIDs;
  persistent positiveIonIDs;
  persistent negativeIonIDs;
  persistent charDiffLength;
  persistent sigma;
  persistent reducedMass;
  persistent dependentInfo;

  % --- performance sensitive calculations ---

  % calculations performed only once per simulationm when this function is called for the first time
  if isempty(dependentInfo)
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % evaluate geometrical parameters for a cylinder (including infinitely long cylinder and slab limiting cases)
    if chemistry.workCond.chamberLength == 0        % infinitely long cylinder
      charDiffLength = chemistry.workCond.chamberRadius/firstBesselZero;
    elseif chemistry.workCond.chamberRadius == 0    % infinitely wide cylinder (slab)
      charDiffLength = chemistry.workCond.chamberLength/pi;
    else                                  % finite cylinder
      charDiffLength = 1/sqrt((firstBesselZero/chemistry.workCond.chamberRadius)^2+(pi/chemistry.workCond.chamberLength)^2);
    end
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % find the IDs of the (volume) positive and negative ions (singly ionized) as well as neutral (volume) species
    positiveIonIDs = [];
    negativeIonIDs = [];
    neutralSpeciesIDs = [];
    for i = 1:length(chemistry.stateArray)
      state = chemistry.stateArray(i);
      if strcmp(state.ionCharg, '+') && state.isVolumeSpecies
        positiveIonIDs(end+1) = state.ID;
      elseif strcmp(state.ionCharg, '-') && state.isVolumeSpecies
        negativeIonIDs(end+1) = state.ID;
      elseif isempty(chemistry.stateArray(i).ionCharg) && isempty(chemistry.stateArray(i).childArray) && chemistry.stateArray(i).isVolumeSpecies
        neutralSpeciesIDs(end+1) = i;
      end
    end
    % evaluate parameters needed for the calculation of the inverse of the binary diffusion coefficients
    reducedMass = zeros(length(chemistry.stateArray));
    sigma = zeros(length(chemistry.stateArray));
    for i = positiveIonIDs
      % parameters for positive_ion-neutral collisions
      for j = neutralSpeciesIDs
        % error checking
        if isempty(chemistry.stateArray(i).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(j).gas.name);
        end
        % evaluation of the hard-spheres cross section for positive_ion-neutral collisions
        sigma(i,j) = (chemistry.stateArray(i).gas.lennardJonesDistance+chemistry.stateArray(j).gas.lennardJonesDistance)^2;
      end
      % parameters for ion-ion collisions
      for j = [positiveIonIDs negativeIonIDs]
        if i==j
          continue;
        end
        % error checking
        if isempty(chemistry.stateArray(i).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(j).gas.name);
        end
        % evaluation of the reduced mass of the colliding particles
        reducedMass(i,j) = chemistry.stateArray(i).gas.mass*chemistry.stateArray(j).gas.mass/(chemistry.stateArray(i).gas.mass+chemistry.stateArray(j).gas.mass);
      end
    end
  end

  % --- regular calculations ---

  % local save of the ID of the reactant species
  reactantID = reaction.reactantArray(1).ID;

  % check that the reactant species is a (volume) positive ion
  if ~strcmp(chemistry.stateArray(reactantID).ionCharg, '+') || chemistry.stateArray(reactantID).isSurfaceSpecies
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nReactant %s is not a (volume) ' ...
      'positive ion.\nPlease, fix the problem and run the code again.'], reaction.type, ...
      reaction.description, chemistry.stateArray(reactantID).name);
  end

  % Calculation of the ions temperature (both positive and negative) and gamma parameter
  if chemistry.workCond.gasPressure > 0.133
    ionTemperature = (5800 - chemistry.workCond.gasTemperature) * 0.133/chemistry.workCond.gasPressure + chemistry.workCond.gasTemperature;
  else
    ionTemperature = 5800;
  end
  electronTemperatureInKelvin = chemistry.workCond.electronTemperature / Constant.boltzmannInEV;
  gamma = electronTemperatureInKelvin/ionTemperature;

  % evaluate electronegativity
  alpha = sum(densitiesAll(negativeIonIDs))/chemistry.workCond.electronDensity;

  % evaluate density-weighted diffusion coefficient for positive ions
  weightedDiffCoeff = 0;
  includeIonCorrectionDiffCoeff = rateCoeffParams{2}; % boolean rate coefficient parameter (true or false)
  if includeIonCorrectionDiffCoeff
    lambdaDebye = sqrt(Constant.vacuumPermittivity*chemistry.workCond.electronTemperature/ ...
      (Constant.electronCharge*chemistry.workCond.electronDensity));
  end

  modelCase = rateCoeffParams{3}; % string rate coefficient parameter ("transport" or "diffusion")
  if strcmp(modelCase, 'diffusion')
      positiveIonValues = reactantID;
  elseif strcmp(modelCase, 'transport')
      positiveIonValues = positiveIonIDs;
  end
  for i = positiveIonValues
    % evaluate the thermal velocity for the ion "i"
    thermalVelocity = sqrt(8*Constant.boltzmann*ionTemperature/(pi*chemistry.stateArray(i).gas.mass));
    % evaluate the inverse of the mean free path for the ion "i"
    inverseLambda = 0;
    % evaluate contribution to the mean free path from ion-neutral collisions
    for j = neutralSpeciesIDs
      inverseLambda = inverseLambda + densitiesAll(j)*sigma(i,j);
    end
    % evaluate contribution to the mean free path from ion-ion collisions
    if includeIonCorrectionDiffCoeff
      for j = [positiveIonIDs negativeIonIDs]
        if j==i
          continue;
        end
        b0 = Constant.electronCharge^2/(2*pi*Constant.vacuumPermittivity*reducedMass(i,j)*thermalVelocity^2);
        sigmaij = pi*b0^2*log(2*lambdaDebye/b0);
        inverseLambda = inverseLambda + densitiesAll(j)*sigmaij;
      end
    end
    % evaluate contribution of ion "i" to the density-weighted diffusion coefficient for positive ions
    if strcmp(modelCase, 'transport')
      weightedDiffCoeff = weightedDiffCoeff + densitiesAll(i)*thermalVelocity/inverseLambda;
    elseif strcmp(modelCase, 'diffusion')
      weightedDiffCoeff = weightedDiffCoeff + thermalVelocity/inverseLambda;
    end
  end

  if strcmp(modelCase, 'transport')
    if weightedDiffCoeff ~= 0
      weightedDiffCoeff = (pi/8)*weightedDiffCoeff/sum(densitiesAll(positiveIonIDs));
    end
  elseif strcmp(modelCase, 'diffusion')
    weightedDiffCoeff = (pi/8)*weightedDiffCoeff;
  end

  % evaluate the ambipolar diffusion coefficient for positive ions
  ambipolarDiffCoeff = weightedDiffCoeff*((1 + gamma*(1+2*alpha))/(1+alpha*gamma));

  % evaluate rate coefficient
  if strcmp(modelCase, 'transport')

    % evaluate the thermal velocity for the reactant species
    thermalVelocity = sqrt(8*Constant.boltzmann*ionTemperature/(pi*chemistry.stateArray(reactantID).gas.mass));

    % evaluate wall probability
    wallProbability = rateCoeffParams{1};

    % rate coefficient for case 'transport'
    rateCoeff = chemistry.workCond.areaOverVolume*ambipolarDiffCoeff*wallProbability/...
        (charDiffLength*wallProbability+4*ambipolarDiffCoeff/thermalVelocity);

  elseif strcmp(modelCase, 'diffusion')

    % rate coefficient for case 'diffusion'
    rateCoeff = ambipolarDiffCoeff/charDiffLength^2;

  end

  % set function dependencies
  dependent = dependentInfo;

end
