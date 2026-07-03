function [rateCoeff, dependent] = classicalAmbipolarDiff(time, densitiesAll, totalGasDensity, reaction, ...
  rateCoeffParams, chemistry)
% classicalAmbipolarDiff evaluates the diffusion rate coefficient of a particular positive ion, 
%  for a plasma with multiple positive ions, under the classical ambipolar diffusion approximation

% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Coche P, Guerra V and Alves L L 2016 J. Phys. D: Appl. Phys. 49 235207
% Rogoff G L 1985 J. Phys. D: Appl. Phys. 18 1533-45

% For more info check documentation.
  
  persistent firstBesselZero;
  persistent positiveIonIDs;
  persistent initialTime;
  persistent characteristicLengthSquared;
  persistent correctionFactor;
  persistent dependentInfo;
  
  % --- performance sensitive calculations ---
  
  % initialize variables the first time the classicalAmbipolar function is called
  if isempty(firstBesselZero)
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % save initial time of the simulation
    initialTime = time;
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', false, 'onElectronKinetics', false);
    % evaluate the IDs of the positive ions (singly ionized and gas phase) and check dependencies on gas temperature
    for state = chemistry.stateArray
      if strcmp(state.ionCharg, '+') && state.isVolumeSpecies
        positiveIonIDs(end+1) = state.ID;
        if ~dependentInfo.onGasTemperature && ~isempty(state.reducedDiffCoeffFunc)
          for param = state.reducedDiffCoeffParams
            if strcmp(param{1}, 'gasTemperature')
              dependentInfo.onGasTemperature = true;
              break;
            end
          end
        end
      end
    end
  end
  
  % --- time independent calculations (also performance sensitive) ---
  workCond = chemistry.workCond;
  if time == initialTime
    % evaluate the squared characteristic length
    if workCond.chamberLength == 0
      characteristicLengthSquared = (workCond.chamberRadius/firstBesselZero)^2;
    elseif workCond.chamberRadius == 0
      characteristicLengthSquared = (workCond.chamberLength/pi)^2;
    else
      characteristicLengthSquared = 1/((workCond.chamberRadius/firstBesselZero)^-2+(workCond.chamberLength/pi)^-2);
    end
  end
  
  % --- calculations that are equal for every "classicalAmbipolar" reaction (at a given time) ---
  
  % evaluate weighted sums of the diffusion coefficient and reduced mobility of all positive ions (singly ionized)
  ionRedDiffCoeffSum = 0;
  ionRedMobCoeffSum = 0;
  for i = 1:length(positiveIonIDs)
    ionID = positiveIonIDs(i);
    ionDensity = densitiesAll(ionID);
    ionRedDiffCoeffSum = ionRedDiffCoeffSum + ionDensity*chemistry.stateArray(ionID).evaluateReducedDiffCoeff(workCond);
    ionRedMobCoeffSum = ionRedMobCoeffSum + ionDensity*chemistry.stateArray(ionID).evaluateReducedMobility(workCond);
  end
  % obtain electron data
  electronDensity = workCond.electronDensity;
  electronRedDiffCoeff = chemistry.electronTransportProperties.reducedDiffCoeff;
  electronRedMobCoeff = chemistry.electronTransportProperties.reducedMobility;
  % evaluate correction factor
  correctionFactor = (ionRedDiffCoeffSum-electronDensity*electronRedDiffCoeff)/...
    (ionRedMobCoeffSum+electronDensity*electronRedMobCoeff);

  
  % --- regular calculations ---
  
  % evaluate diffusion rate coefficient
  if workCond.electronDensity==0 && all(densitiesAll(positiveIonIDs)==0)
    rateCoeff = 0;
  else
    wallProbability = rateCoeffParams{1};
    redDiffCoeff = reaction.reactantArray.evaluateReducedDiffCoeff(workCond);
    redMobility = reaction.reactantArray.evaluateReducedMobility(workCond);
    rateCoeff = (redDiffCoeff-redMobility*correctionFactor)/(characteristicLengthSquared*totalGasDensity);
    rateCoeff = wallProbability*rateCoeff;
  end
  % set function dependencies
  dependent = dependentInfo;

end
