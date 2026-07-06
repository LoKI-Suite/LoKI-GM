function [rateCoeff, dependent] = classicalAmbipolarDiffNegIon(time, densitiesAll, totalGasDensity, reaction, ~, chemistry)
% classicalAmbipolarDiffNegIon evaluates the diffusion rate coefficient of a particular positive ion, 
%  for a plasma with multiple positive ions under the classical ambipolar diffusion approximation,
%  considering the effect of several negative ions

% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Rogoff G L 1985 J. Phys. D: Appl. Phys. 18 1533-45
% Ferreira C M, Gousset G and Touzeau M 1988 J. Phys. D: Appl. Phys. 21 1403-13

% Generalized from:
% 
% Guerra V and J Loureiro J 1999 Plasma Sources Sci. Technol. 8 110-124
% Main approximations:
% - pressure sufficiently high, so as to verify the classical ambipolar regime
% - Cold ion approximation: T_char_ions << T_char_electrons
% - Ion mobilities much smaller than the electron mobility
% - Infinite cylinder
% Made by T C Dias and V Guerra (April 2022)

% For more info check documentation.

  persistent firstBesselZero;
  persistent positiveIonIDs;
  persistent negativeIonIDs;
  persistent electronCharEnergy;
  persistent characteristicLengthSquared;
  persistent previousTime;
  persistent rateCoeffs;
  persistent rateCoeffIsTimeDep;
  persistent ionReactions;

  % --- performance sensitive calculations ---

  % initialize variables the first time the classicalAmbipolar function is called
  if isempty(firstBesselZero)
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % evaluate the IDs of the positive ions (singly ionized) and check dependencies on gas temperature
    for i = 1:length(chemistry.stateArray)
      state = chemistry.stateArray(i);
      if strcmp(state.ionCharg, '+') && state.isVolumeSpecies
        positiveIonIDs(end+1) = state.ID;
        ionReactions = [ionReactions state.reactionsDestruction state.reactionsCreation];
      elseif strcmp(state.ionCharg, '-') && state.isVolumeSpecies
        negativeIonIDs(end+1) = state.ID;
        ionReactions = [ionReactions state.reactionsDestruction state.reactionsCreation];
      end
    end

    % eliminate repeated ion reactions
    ionReactions = unique(ionReactions);
    % eliminate diffusion reactions
    i = 1;
    while i <= length(ionReactions)
      if ionReactions(i).isTransport
        ionReactions(i) = [];
      else
        i = i + 1;
      end
    end

    % initialize the rate-coefficients array
    rateCoeffs = zeros(size(chemistry.reactionArray));

    % initialize previousTime variable
    previousTime = -1;
  end

  if time == 0
    rateCoeffIsTimeDep = ones(size(chemistry.reactionArray));
  end

  % --- calculations that are equal for every "classicalAmbipolarDiffNegIon" reaction (at a given time) ---

  if previousTime ~= time
    % actualize previousTime variable
    previousTime = time;
    % obtain electron data
    electronDensity = chemistry.workCond.electronDensity;
    electronRedDiffCoeff = chemistry.electronTransportProperties.reducedDiffCoeff;
    electronRedMobCoeff = chemistry.electronTransportProperties.reducedMobility;
    electronCharEnergy = electronRedDiffCoeff/electronRedMobCoeff;

    % ----- calculate the correction in the characteristic length due to  ---- %
    % -----       the presence of the negative ion                        ---- %

    if ~isempty(negativeIonIDs)
      % calculate all the rate coefficients that will be needed (the ones
      % that involve an ion)
      for rxn = ionReactions
        if rateCoeffIsTimeDep(rxn.ID)
          [rateCoeffs(rxn.ID), dep] = rxn.rateCoeffFuncHandle(time, densitiesAll, totalGasDensity, rxn, rxn.rateCoeffParams, chemistry);
          rateCoeffIsTimeDep(rxn.ID) = dep.onTime || dep.onDensities || dep.onGasTemperature;
        end
      end

      % calculate the total density of negative ions
      totalNegIonDensity = sum(densitiesAll(negativeIonIDs));

      % calculate the sum of the attachment (detachment) frequencies, weighted by the
      % inverse of the negative-ion mobilities
      totalAttachFrequencyWeighted = 0;
      totalDetachFrequencyWeighted = 0;
      for negIonID = negativeIonIDs
        negIon = chemistry.stateArray(negIonID);
        negIonDensity = densitiesAll(negIonID);
        % attachment
        if isempty(negIon.evaluateReducedMobility(chemistry.workCond))
            error(['The ''reducedMobility'' of the ion ''%s'' was not found.'],negIon.name);
        end
        attachmentDensityRate = 0;
        for rxn = negIon.reactionsCreation
          partialAttachmentDensityRate = rateCoeffs(rxn.ID)*electronDensity^rxn.reactantElectrons;
          for i = 1:length(rxn.reactantArray)
            partialAttachmentDensityRate = partialAttachmentDensityRate*...
              densitiesAll(rxn.reactantArray(i).ID)^rxn.reactantStoiCoeff(i);
          end
          for i = 1:length(rxn.catalystArray)
            partialAttachmentDensityRate = partialAttachmentDensityRate*...
              densitiesAll(rxn.catalystArray(i).ID)^rxn.catalystStoiCoeff(i);
          end
          attachmentDensityRate = attachmentDensityRate + partialAttachmentDensityRate;
        end
        totalAttachFrequencyWeighted = totalAttachFrequencyWeighted + attachmentDensityRate/electronDensity/negIon.evaluateReducedMobility(chemistry.workCond);
        % detachment
        if negIonDensity == 0
          continue;
        end
        detachmentDensityRate = 0;
        for rxn = negIon.reactionsDestruction
          partialDetachmentDensityRate = rateCoeffs(rxn.ID)*electronDensity^rxn.reactantElectrons;
          for i = 1:length(rxn.reactantArray)
            partialDetachmentDensityRate = partialDetachmentDensityRate*...
              densitiesAll(rxn.reactantArray(i).ID)^rxn.reactantStoiCoeff(i);
          end
          for i = 1:length(rxn.catalystArray)
            partialDetachmentDensityRate = partialDetachmentDensityRate*...
              densitiesAll(rxn.catalystArray(i).ID)^rxn.catalystStoiCoeff(i);
          end
          detachmentDensityRate = detachmentDensityRate + partialDetachmentDensityRate;
        end
        totalDetachFrequencyWeighted = totalDetachFrequencyWeighted + detachmentDensityRate/negIonDensity/negIon.reducedMobility;
      end

      % calculate the sum of the net ionization frequencies, weighted by the
      % inverse of the positive-ion mobilities
      totalIonizFrequencyWeighted = 0;
      for posIonID = positiveIonIDs
        ionizDensityRate = 0;
        posIon = chemistry.stateArray(posIonID);
        for rxn = posIon.reactionsCreation
          if rxn.isTransport
            continue;
          end
          partialIonizDensityRate = rateCoeffs(rxn.ID)*electronDensity^rxn.reactantElectrons;
          for i = 1:length(rxn.reactantArray)
            partialIonizDensityRate = partialIonizDensityRate*...
              densitiesAll(rxn.reactantArray(i).ID)^rxn.reactantStoiCoeff(i);
          end
          for i = 1:length(rxn.catalystArray)
            partialIonizDensityRate = partialIonizDensityRate*...
              densitiesAll(rxn.catalystArray(i).ID)^rxn.catalystStoiCoeff(i);
          end
          ionizDensityRate = ionizDensityRate + partialIonizDensityRate;
        end
        for rxn = posIon.reactionsDestruction
          if rxn.isTransport
            continue;
          end
          partialIonizDensityRate = -rateCoeffs(rxn.ID)*electronDensity^rxn.reactantElectrons;
          for i = 1:length(rxn.reactantArray)
            partialIonizDensityRate = partialIonizDensityRate*...
              densitiesAll(rxn.reactantArray(i).ID)^rxn.reactantStoiCoeff(i);
          end
          for i = 1:length(rxn.catalystArray)
            partialIonizDensityRate = partialIonizDensityRate*...
              densitiesAll(rxn.catalystArray(i).ID)^rxn.catalystStoiCoeff(i);
          end
          ionizDensityRate = ionizDensityRate + partialIonizDensityRate;
        end
        totalIonizFrequencyWeighted = totalIonizFrequencyWeighted + ionizDensityRate/electronDensity/posIon.evaluateReducedMobility(chemistry.workCond);
      end

      % parameter measuring the intensity of attachment
      P = totalAttachFrequencyWeighted/totalIonizFrequencyWeighted;

      % parameter measuring the intensity of detachment
      Q = totalDetachFrequencyWeighted/totalIonizFrequencyWeighted;

      % use an exponential fit to calculate the slope of lambda(P)
      % The fit was obtained so as to reproduce the results of figure 6 of Ferreira 1988
      a1 = 2.614; t1 = 0.5662;
      a2 = 1.34295; t2 = 9.34383;
      a3 = 0.114319; t3 = 139.119;
      slope = a1*exp(-Q/t1)+a2*exp(-Q/t2)+a3*exp(-Q/t3);

      lambda = slope*P + firstBesselZero^2;

      % at initial stages of the simulation, there can be a negative net
      % production of positive ions, which would lead to negative parameters P,Q
      % In these cases, we use the classical solution, with no influence of
      % negative ions
      if totalIonizFrequencyWeighted < 0 || totalAttachFrequencyWeighted < 0 || totalDetachFrequencyWeighted < 0
        characteristicLengthSquared = (chemistry.workCond.chamberRadius/firstBesselZero)^2;
      else
        characteristicLengthSquared = chemistry.workCond.chamberRadius^2/lambda/(1-totalNegIonDensity/(totalNegIonDensity+electronDensity));
      end
    else
      % evaluate the squared characteristic length
      if chemistry.workCond.chamberLength == 0
        characteristicLengthSquared = (chemistry.workCond.chamberRadius/firstBesselZero)^2;
      elseif chemistry.workCond.chamberRadius == 0
        characteristicLengthSquared = (chemistry.workCond.chamberLength/pi)^2;
      else
        characteristicLengthSquared = 1/((chemistry.workCond.chamberRadius/firstBesselZero)^-2+(chemistry.workCond.chamberLength/pi)^-2);
      end
    end
  end

  % --- regular calculations ---

  % evaluate diffusion rate coefficient
  redMobility = reaction.reactantArray.evaluateReducedMobility(chemistry.workCond);
  if isempty(redMobility)
    error(['Error found when evaluating the classicalAmbipolarDiffNegIon rate-coefficient function for\n%s\n'...
          'The ''reducedMobility'' of the state ''%s'' was not found.'],...
          reaction.descriptionExtended, reaction.reactantArray.name);
  end
  rateCoeff = redMobility*electronCharEnergy/(characteristicLengthSquared*totalGasDensity);

  if chemistry.workCond.electronDensity == 0 || sum(densitiesAll(positiveIonIDs)) == 0
    rateCoeff = 0;
  end

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', true);
end
