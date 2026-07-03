function [rateCoeff, dependent] = hFactorTransport(time, densitiesAll, totalGasDensity, reaction, rateCoeffParams, chemistry)
% hFactorTransport evaluates the transport rate coefficient of a particular positive ion, adopting h-factor models
%  for drift-dominated scenarios in electropositive plasmas, solved numerically at low-pressure
%  for drift-diffusion scenarios, featuring the transition from electropositive to moderate or strong electronegative plasmas
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Godyak V 1986 ed V D Associates (Fall Church, VA: Delphic Associates)
% Lichtenberg A~J, Vahedi V, Lieberman M~A and Rognlien T 1994 J. Appl. Phys. 75 2339
% Lee C and Lieberman M~A 1995 J. Vac. Sci. Technol. A 13 368
% Kim S, Lieberman M~A, Lichtenberg A~J and Gudmundsson J~T 2006 J. Vac. Sci. Technol. A 24 2025
% Thorsteinsson E~G and Gudmundsson J~T 2010 Plasma Sources Sci. Technol. 19 015001
% Thorsteinsson E~G and Gudmundsson J~T 2010 Plasma Sources Sci. Technol. 19 055008
% Chabert P 2016 Plasma Sources Sci. Technol. 25 025010

  persistent firstBesselZero;
  persistent firstOrderBessel;
  persistent kB;
  persistent elecCharge;
  persistent parameterArray;
  persistent L;
  persistent R;
  persistent positiveIonIDs;
  persistent negativeIonIDs;
  persistent recombinationReactionIDs;

  % --- performance sensitive calculations ---

  % initialize variables the first time the hFactorTransport function is called
  if isempty(firstBesselZero)
    % evaluate the first zero of the zero order bessel function and first order bessel function at that point
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    firstOrderBessel = besselj(1,firstBesselZero);
    % save local copies of physical constants
    kB = Constant.boltzmann;
    elecCharge = Constant.electronCharge;
    % parameters for the fit for alphaS as done in Thorsteinsson 2010 (first paper)
    parameterArray = [0.607, 5.555, -11.16, 1.634, 12*1e-3, -107*1e-3];
    % save local copies of chamber dimensions (with error checking)
    L = chemistry.workCond.chamberLength;
    R = chemistry.workCond.chamberRadius;
    % find the IDs of the positive ions and negative ion (singly ionized and gas phase, with error checking)
    for i = 1:length(chemistry.stateArray)
      state = chemistry.stateArray(i);
      if strcmp(state.ionCharg, '+') && state.isVolumeSpecies
        positiveIonIDs(end+1) = state.ID;
      elseif strcmp(state.ionCharg, '-')  && state.isVolumeSpecies
        negativeIonIDs(end+1) = state.ID;
      end
    end
    % find the IDs of positive-negative ion recombination reactions (with error checking)
    for negativeIonID = negativeIonIDs
      for rxn = chemistry.stateArray(negativeIonID).reactionsDestruction
        if ~rxn.isTransport && rxn.reactantElectrons == 0 && rxn.productElectrons == 0
          positiveReactantFound = false;
          negativeReactantFound = false;
          chargedProductFound = false;
          chargedCatalystFound = false;
          for reactant = rxn.reactantArray
            if ~positiveReactantFound && any(reactant.ID == positiveIonIDs)
              positiveReactantFound = true;
            elseif ~negativeReactantFound && any(reactant.ID == negativeIonIDs)
              negativeReactantFound = true;
            elseif positiveReactantFound && any(reactant.ID == positiveIonIDs)
              positiveReactantFound = false;
              negativeReactantFound = false;
              break
            elseif negativeReactantFound && any(reactant.ID == negativeIonIDs)
              positiveReactantFound = false;
              negativeReactantFound = false;
              break
            end
          end
          for product = rxn.productArray
            if ~isempty(product.ionCharg)
              chargedProductFound = true;
              break
            end
          end
          for catalist = rxn.catalystArray
            if ~isempty(catalist.ionCharg)
              chargedCatalystFound = true;
              break
            end
          end
          if positiveReactantFound && negativeReactantFound && ~chargedProductFound && ~chargedCatalystFound
            recombinationReactionIDs(end+1) = rxn.ID;
          end
        end
      end
    end
  end

  % --- regular calculations ---

  % evaluate temperatures (assuming both positive and negative ion temperatures equal to gas temperature)
  elecTemperature = chemistry.workCond.electronTemperature*elecCharge/kB; % in K
  gasTemperature = chemistry.workCond.gasTemperature;                     % in K
  gamma = elecTemperature/gasTemperature;

  % evaluate ion mean free path
  sigma = rateCoeffParams{1};           % total ion collision cross section (constant value)
  lambda = 1/(totalGasDensity*sigma);

  % evaluate ion sound velocity
  ionMass = reaction.reactantArray(1).mass;
  bohmVel = sqrt(kB*elecTemperature/ionMass);

  % check for the need of corrections due to presence of negative ions
  if isempty(negativeIonIDs)
    alpha0 = 0;
    alphaS = 0;
    hc = 0;
  else
    % evaluate total positive ion density
    totalPositiveIonDensity = sum(densitiesAll(positiveIonIDs));
    % evaluate total negative ion density and electronegativities (mean, axis and sheath values)
    totalNegativeIonDensity = sum(densitiesAll(negativeIonIDs));
    alpha = totalNegativeIonDensity/chemistry.workCond.electronDensity;
    alpha0 = 1.5*alpha;
    useFitToAlphaS = rateCoeffParams{3};
    if useFitToAlphaS
      % fit for alphaS from Thorsteinsson 2010 (1st paper) only works for gamma > 10
      if gamma <= 5
        error(['Ratio between electron temperature and gas temperature should be > 10 for high pressure regime ' ...
          'calculations in hFactorTransport']);
      end
      % expression from Thorsteinsson 2010 (first paper)
      rho = abs(alpha + parameterArray(5)*(exp(parameterArray(6)*(gamma-50))-1));
      alphaS = alpha*((parameterArray(1)*erf(parameterArray(2)*rho + parameterArray(3))*...
        exp(-parameterArray(4)/rho^(1.35)))/(exp((gamma-1)/2*gamma - 0.49)));
      if alphaS < 0
        alphaS = 0;
      end
    else
      if alpha > 0
        alphaS = fzero(@(x) x/alpha-exp(0.5*(1+x)*(1-gamma)/(1+x*gamma)), alpha);
      else
        alphaS = 0;
      end
    end

    % evaluate effective recombination rate coefficient (except division by ion densities, both positive and negative)
    effectiveRecombinationRate = 0;
    for ID = recombinationReactionIDs
      rxn = chemistry.reactionArray(ID);
      rateCoeff = rxn.rateCoeffFuncHandle(time, densitiesAll, totalGasDensity, rxn, rxn.rateCoeffParams, chemistry);
      effectiveRecombinationRate = effectiveRecombinationRate + rateCoeff * ...
        prod(densitiesAll([rxn.reactantArray.ID]).^([rxn.reactantStoiCoeff]')) * ...
        prod(densitiesAll([rxn.catalystArray.ID]).^([rxn.catalystStoiCoeff]')) * ...
        totalGasDensity.^rxn.isGasStabilised;
    end
    % evaluate high electronegativity contribution to h-factors
    if effectiveRecombinationRate == 0
      hc = 0;
    else
      effectiveRecombinationRate = effectiveRecombinationRate/(totalPositiveIonDensity*totalNegativeIonDensity);
      meanVel = sqrt(8*kB*gasTemperature/(pi*ionMass));
      nAst = (15/56)*meanVel/(effectiveRecombinationRate*lambda);
      hc = 1/(sqrt(gamma)*(1+sqrt(nAst)*totalPositiveIonDensity/totalNegativeIonDensity^(3/2)));
    end
  end

  % evaluate the ambipolar diffusion coefficient as function of the free reduced diffusion coefficient
  Da = reaction.reactantArray.evaluateReducedDiffCoeff(chemistry.workCond)/totalGasDensity*...
    (1+gamma*(1+2*alphaS))/(1+gamma*alphaS);

  % evaluate axial h-factor (axial edge-to-center positive ion density ratio)
  if L ~= 0
    useLowPressureCorrectionForhL = rateCoeffParams{2};
    if useLowPressureCorrectionForhL
      % expression from Chabert 2016 (including corrections for low pressure)
      hL = sqrt(((gamma-1)/((1+alpha0)^2)+1)/gamma)*0.86/sqrt(3+L/(2*lambda)+sqrt(1+alpha0)*(L/lambda)^2/(5*gamma));
    else
      % expression from Thorsteinsson 2010 (second paper)
      h0L = 0.86/sqrt(3 + L/(2*lambda) + (0.86*L*bohmVel/(pi*Da))^2);
      hL = sqrt((h0L/(1+alpha0))^2 + hc^2);
    end
    hLoL = hL/L;
  else
    hLoL = 0;
  end

  % evaluate radial h-factor (radial edge-to-center positive ion density ratio)
  if R~=0
    % expression from Thorsteinsson 2010 (second paper)
    hR0 = 0.8/sqrt(4 + R/lambda + (0.8*R*bohmVel/(firstBesselZero*firstOrderBessel*Da))^2);
    hR = sqrt((hR0/(1+alpha0))^2 + hc^2);
    hRoR = hR/R;
  else
    hRoR = 0;
  end

  % evaluate de rate coefficient (loss frequency)
  rateCoeff = 2*bohmVel*(hLoL + hRoR);

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', true);
end
