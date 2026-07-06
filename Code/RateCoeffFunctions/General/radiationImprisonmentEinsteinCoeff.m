function [rateCoeff, dependent] = radiationImprisonmentEinsteinCoeff(~, densitiesAll, ~, reaction, rateCoeffParams, chemistry)
  % radiationImprisonmentEinsteinCoeff evaluates the Einstein coefficient multiplicated by the escape factor, using the expression:
  %  C = Aij*lambda
  % M.Santos et al 2014 J. Phys. D: Appl. Phys, 47, 265201
  
  % values that will not be updated at each time
  persistent  currentReactionID;        % identifies the transition position in the arrays used
  persistent  einsteinProbability;      % stores the einstein coefficients of each transition
  persistent  amplitudeVoigtProfileAux; % stores the independent factors in the Voigt amplitude profile
  persistent  dampingCoeff;             % stores the independent factors of the damping coefficient
  persistent  frequencyTransition;      % stores the photon frequency of a transitions
  persistent  dependentInfo;
  
  kb = Constant.boltzmann; % K
  Tg = chemistry.workCond.gasTemperature; % K
  mass = reaction.reactantArray(1).mass; % kg
  e = Constant.electronCharge; % C
  
  h = Constant.planck; % J s
  me = Constant.electronMass; % kg
  c = Constant.speedOfLight; % m s-1
  R = chemistry.workCond.chamberRadius; % m
  mu0 = Constant.vacuumPermeability; % SI
  
  % initiates the arrays the first time that this function is called
  if isempty(currentReactionID)
    currentReactionID = [];
    einsteinProbability = [];
    amplitudeVoigtProfileAux = [];
    dampingCoeff = [];
    frequencyTransition = [];
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
  end
  
  % defined the transition ID as zero
  localID = 0;
  
  % assignes the reaction ID to the localID variable
  % this will only happen after this function is called for every transition
  % in which imprisonment effects are included
  for i=1:length(currentReactionID)
    if currentReactionID(i) == reaction.ID
      localID = i;
    end
  end
  
  % calculates time independent quantities and stores them in the respective arrays
  % this will only happen the first time that this function is called for a given transition
  if localID==0
    Aij = rateCoeffParams{1}; % Einstein coefficients for radiative transitions
  
    ui = reaction.reactantArray.energy;
    gi = reaction.reactantArray.statisticalWeight;
    uj = reaction.productArray.energy;
    gj = reaction.productArray.statisticalWeight;
  
    du_ij = ui - uj; % eV
  
    % photon frequency
    nuij = du_ij*e/h; % s s-1
  
    % common factor between AAijaux and k0aux
    C = Aij*(c/nuij)^3*gi/gj;
  
    currentReactionID(end+1) = reaction.ID;
    frequencyTransition(end+1) = nuij;
    AAijaux = C/(6*pi^2);
    dampingCoeff(end+1) = AAijaux;
    k0aux = C/(8*sqrt(pi)*pi);
    amplitudeVoigtProfileAux(end+1) =  k0aux;
    einsteinProbability(end+1) = Aij;
  else
    k0aux = amplitudeVoigtProfileAux(localID);
    Aij = einsteinProbability(localID);
    AAijaux = dampingCoeff(localID);
    nuij = frequencyTransition(localID);
  end

  % most probable speed of the Maxwellian distribution with the atomic species
  v0 = sqrt(2*kb*Tg/mass);
  nj = densitiesAll(reaction.productArray.ID);

  % Voigt spectral profile amplitude
  k0 = k0aux * nj / v0;
  
  k0R = k0*R;
  
  % damping coefficient
  AAij = AAijaux * nj;
  % absorption coefficient
  a = c *(Aij + AAij)/(4 * pi)/( nuij*v0);
  
  % calculation of partial escape factors
  Ld = 1.60/(k0R*sqrt(pi*log(k0R)));
  Lc = 2/sqrt(pi)*sqrt( sqrt(pi)*a / (pi*k0R) ) ;
  Lcd = 2*a/(pi*sqrt(log(k0R)));
  
  % in the article of M.Santos, it is said that k0R > 8 is a necessary
  % condition for the use of these formulas
  if k0R<8
      lambda = 1;
  else
      lambda = Ld * exp(-Lcd^2/Lc^2) + Lc * erf(Lcd/Lc);
  end
  
  rateCoeff = Aij*lambda;
  
  dependent = dependentInfo; 
  
end
