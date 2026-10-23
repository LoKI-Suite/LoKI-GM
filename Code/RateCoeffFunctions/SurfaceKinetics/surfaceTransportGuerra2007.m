function [rateCoeff, dependent] = surfaceTransportGuerra2007(~, ~, ~, ~, rateCoeffParams, chemistry)
% surfaceTransportGuerra2007 evaluates a reaction rate coefficient following eqs. 58, 60 and 61 of Guerra 2007
% The rate coefficient is modified to consider the collision of the volume species with ONE surface species 
%  (dividing by surfaceSiteDensity*AreaOverVolume).
%  k_r ... stericFactor
%  E_r ... activationEnergy
%  nu_D ... freqDiff
%  E_D ... energyBarrier
% 
% Guerra V, Tejero-del-Caz A, Pintassilgo C D and Alves L L 2019 Plasma Sources Sci. Technol. 28 073001
% Guerra V 2007 IEEE Trans. Plasma Sci. 35 1397
% Marinov D, Teixeira C and Guerra V 2017Plasma Processes Polym. 14 1600175

% For more info check documentation.
  
  persistent dependentInfo;
  
  % --- performance sensitive calculations ---
  
  % calculations performed once per simulation 
  if isempty(dependentInfo)
    % define dependencies of the rate coefficient
    if isempty(chemistry.workCond.wallTemperature)
      dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
    else
      dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);
    end
  end
    
  % --- regular calculations ---
  
  % evaluate diffusion frequency, energy barrier for diffusion, steric factor and activation energy
  freqDiff = rateCoeffParams{1};
  energyBarrier = rateCoeffParams{2};
  stericFactor = rateCoeffParams{3};
  activationEnergy = rateCoeffParams{4};
  includeBackDiffusionCorrection = rateCoeffParams{5};

  % extra parameters to calculate the rate coefficient, following Guerra 2007
  sFraction = rateCoeffParams{6};
  fFraction = 1-sFraction;
  freqDes = rateCoeffParams{7};
  energyDes = rateCoeffParams{8};

  % select temperature to use
  if ~isempty(chemistry.workCond.wallTemperature)
    T = chemistry.workCond.wallTemperature;
  else
    T = chemistry.workCond.gasTemperature;
  end

  % evaluate rate coefficient, following eqs. 58, 60 and 61 of Guerra 2007
  kDp=(sFraction/fFraction)*...
      ((freqDiff/freqDes)*(exp((energyDes-energyBarrier)/(Constant.idealGas*T)))-0.25);
  k4p=stericFactor*exp(-activationEnergy/(Constant.idealGas*T));
  kDp=min(1.0,max(0.0,kDp));
  k5 = kDp * k4p * freqDes*exp(-energyDes/(Constant.idealGas*T))/...
    (sFraction*chemistry.workCond.surfaceSiteDensity);
  rateCoeff = k5 / chemistry.workCond.areaOverVolume;
  if includeBackDiffusionCorrection
    rateCoeff = rateCoeff*3/4;
  end
  
  % set function dependencies
  dependent = dependentInfo;

end