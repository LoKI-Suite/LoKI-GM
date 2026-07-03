function [rateCoeff, dependent] = surfaceTransport(~, ~, ~, ~, rateCoeffParams, chemistry)
% surfaceTransport evaluates a reaction rate coefficient using the expression:
%  C = k_r * exp(E_r/(RTw)) (3/4) * nu_D exp(E_D/(RTw)) / {([S]_surface+[F]_surface)*(A/V)}
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

  % select temperature to use
  if ~isempty(chemistry.workCond.wallTemperature)
    T = chemistry.workCond.wallTemperature;
  else
    T = chemistry.workCond.gasTemperature;
  end

  % evaluate rate coefficient
  rateCoeff = freqDiff*exp(-energyBarrier/(Constant.idealGas*T))*...
    stericFactor*exp(-activationEnergy/(Constant.idealGas*T))/(chemistry.workCond.surfaceSiteDensity*chemistry.workCond.areaOverVolume);
  if includeBackDiffusionCorrection
    rateCoeff = rateCoeff*3/4;
  end
  
  % set function dependencies
  dependent = dependentInfo;

end