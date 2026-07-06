function [rateCoeff, dependent] = thermalDesorption(~, ~, ~, ~, rateCoeffParams, chemistry)
% thermalDesorption evaluates a reaction rate coefficient using the expression:
%  C = nu_d exp(E_d/(RTw))
%  E_d ... activationEnergy
%  nu_d ... freqFactor
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
  
  % evaluate frequency factor and activation energy
  freqFactor = rateCoeffParams{1};
  activationEnergy = rateCoeffParams{2};

  % select temperature to use
  if ~isempty(chemistry.workCond.wallTemperature)
    T = chemistry.workCond.wallTemperature;
  else
    T = chemistry.workCond.gasTemperature;
  end

  % evaluate rate coefficient
  rateCoeff = freqFactor*exp(-activationEnergy/(Constant.idealGas*T));
  
  % set function dependencies
  dependent = dependentInfo;

end