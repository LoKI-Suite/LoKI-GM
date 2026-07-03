function [rateCoeff, dependent] = powerElectronTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% powerElectronTemp evaluates a reaction rate coefficients using the following expression:

  Te = chemistry.workCond.electronTemperature; % in eV
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};

  rateCoeff = a * Te^b;
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end