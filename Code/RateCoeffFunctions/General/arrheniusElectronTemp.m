function [rateCoeff, dependent] = arrheniusElectronTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% arrheniusElectronTemp evaluates a reaction rate coefficient using the expression:
%  C = a * exp(b/Te)

  Te = chemistry.workCond.electronTemperature; % in eV
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};

  rateCoeff = a * exp(b/Te);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end