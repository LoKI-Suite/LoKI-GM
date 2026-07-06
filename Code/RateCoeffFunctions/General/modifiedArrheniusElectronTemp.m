function [rateCoeff, dependent] = modifiedArrheniusElectronTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% modifiedArrheniusElectronTemp evaluates a reaction rate coefficient using the expression:
%  C = a * Te^b * exp(c/Te)

  Te = chemistry.workCond.electronTemperature; % in eV
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};

  rateCoeff = a * Te^b * exp(c/Te);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end