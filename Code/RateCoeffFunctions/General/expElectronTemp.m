function [rateCoeff, dependent] = expElectronTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% expElectronTemp evaluates a reaction rate coefficient using the expression:
%  C = a * exp(Te/b)

  Te = chemistry.workCond.electronTemperature; % in eV
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};

  rateCoeff = a * exp(Te/b);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end