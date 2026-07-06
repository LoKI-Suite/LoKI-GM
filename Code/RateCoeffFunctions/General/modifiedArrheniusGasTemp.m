function [rateCoeff, dependent] = modifiedArrheniusGasTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% modifiedArrheniusGasTemp evaluates a reaction rate coefficient using the expression:
%  C = a * Tg^b * exp(c/Tg)

  Tg = chemistry.workCond.gasTemperature;
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};

  rateCoeff = a * Tg^b * exp(c/Tg);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end