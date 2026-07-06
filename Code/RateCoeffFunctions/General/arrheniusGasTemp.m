function [rateCoeff, dependent] = arrheniusGasTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% arrheniusGasTemp evaluates a reaction rate coefficient using the expression:
%  C = a * exp(b/Tg)

  Tg = chemistry.workCond.gasTemperature;
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};

  rateCoeff = a * exp(b/Tg);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end