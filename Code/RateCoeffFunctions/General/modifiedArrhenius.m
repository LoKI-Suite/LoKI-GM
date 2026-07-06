function [rateCoeff, dependent] = modifiedArrhenius(~, ~, ~, ~, rateCoeffParams, ~)
% modifiedArrhenius evaluates a reaction rate coefficient using the expression:
%  C = a * T^b * exp(c/T)

  T = rateCoeffParams(1);
  a = rateCoeffParams(2);
  b = rateCoeffParams(3);
  c = rateCoeffParams(4);

  rateCoeff = a * T^b * exp(c/T);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);

end