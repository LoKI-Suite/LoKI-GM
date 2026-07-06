function [rateCoeff, dependent] = arrhenius(~, ~, ~, ~, rateCoeffParams, ~)
% arrhenius evaluates a reaction rate coefficient using the expression:
%  C = a * exp(b/T)

  T = rateCoeffParams{1};
  a = rateCoeffParams{2};
  b = rateCoeffParams{3};

  rateCoeff = a * exp(b/T);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);

end