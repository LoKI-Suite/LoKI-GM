function [rateCoeff, dependent] = expTemp(~, ~, ~, ~, rateCoeffParams, ~)
% expTemp evaluates a reaction rate coefficient using the expression:
%  C = a * exp(T/b)

  T = rateCoeffParams{1};
  a = rateCoeffParams{2};
  b = rateCoeffParams{3};

  rateCoeff = a * exp(T/b);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);

end