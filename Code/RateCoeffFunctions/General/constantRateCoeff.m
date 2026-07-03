function [rateCoeff, dependent] = constantRateCoeff(~, ~, ~, ~, rateCoeffParams, ~)
% constantRateCoeff evaluates a reaction rate coefficient using the expression:
%  C = a

  rateCoeff = rateCoeffParams{1};
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);

end