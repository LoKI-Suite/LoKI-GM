function [rateCoeff, dependent] = offsetArrheniusGasTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% offsetArrheniusGasTemp evaluates a reaction rate coefficient using the expression:
%  C = a0 + a * exp(b/Tg)

  Tg = chemistry.workCond.gasTemperature;
  a0 = rateCoeffParams{1};
  a = rateCoeffParams{2};
  b = rateCoeffParams{3};

  rateCoeff = a0 + a * exp(b/Tg);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end
