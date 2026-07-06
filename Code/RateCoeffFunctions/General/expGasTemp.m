function [rateCoeff, dependent] = expGasTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% expGasTemp evaluates a reaction rate coefficient using the expression:
%  C = a * exp(Tg/b)

  Tg = chemistry.workCond.gasTemperature; % in K
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};

  rateCoeff = a * exp(Tg/b);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end