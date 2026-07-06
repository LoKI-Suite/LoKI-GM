function [rateCoeff, dependent] = powerlnGasTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% powerlnGasTemp evaluates a reaction rate coefficient using the expression:
%  C = (a + b * log(c/Tg) ) * Tg^d

  Tg = chemistry.workCond.gasTemperature; % in K
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};
  d = rateCoeffParams{4};

  rateCoeff = (a + b * log(c/Tg) ) * Tg^d;
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end