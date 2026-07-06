function [rateCoeff, dependent] = powerGasElecTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% powerGasElecTemp evaluates a reaction rate coefficient using the expression:
%  C = a * Tg^b * Te^c

  Tg = chemistry.workCond.gasTemperature; % in K
  Te = chemistry.workCond.electronTemperature; %in eV
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};

  rateCoeff = a * Tg^b * Te^c;
  
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', true);

end
