function [rateCoeff, dependent] = modifiedArrheniusCharEnergyTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
% modifiedArrheniusCharEnergyTemp evaluates a reaction rate coefficient using the expression:
%  C = a * charEnergy^b * exp(c/charEnergy)

  charEnergy = chemistry.electronTransportProperties.reducedDiffCoeff / ...
    chemistry.electronTransportProperties.reducedMobility;
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};

  rateCoeff = a * charEnergy^b * exp(c/charEnergy);
  
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end
