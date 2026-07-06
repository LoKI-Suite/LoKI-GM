function [rateCoeff, dependent] = powerGasCharEnergyTemp(~, ~, ~, ~, rateCoeffParams, chemistry)
  % power Gas char energy temp evaluates a reaction rate coefficients using the following expression:
  % k = a * Tg^b * charEnergy^c
  % where charEnergy is the characteristic energy of the gas, defined as the ratio between the reduced diffusion coefficient and the reduced mobility of electrons in the gas (in eV) and Tg is the gas temperature (in K)
% power Gas char energy temp evaluates a reaction rate coefficients using the following expression:

  Tg = chemistry.workCond.gasTemperature; % in K
  charEnergy = chemistry.electronTransportProperties.reducedDiffCoeff / ...
    chemistry.electronTransportProperties.reducedMobility / Constant.boltzmannInEV; % in K
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  c = rateCoeffParams{3};

  rateCoeff = a * Tg^b * charEnergy^c;

  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', true);

end
