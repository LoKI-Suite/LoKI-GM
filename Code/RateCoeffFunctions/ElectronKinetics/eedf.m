function [rateCoeff, dependent] = eedf(~, ~, ~, reaction, ~, ~)
% eedf evaluates a reaction rate coefficient by integrating an electron-impact cross section over an eedf 
%  (either prescribed or obtained from the solution of the Boltzmann equation)

  rateCoeff = reaction.eedfEquivalent.ineRateCoeff;
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end