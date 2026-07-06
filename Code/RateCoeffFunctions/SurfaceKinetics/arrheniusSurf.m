function [rateCoeff, dependent] = arrheniusSurf(~, ~, ~, ~, rateCoeffParams, chemistry)
% arrheniusSurf evaluates a reaction rate coefficient using the expression:
%  C = a * exp(b/Tg) / {([S]_surface+[F]_surface)*(A/V)}
% The rate coefficient is modified to consider the collision of the volume species with ONE surface species 
%  (dividing by surfaceSiteDensity*AreaOverVolume).

  Tg = chemistry.workCond.gasTemperature;   %This should be surface temperature
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
   
  rateCoeff = a * exp(b/Tg)/(chemistry.workCond.surfaceSiteDensity*chemistry.workCond.areaOverVolume);
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end