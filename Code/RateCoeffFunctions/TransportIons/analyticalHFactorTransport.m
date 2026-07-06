function [rateCoeff, dependent] = analyticalHFactorTransport(~, ~, totalGasDensity, reaction, rateCoeffParams, chemistry)
% analyticalHFactorTransport evaluates the transport rate coefficient of a particular positive ion 
%  for a plasma with a single type of positive ions, adopting an h-factor drift-dominated model, 
%  solved analytically for plane-parallel electropositive discharges
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Czarnetzki, U and Alves, L L 2022 Reviews of Modern Plasma Physics, 6(1) 31. https://doi.org/10.1007/s41614-022-00086-0

% For more info check documentation.

  persistent kB;
  persistent elecCharge;
  persistent L;

  % --- performance sensitive calculations ---

  % initialize variables the first time the analyticalHFactorTransport function is called
  if isempty(kB)
    % save local copies of physical constants
    kB = Constant.boltzmann;
    elecCharge = Constant.electronCharge;
    % save local copies of chamber dimensions (with error checking)
    L = chemistry.workCond.chamberLength;
  end

  % --- regular calculations ---
  
  % evaluate electron temperature in Kelvin
  elecTemperature = chemistry.workCond.electronTemperature*elecCharge/kB; % in K

  % evaluate ion mean free path
  sigma = rateCoeffParams{1};           % total ion collision cross section (constant value)
  lambda = 1/(totalGasDensity*sigma); 

  % evaluate ion sound velocity
  ionMass = reaction.reactantArray(1).mass; 
  bohmVel = sqrt(kB*elecTemperature/ionMass);
  
  % evaluate axial h-factor (axial edge-to-center positive ion density ratio)
  hLoL = (pi/2-1)/(L*(1+2/pi^2*(pi/2-1)*L/lambda));
  
  % evaluate de rate coefficient (loss frequency)
  rateCoeff = 2*bohmVel*hLoL;

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', true);
end
