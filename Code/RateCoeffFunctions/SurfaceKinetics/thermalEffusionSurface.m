function [rateCoeff, dependent] = thermalEffusionSurface(~, ~, ~, reaction, rateCoeffParams, chemistry)
% thermalEffusionSurface evaluates a reaction rate coefficient using the expression:
%  C = k_r * exp(E_r/(RT)) (1/tau_r) / {([S]_surface+[F]_surface)*(A/V)}
% The rate coefficient is modified to consider the collision of the volume species with ONE surface species 
%  (dividing by surfaceSiteDensity*AreaOverVolume).
% The current version of this function adopts the wall recombination model, with 
%  tau_r = tau_wall = (4/<v>) (V/A)
% Therefore
%  C = k_r * exp(E_r/(RT)) (<v>/4) / ([S]_surface+[F]_surface)
%  k_r ... stericFactor
%  E_r ... activationEnergy
%  T ... temperature
% 
% Guerra V, Tejero-del-Caz A, Pintassilgo C D and Alves L L 2019 Plasma Sources Sci. Technol. 28 073001
% Guerra V 2007 IEEE Trans. Plasma Sci. 35 1397
% Marinov D, Teixeira C and Guerra V 2017Plasma Processes Polym. 14 1600175

% For more info check documentation.

  persistent dependentInfo;

  % local save of the ID of the volume reactant species
  for i = 1:length(reaction.reactantArray)
    if reaction.reactantArray(i).isVolumeSpecies
      volumeReactantID = reaction.reactantArray(i).ID;
    end
  end

  % --- performance sensitive calculations ---

  % calculations performed once per simulation
  if isempty(dependentInfo)
    % define dependencies of the rate coefficient
      dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  end

  % --- regular calculations ---

  % evaluate steric factor and activation energy
  stericFactor = rateCoeffParams{1};
  activationEnergy = rateCoeffParams{2};

  temperatureStr = rateCoeffParams{3};
  if ~any(strcmp(temperatureStr, {'gasTemperature' 'nearWallTemperature' 'wallTemperature'}))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' is not a valid temperature ' ...
      'to evaluate the thermal velocity with.\nChoose one of the following: ''gasTemperature'', ''nearWallTemperature'' ' ...
      'or ''wallTemperature''.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{3});
  elseif isempty(chemistry.workCond.(temperatureStr))
    error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n''%s'' not found in the working ' ...
      'conditions object.\nPlease, fix the problem and run the code again'], reaction.type, ...
      reaction.description, rateCoeffParams{3});
  else
    T = chemistry.workCond.(temperatureStr);
  end

  % evaluate thermal velocity
  thermalVelocity = sqrt(8*Constant.boltzmann*T/(pi*chemistry.stateArray(volumeReactantID).gas.mass));

  % evaluate rate coefficient
  rateCoeff = stericFactor*exp(-activationEnergy/(Constant.idealGas*T))*thermalVelocity/4/chemistry.workCond.surfaceSiteDensity;

  % set function dependencies
  dependent = dependentInfo;

end
