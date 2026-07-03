function [rateCoeff, dependent] = multicomponentDiffusion(~, densitiesAll, ~, reaction, ~, chemistry)
% multicomponentDiffusion evaluates the diffusion rate coefficient of a particular neutral species interacting with the wall, 
%  adopting multicomponent diffusion (using Wilke's formula) and neglecting the wall reaction time,
%  assuming a density profile that is zero at the walls and equal probability for each wall-recombination channel
%
% Alves L L and Tejero-del-Caz A 2023 Plasma Sources Sci. Technol. 32 054003
% Chantry P J 1987 Journal of Applied Physics, 62(4) 1141-1148. https://doi.org/10.1063/1.339662
% Wilke C R 1950 Chem. Eng. Prog. 46 95-104

% For more info check documentation.

  persistent neutralSpeciesIDs;
  persistent reducedMass;
  persistent sigma;
  persistent epsilon;
  persistent numberOfChannelsArray;
  persistent charDiffLengthSquared;
  persistent dependentInfo;

  % local save of the ID of the reactant species
  reactantID = reaction.reactantArray(1).ID;

  % --- performance sensitive calculations ---

  % calculations performed once per simulation
  if isempty(dependentInfo)
    % evaluate the first zero of the zero order bessel function
    firstBesselZero = fzero(@(x) besselj(0,x), [2.4 2.5]);
    % evaluate geometrical parameters for a cylinder (including infinitely long cylinder and slab limiting cases)
    if chemistry.workCond.chamberLength == 0        % infinitely long cylinder
      charDiffLengthSquared = (chemistry.workCond.chamberRadius/firstBesselZero)^2;
    elseif chemistry.workCond.chamberRadius == 0    % infinitely wide cylinder (slab)
      charDiffLengthSquared = (chemistry.workCond.chamberLength/pi)^2;
    else                                  % finite cylinder
      charDiffLengthSquared = 1/((firstBesselZero/chemistry.workCond.chamberRadius)^2+(pi/chemistry.workCond.chamberLength)^2);
    end
    % define dependencies of the rate coefficient
    dependentInfo = struct('onTime', false, 'onDensities', true, 'onGasTemperature', true, 'onElectronKinetics', false);
    % find IDs of volume neutral species that are going to be taken into account for the multicomponent transport
    % (childless neutral volume species)
    neutralSpeciesIDs = [];
    for i = 1:length(chemistry.stateArray)
      if isempty(chemistry.stateArray(i).ionCharg) && isempty(chemistry.stateArray(i).childArray) && chemistry.stateArray(i).isVolumeSpecies
        neutralSpeciesIDs(end+1) = i;
      end
    end
    % evaluate parameters needed for the calculation of the inverse of the binary diffusion coefficients
    reducedMass = zeros(length(chemistry.stateArray));
    sigma = zeros(length(chemistry.stateArray));
    epsilon = zeros(length(chemistry.stateArray));
    for idx = 1:length(neutralSpeciesIDs)
      for jdx = idx+1:length(neutralSpeciesIDs)
        i = neutralSpeciesIDs(idx);
        j = neutralSpeciesIDs(jdx);
        % error checking
        if isempty(chemistry.stateArray(i).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.mass)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\nMass of %s not found.\n' ...
            'Please, fix the problem and run the code again'], reaction.type, ...
            reaction.description, chemistry.stateArray(j).gas.name);
        elseif isempty(chemistry.stateArray(i).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.lennardJonesDistance)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDistance'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(j).gas.name);
        elseif isempty(chemistry.stateArray(i).gas.lennardJonesDepth)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDepth'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(i).gas.name);
        elseif isempty(chemistry.stateArray(j).gas.lennardJonesDepth)
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            '''lennardJonesDepth'' property of %s not found.\nPlease, fix the problem and run the code again'], ...
            reaction.type, reaction.description, chemistry.stateArray(j).gas.name);
        end
        % evaluation of different auxiliary variables needed
        reducedMass(i,j) = chemistry.stateArray(i).gas.mass*chemistry.stateArray(j).gas.mass/(chemistry.stateArray(i).gas.mass+chemistry.stateArray(j).gas.mass);
        sigma(i,j) = 0.5*(chemistry.stateArray(i).gas.lennardJonesDistance+chemistry.stateArray(j).gas.lennardJonesDistance);
        epsilon(i,j) = sqrt(chemistry.stateArray(i).gas.lennardJonesDepth*chemistry.stateArray(j).gas.lennardJonesDepth);
        reducedMass(j,i) = reducedMass(i,j);
        sigma(j,i) = sigma(i,j);
        epsilon(j,i) = epsilon(i,j);
      end
    end
    % initialize array of numberOfChannels (number of destruction channels at the wall for each species)
    numberOfChannelsArray = zeros(size(chemistry.stateArray));
  end

  % evaluate number of destruction channels at the wall for the diffused species (only done once per lost species)
  if numberOfChannelsArray(reactantID) == 0
    % evaluate total wall reaction coefficient
    for rxn = chemistry.stateArray(reactantID).reactionsDestruction
      if rxn.isTransport
        if strcmp(rxn.type, reaction.type)
          numberOfChannelsArray(reactantID) = numberOfChannelsArray(reactantID) + 1;
        else % error checking
          error(['Error found when evaluating ''%s'' rate coefficient for reaction:\n%s.\n' ...
            'It has been found another transport reaction of a different type for the same species:\n%s.\nPlease, ' ...
            'fix the problem and run the code again'], reaction.type, ...
            reaction.description, rxn.description);
        end
      end
    end
  end

  % --- regular calculations ---

  % evaluate total neutral density (neutralSpeciesIDs only include childless neutral species)
  neutralDens = sum(densitiesAll(neutralSpeciesIDs));

  % evaluation of the inverse of the reduced binary diffusion coeficients
  inverseReducedBinaryDiffCoeff = zeros(size(densitiesAll));
  for i = neutralSpeciesIDs
    if i == reactantID
      continue;
    end
    Taux = chemistry.workCond.gasTemperature*Constant.boltzmann/epsilon(i,reactantID);
    % fit to first order collision integral from (R.J. Kee, M.E. Coltrin, P. Glarborg "Chemically Reacting Flow:
    % Theory and Practice" John Wiley (2003), Pag. 492)
    collIntegral = 1.0548*Taux^-0.15504+(Taux+0.55909)^-2.1705;
    % actual evaluation of the inverse of the reduced binary diffusion coeficient between species i and j
    inverseReducedBinaryDiffCoeff(i) = 16/3*sqrt(pi*reducedMass(i,reactantID)/(2*Taux*epsilon(i,reactantID)))*...
      sigma(i,reactantID)^2*collIntegral;
  end

  % evaluate effective diffusion coefficient (using Wilke's formula)
  effDiffCoeff = (1-densitiesAll(reactantID)/neutralDens)/sum(densitiesAll.*inverseReducedBinaryDiffCoeff);

  % evaluate rate coefficient
  rateCoeff = (1/numberOfChannelsArray(reactantID)) * (effDiffCoeff/charDiffLengthSquared);

  % set function dependencies
  dependent = dependentInfo;

end
