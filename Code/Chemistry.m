% LoKI-GM comprises two modules, that can run self-consistently coupled 
% or as standalone tools.
% 
% LoKI-B, which solves the space independent form of the two-term 
% electron Boltzmann equation (EBE) to calculate the isotropic and the 
% anisotropic parts of the electron distribution function, 
% and the associated electron macroscopic parameters. 
% LoKI-B applies to non-magnetised non-equilibrium LTPs, excited by 
% DC/HF electric fields or time-dependent (non-oscillatory) electric fields 
% from different gases or gas mixtures. 
% The tool uses a stationary description for DC fields, 
% a Fourier time-expansion description for HF fields, 
% and a time-dependent description for time-varying fields.
% 
% LoKI-C, which solves the system of zero-dimensional (volume average) 
% rate balance equations for the most relevant 
% charged and neutral species in the plasma. 
% LoKI-C receives as input data the kinetic schemes for the gas/plasma/
% surface system under study, via an intuitive csv-like input file, 
% and gives as output the particle densities of the different gas/plasma/
% surface species, the corresponding creation/destruction reaction rates, 
% and the reduced electric field (and any related quantity, such as 
% the discharge current or the discharge power-density).
% The tool uses several modules to describe the mechanisms 
% (collisional, radiative and transport) controlling the
% creation/destruction of species, namely various transport models 
% for the charged particles and for the neutral particles. 
% LoKI-C includes also a gas/plasma thermal model, for the self-consistent 
% calculation of the gas temperature, and supports multicomponent 
% mean-field microkinetic mesoscopic models to handle surface kinetics 
% in a fully coupled way with volume kinetics.
%
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% Copyright (C) 2026 L. L. Alves, A. Tejero-del-Caz, T. C. Dias, 
% A. Gonçalves, L. Marques, P. Pereira, N. Pinhão, C. D. Pintassilgo, 
% T. Silva, P. Viegas and V. Guerra
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

classdef Chemistry < handle
  %Chemistry Class that solves the heavy species rate balance equations under certain conditions to
  %time evoluation of their densities.
  %
  %   Chemistry uses an stiff ODE solver (ode15s Matlab solver) to solve the rate balance equations corresponding
  %   to the set of reactions specified in the ".chem" files provided by the user in the setup file.


  properties
    
    convergenceParameter = '';

    gasArray = ChemGas.empty;           % handle array to the gas mixture (to perform operations in a per gas basis)
    stateArray = ChemState.empty;       % handle array to the states (species) taken into account
    reactionArray = Reaction.empty;     % handle array to the reactions taken into account
    numberOfSpecies;
    initialDensities;
    gasIDs;
    childIDs;
    volumePhaseSpeciesIDs;
    childlessVolumePhaseSpeciesIDs;
    surfacePhaseSpeciesIDs;
    numberOfReactions;
    rateCoeffFuncHandles;
    rateCoeffParams;
    reactantElectrons;
    reactantIDs;
    reactantStoiCoeffs;
    catalystIDs;
    catalystStoiCoeffs;
    productElectrons;
    productIDs;
    productStoiCoeffs;
    gasStabilisedReactionIDs;
    transportReactionIDs;
    volumeReactionIDs;
    surfaceReactionIDs;

    workCond = WorkingConditions.empty; % handle to the working conditions of the simulation
    ensureIsobaric = false;             % boolean to select whether the 'ensureIsobaric' flow model is active or not
    initialGasPressure = [];            % initial value of gas pressure for the integration of the rate balance equations
    targetGasPressure = [];             % target value of gas pressure precribed in the working conditions
    
    isPulsed = false;                   % boolean to select whether the simulation is pulsed or discharge + post-discharge
    pulseFunction;                      % function handle of the pulse function (if pulsed simulation)
    pulseFirstStep;                     % first time step of the pulse (if pulsed simulation)
    pulseFinalTime;                     % final time of the pulse (if pulsed simulation)
    pulseSamplingType;                  % type of time-sampling: logscale or linspace (if pulsed simulation)
    pulseSamplingPoints;                % number of time points for the sampling of the pulse (if pulsed simulation)
    pulseFunctionParameters;            % parameters of the pulse function (if pulsed simulation)
    lookUpMethod;                       % method to interpolate from look up tables (if pulsed simulation)
    lookUpTablePower;                   
    lookUpTableSwarm;
    lookUpTableRateCoeff;
    lookUpTableRedFieldValues;
    lookUpTableEleTempValues;
    lookUpTableRedMobValues;
    lookUpTableRedDiffValues;
    lookUpTableDriftVelocityValues;
    lookUpTablePowerFieldValues;
    lookUpTablePowerElasticNetValues;
    lookUpTablePowerVibrationalNetValues;
    lookUpTablePowerRotationalNetValues;
    lookUpTablePowerCARNetValues;

    electronKinetics = [];              % handle to the electron kinetics object of the simulation
    gasesIDsToUpdateInElectronKinetics = [];
    statesIDsToUpdateInElectronKinetics = [];
    electronTransportProperties = struct('reducedDiffCoeff', [], 'reducedMobility', []);

    includeThermalModel = false;        % boolean to select whether the thermal model is active or not
    thermalModelBoundary = '';          % wall or external (location of the model boundary condition)
    thermalModelWallFraction = 0;       % fraction of wall-released power heating the plasma (see fw in documentation)
    intConvCoeff = [];                  % internal convection coefficient (of the cylindrical wall)
    extConvCoeff = [];                  % external convection coefficient (of the cylindrical wall)
    initialGasTemperature = [];         % initial value of gas temperature for the integration of the rate balance equations
    initialNearWallTemperature = [];
    initialWallTemperature = [];
    includeJouleHeating = false;        % boolean to select if Joule heating is used as source term in the thermal model
    electronKineticsGasTemperature;     % gas temperature used in the electron kinetics module (updated in the global cycle)
    includeGasRotHeating;
    includeGasVibHeating;

    solveEedf = false;                          % boolean to select whether to solve the EEDF when solving the electron density
    imposeQuasiNeutrality = false;              % boolean to select whether to impose quasi-neutrality should be imposed during temporal integration
    electronKineticsDependence = 'steadyState'; % string to select the dependence of the electron kinetics during temporal integration: 'steadyState', 'postDischarge' or 'quasiStationary' 

    odeSolver;
    odeOptions;
    odeDischargeTime = [];
    odePostDischargeTime = [];

    pressureRelError = [];
    pressureMaxIterations = [];
    pressureIterationCurrent = [];
    pressureRelErrorCurrent = [];
    neutralityRelError = [];
    neutralityMaxIterations = [];
    neutralityIterationCurrent = [];
    neutralityRelErrorCurrent = [];
    globalRelError = [];
    globalMaxIterations = [];
    globalIterationCurrent = [];
    globalRelErrorCurrent = [];

    solution;

  end

  events
    genericStatusMessage;
    newPressureCycleIteration;
    newNeutralityCycleIteration;
    newGlobalCycleIteration;
    obtainedNewChemistrySolution;
  end

  methods (Access = public)

    function chemistry = Chemistry(setup)

      % store variable which defines the type of neutrality cycle
      chemistry.convergenceParameter = setup.info.chemistry.convergenceParameter;

      % store the gas array (in order to be able to perform operations in a per gas basis)
      chemistry.gasArray = setup.chemistryGasArray;
      % find and store IDs of gases with electron kinetics equivalents that needs to be updated
      gasesIDsToUpdateInElectronKinetics = [];
      for i = 1:length(chemistry.gasArray)
        gas = chemistry.gasArray(i);
        if ~isempty(gas.eedfEquivalent) && ~isempty(gas.eedfEquivalent.collisionArray)
          gasesIDsToUpdateInElectronKinetics(end+1) = i;
        end
      end
      chemistry.gasesIDsToUpdateInElectronKinetics = gasesIDsToUpdateInElectronKinetics;

      % store the state array (with all the species considered in the chemistry's gas mixture)
      chemistry.stateArray = setup.chemistryStateArray;

      % store the reaction array (with all the reactions considered in the chemistry)
      chemistry.reactionArray = setup.chemistryReactionArray;

      % convert species info into simple and fast numeric/cell arrays (ChemState objects -> arrays)
      numberOfSpecies = length(chemistry.stateArray);
      initialDensities = zeros(1, numberOfSpecies);
      childIDs = cell(1, numberOfSpecies);
      gasIDs = zeros(1, numberOfSpecies);
      volumePhaseSpeciesIDs = [];
      childlessVolumePhaseSpeciesIDs = [];
      surfacePhaseSpeciesIDs = [];
      for i = 1:numberOfSpecies
        initialDensities(i) = chemistry.stateArray(i).density;
        gasIDs(i) = chemistry.stateArray(i).gas.ID;
        if ~isempty(chemistry.stateArray(i).childArray)
          childIDs{i} = [chemistry.stateArray(i).childArray.ID];
        end
        if chemistry.stateArray(i).isVolumeSpecies
          volumePhaseSpeciesIDs(end+1) = i;
          if isempty(chemistry.stateArray(i).childArray)
            childlessVolumePhaseSpeciesIDs(end+1) = i;
          end
        else
          surfacePhaseSpeciesIDs(end+1) = i;
        end
      end
      chemistry.numberOfSpecies = numberOfSpecies;
      chemistry.initialDensities = initialDensities;
      chemistry.gasIDs = gasIDs;
      chemistry.childIDs = childIDs;
      chemistry.volumePhaseSpeciesIDs = volumePhaseSpeciesIDs;
      chemistry.childlessVolumePhaseSpeciesIDs = childlessVolumePhaseSpeciesIDs;
      chemistry.surfacePhaseSpeciesIDs = surfacePhaseSpeciesIDs;

      % convert reactions info into simple and fast numeric/cell arrays (Reaction objects -> arrays)
      numberOfReactions = length(chemistry.reactionArray);
      rateCoeffFuncHandles = cell(1, numberOfReactions);
      rateCoeffParams = cell(1, numberOfReactions);
      reactantElectrons = zeros(1, numberOfReactions);
      reactantIDs = cell(1, numberOfReactions);
      reactantStoiCoeffs = zeros(numberOfSpecies, numberOfReactions);
      catalystIDs = cell(1, numberOfReactions);
      catalystStoiCoeffs = zeros(numberOfSpecies, numberOfReactions);
      productElectrons = zeros(1, numberOfReactions);
      productIDs = cell(1, numberOfReactions);
      productStoiCoeffs = zeros(numberOfSpecies, numberOfReactions);
      gasStabilisedReactionIDs = [];
      transportReactionIDs = [];
      volumeReactionIDs = [];
      surfaceReactionIDs = [];
      for i = 1:numberOfReactions
        rateCoeffFuncHandles{i} = chemistry.reactionArray(i).rateCoeffFuncHandle;
        rateCoeffParams{i} = chemistry.reactionArray(i).rateCoeffParams;
        reactantElectrons(i) = chemistry.reactionArray(i).reactantElectrons;
        reactantIDs{i} = [chemistry.reactionArray(i).reactantArray.ID];
        reactantStoiCoeffs(reactantIDs{i}, i) = chemistry.reactionArray(i).reactantStoiCoeff;
        catalystIDs{i} = [chemistry.reactionArray(i).catalystArray.ID];
        catalystStoiCoeffs(catalystIDs{i}, i) = chemistry.reactionArray(i).catalystStoiCoeff;
        productElectrons(i) = chemistry.reactionArray(i).productElectrons;
        productIDs{i} = [chemistry.reactionArray(i).productArray.ID];
        productStoiCoeffs([productIDs{i}], i) = chemistry.reactionArray(i).productStoiCoeff;
        if chemistry.reactionArray(i).isGasStabilised
          gasStabilisedReactionIDs(end+1) = i;
        end
        isSurface = false;
        for ID = [reactantIDs{i} productIDs{i}]
          if chemistry.stateArray(ID).isSurfaceSpecies
            isSurface = true;
            break;
          end
        end
        if isSurface
          surfaceReactionIDs(end+1) = i;
        elseif chemistry.reactionArray(i).isTransport
          transportReactionIDs(end+1) = i;
        else
          volumeReactionIDs(end+1) = i;
        end
      end
      chemistry.numberOfReactions = numberOfReactions;
      chemistry.rateCoeffFuncHandles = rateCoeffFuncHandles;
      chemistry.rateCoeffParams = rateCoeffParams;
      chemistry.reactantElectrons = reactantElectrons;
      chemistry.reactantIDs = reactantIDs;
      chemistry.reactantStoiCoeffs = sparse(reactantStoiCoeffs);
      chemistry.catalystIDs = catalystIDs;
      chemistry.catalystStoiCoeffs = sparse(catalystStoiCoeffs);
      chemistry.productElectrons = productElectrons;
      chemistry.productIDs = productIDs;
      chemistry.productStoiCoeffs = sparse(productStoiCoeffs);
      chemistry.gasStabilisedReactionIDs = gasStabilisedReactionIDs;
      chemistry.transportReactionIDs = transportReactionIDs;
      chemistry.volumeReactionIDs = volumeReactionIDs;
      chemistry.surfaceReactionIDs = surfaceReactionIDs;

      % store working conditions and add corresponding listeners (if needed)
      chemistry.workCond = setup.workCond;

      % store configuration for the 'ensureIsobaric' flow model
      chemistry.ensureIsobaric = strcmp(setup.info.workingConditions.totalSccmOutFlow,'ensureIsobaric');

      % store the initial value of the gas pressure for the temporal integration of the rate balance equations
      % (this initial value is updated in the iterations of the pressure cycle)
      chemistry.initialGasPressure = setup.workCond.gasPressure;

      % store the target value of the gas pressure precribed in the working conditions
      chemistry.targetGasPressure = setup.workCond.gasPressure;

      % check if the the simulation is pulsed (oposed to discharge + post-discharge)
      if setup.pulsedSimulation
        % store information about pulsed simulation in case it is activated
        chemistry.isPulsed = true;
        chemistry.pulseFunction = setup.pulseInfo.function;
        chemistry.pulseFirstStep = setup.pulseInfo.firstStep;
        chemistry.pulseFinalTime = setup.pulseInfo.finalTime;
        chemistry.pulseSamplingType = setup.pulseInfo.samplingType;
        chemistry.pulseSamplingPoints = setup.pulseInfo.samplingPoints;
        chemistry.pulseFunctionParameters = setup.pulseInfo.functionParameters;
        chemistry.pulseFunctionParameters{end+1} = chemistry;
        
        % save look up tables (THIS SHOULD BE MOVED TO PARSE CLASS AND BE DONE AT SETUP TIME)
        chemistry.lookUpMethod = setup.info.chemistry.lookUpTables.lookUpMethod;
        chemistry.lookUpTablePower = readtable([Parse.inputFolder filesep setup.info.chemistry.lookUpTables.power], ...
          ReadVariableNames=false, VariableDescriptionsLine=1);
        chemistry.lookUpTableSwarm = readtable([Parse.inputFolder filesep setup.info.chemistry.lookUpTables.swarm], ... 
          ReadVariableNames=false, VariableDescriptionsLine=1);
        % Extract the original unsanitized headers from the file
        lines = readlines([Parse.inputFolder filesep setup.info.chemistry.lookUpTables.rateCoeff]);
        reactionsDescriptions = {};
        for i=1:length(lines)
          aux = split(lines(i));
          if length(aux)>1 && ~isnan(str2double(aux(2)))
            reactionDescriptionParts = split(aux(3),',');
            reactionDescription = reactionDescriptionParts(1);
            for j = 2:length(reactionDescriptionParts)-1
              reactionDescription = sprintf('%s,%s',reactionDescription, reactionDescriptionParts(j));
            end
            reactionsDescriptions{str2double(aux(2))} = char(reactionDescription);
          elseif strncmp(lines(i), 'RedField(Td)', 12)
            startLineIdx = i;
            break;
          else
            continue;
          end
        end
        chemistry.lookUpTableRateCoeff = readtable([Parse.inputFolder filesep setup.info.chemistry.lookUpTables.rateCoeff], ...
          ReadVariableNames=false, CommentStyle='#', VariableDescriptionsLine=startLineIdx);
          variableDescriptions = chemistry.lookUpTableRateCoeff.Properties.VariableDescriptions;
        for i=2:length(variableDescriptions)
          description = variableDescriptions(i);
          tokens = regexp(description,'R(\d+)_(ine|sup)','tokens');
          idx = str2double(tokens{1}{1}{1});
          direction = tokens{1}{1}{2};
          variableDescriptions{i} = sprintf('%s|%s',reactionsDescriptions{idx},direction);
        end
        chemistry.lookUpTableRateCoeff.Properties.VariableDescriptions = variableDescriptions;

        % check consistency of the look-up tables (number and values of reduced field should be the same)
        if length(chemistry.lookUpTablePower.Var1) ~= length(chemistry.lookUpTableSwarm.Var1) || ...
          length(chemistry.lookUpTablePower.Var1) ~= length(chemistry.lookUpTableRateCoeff.Var1)
          error("Inconsistency detected in the provided look-up tables." + newline + ...
            "The number of reduced field values in the provided tables are different." + newline + ...
            "Please check the look-up tables and run the code again.");
        elseif any(chemistry.lookUpTablePower.Var1~=chemistry.lookUpTableSwarm.Var1) || ...
          any(chemistry.lookUpTablePower.Var1~=chemistry.lookUpTableRateCoeff.Var1)
          error("Inconsistency detected in the provided look-up tables." + newline + ...
            "The reduced field values in the provided tables are different." + newline + ...
            "Please check the provided look-up tables and run the code again.")
        end

        % extract frequently used columns from the look up tables for easier access
        % load reduced field values
        for idx = 1:length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTableSwarm.Properties.VariableDescriptions{idx}, "RedField(Td)")
            chemistry.lookUpTableRedFieldValues = chemistry.lookUpTableSwarm{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
              error('Reduced field values not found in the provided look-up table for the swarm parameters.');
          end
        end
        % load electron temperature values
        for idx = 1:length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTableSwarm.Properties.VariableDescriptions{idx}, "EleTemp(eV)")
            chemistry.lookUpTableEleTempValues = chemistry.lookUpTableSwarm{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
              error('Electron temperature values not found in the provided look-up table for the swarm parameters.');
          end
        end
        % load reduced mobility values
        for idx = 1:length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTableSwarm.Properties.VariableDescriptions{idx}, "RedMob((msV)^-1)")
            chemistry.lookUpTableRedMobValues = chemistry.lookUpTableSwarm{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
              error('Reduced mobility values not found in the provided look-up table for the swarm parameters.');
          end
        end
        % load reduced diffusion values
        for idx = 1:length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTableSwarm.Properties.VariableDescriptions{idx}, "RedDiff((ms)^-1)")
            chemistry.lookUpTableRedDiffValues = chemistry.lookUpTableSwarm{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
              error('Reduced diffusion values not found in the provided look-up table for the swarm parameters.');
          end
        end
        % load drift velocity values
        for idx = 1:length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTableSwarm.Properties.VariableDescriptions{idx}, "DriftVelocity(ms^-1)")
            chemistry.lookUpTableDriftVelocityValues = chemistry.lookUpTableSwarm{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTableSwarm.Properties.VariableDescriptions)
              error('Drift velocity values not found in the provided look-up table for the swarm parameters.');
          end
        end
        % load power from field values
        for idx = 1:length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTablePower.Properties.VariableDescriptions{idx}, "PowerField(eVm^3s^-1)")
            chemistry.lookUpTablePowerFieldValues = chemistry.lookUpTablePower{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
              error('Power from field values not found in the provided look-up table for the power.');
          end
        end
        % load net elastic power values
        for idx = 1:length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTablePower.Properties.VariableDescriptions{idx}, "PwrElaNet(eVm^3s^-1)")
            chemistry.lookUpTablePowerElasticNetValues = chemistry.lookUpTablePower{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
              error('Net elastic power values not found in the provided look-up table for the power.');
          end
        end
        % load net CAR power values
        for idx = 1:length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTablePower.Properties.VariableDescriptions{idx}, "PwrCARNet(eVm^3s^-1)")
            chemistry.lookUpTablePowerCARNetValues = chemistry.lookUpTablePower{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
              error('Net CAR power values not found in the provided look-up table for the power.');
          end
        end
        % load net vibrational power values
        for idx = 1:length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTablePower.Properties.VariableDescriptions{idx}, "PwrVibNet(eVm^3s^-1)")
            chemistry.lookUpTablePowerVibrationalNetValues = chemistry.lookUpTablePower{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
              error('Net vibrational power values not found in the provided look-up table for the power.');
          end
        end
        % load net rotational power values
        for idx = 1:length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
          if strcmp(chemistry.lookUpTablePower.Properties.VariableDescriptions{idx}, "PwrRotNet(eVm^3s^-1)")
            chemistry.lookUpTablePowerRotationalNetValues = chemistry.lookUpTablePower{:,idx};
            break;
          elseif idx == length(chemistry.lookUpTablePower.Properties.VariableDescriptions)
              error('Net rotational power values not found in the provided look-up table for the power.');
          end
        end

        % set initial value of the reduced electric field in the working conditions object
        chemistry.workCond.update('reducedField', chemistry.pulseFunction(0, chemistry.pulseFunctionParameters));

        % set initial value of electron macroscopic properties (interpolation from look-up tables)
        switch chemistry.lookUpMethod
          case 'localField'
            DN = interp1(chemistry.lookUpTableRedFieldValues, chemistry.lookUpTableRedDiffValues, ...
              chemistry.workCond.reducedField);
            muN = interp1(chemistry.lookUpTableRedFieldValues, chemistry.lookUpTableRedMobValues, ...
              chemistry.workCond.reducedField);
            electronTemperature = interp1(chemistry.lookUpTableRedFieldValues, chemistry.lookUpTableEleTempValues, ...
              chemistry.workCond.reducedField);
            if isnan(DN) || isnan(muN) || isnan(electronTemperature)
              error(['Electron macroscopic properties out of look-up table range.\n'
                'Please check the look-up table and the range of the reduced field during the simulation.'])
            end
            chemistry.workCond.update('electronTemperature', electronTemperature);
          case 'localEnergy'
            DN = interp1(chemistry.lookUpTableEleTempValues, chemistry.lookUpTableRedDiffValues, ... 
              chemistry.workCond.electronTemperature);
            muN = interp1(chemistry.lookUpTableEleTempValues, chemistry.lookUpTableRedMobValues, ... 
              chemistry.workCond.electronTemperature);
            if isnan(DN) || isnan(muN)
              error(['Reduced transport parameters out of look-up table range.\n'
                'Please check the look-up table and the range of the reduced field during the simulation.'])
            end
        end
        chemistry.electronTransportProperties.reducedDiffCoeff = DN;
        chemistry.electronTransportProperties.reducedMobility = muN;
        % in pulsed simulations electron kinetics depends on look up tables (quasi-stationary approximation)
        chemistry.electronKineticsDependence = 'quasiStationary';
        % in pulsed simulations quasi-neutrality is imposed during the integration of the rate balance equations
        chemistry.imposeQuasiNeutrality = true;
      end

      % store electron kinetics object in case it is enabled
      if setup.enableElectronKinetics
        chemistry.electronKinetics = setup.electronKinetics;
        % find and store IDs of states with electron kinetics equivalents that needs to be updated
        statesIDsToUpdateInElectronKinetics = [];
        for i = 1:length(chemistry.stateArray)
          state = chemistry.stateArray(i);
          if ~isempty(state.eedfEquivalent) && state.eedfEquivalent.isTarget
            statesIDsToUpdateInElectronKinetics(end+1) = i;
          end
        end
        chemistry.statesIDsToUpdateInElectronKinetics = statesIDsToUpdateInElectronKinetics;
      end

      % store electron transport properties (in case they are specified in the setup file) (TO BE REMOVED ONCE PULSED SIMULATIONS ARE IMPLEMENTED??)
      if isfield(setup.info.chemistry, 'electronProperties')
        chemistry.electronTransportProperties.reducedDiffCoeff = setup.info.chemistry.electronProperties.reducedDiffCoeff;
        chemistry.electronTransportProperties.reducedMobility = setup.info.chemistry.electronProperties.reducedMobility;
      end

      % store configuration of the thermal model
      chemistry.includeThermalModel = setup.info.chemistry.thermalModel.isOn;
      chemistry.initialGasTemperature = setup.workCond.gasTemperature;
      chemistry.electronKineticsGasTemperature = setup.workCond.gasTemperature;
      if chemistry.includeThermalModel
        % save boundary condition of the thermal model
        chemistry.thermalModelBoundary = setup.info.chemistry.thermalModel.boundary;
        % save machroscopic parameters of the thermal model (according to boundary condition)
        switch chemistry.thermalModelBoundary
          case 'wall'
            chemistry.intConvCoeff = setup.info.chemistry.thermalModel.intConvCoeff;
            chemistry.thermalModelWallFraction = setup.info.chemistry.thermalModel.wallFraction;
          case 'external'
            chemistry.intConvCoeff = setup.info.chemistry.thermalModel.intConvCoeff;
            chemistry.extConvCoeff = setup.info.chemistry.thermalModel.extConvCoeff;
        end
        % save/set initial temperatures
        if isempty(setup.workCond.nearWallTemperature)
          setup.workCond.nearWallTemperature = setup.workCond.gasTemperature;
        end
        if strcmp(chemistry.thermalModelBoundary, 'external') && isempty(setup.workCond.wallTemperature)
          setup.workCond.wallTemperature = setup.workCond.extTemperature;
        end
        chemistry.initialNearWallTemperature = setup.workCond.nearWallTemperature;
        chemistry.initialWallTemperature = setup.workCond.extTemperature;
        % check if Joule heating is or isn't used as source term in the thermal model (false by default)
        if isfield(setup.info.chemistry.thermalModel,'includeJouleHeating')
          chemistry.includeJouleHeating = setup.info.chemistry.thermalModel.includeJouleHeating;
        end
        if ~chemistry.includeJouleHeating
          % check for which gases the rotational and vibrational heating should be considered from the electron kinetics
          nGases = length(chemistry.gasArray);
          chemistry.includeGasRotHeating = false(1,nGases);
          chemistry.includeGasVibHeating = false(1,nGases);
          for i = 1:nGases
            gas = chemistry.gasArray(i);
            eedfEquivalent = gas.eedfEquivalent;
            if ~isempty(eedfEquivalent)
              % include rotational heating from rotational states defined in electron kinetics but not in chemistry
              chemistry.includeGasRotHeating(i) = any(strcmp({eedfEquivalent.stateArray.type},'rot')) && ...
                ~any(strcmp({gas.stateArray.type},'rot'));
              % include rotational heating from vibrational states defined in electron kinetics but not in chemistry
              chemistry.includeGasVibHeating(i) = any(strcmp({eedfEquivalent.stateArray.type},'vib')) && ...
                ~any(strcmp({gas.stateArray.type},'vib'));
            end
          end
        end
      end

      % store configuration of the ODE solver
      chemistry.odeSolver = str2func(setup.info.chemistry.timeIntegrationConf.odeSolver);
      options = odeset();
      for parameter = fields(options)'
        if isfield(setup.info.chemistry.timeIntegrationConf, 'odeSetParameters') && ...
          isfield(setup.info.chemistry.timeIntegrationConf.odeSetParameters, parameter{1})
          options.(parameter{1}) = setup.info.chemistry.timeIntegrationConf.odeSetParameters.(parameter{1});
        else
          options.(parameter{1}) = [];
        end
      end
      options.NonNegative = 1:numberOfSpecies+1; % ensure non negative values for the solution (densities & temperature)
      %options.OutputFcn = @odeProgressBar;      % activate only for debugging purposes
      chemistry.odeOptions = options;
      % set a default value for the relative tolerance of the ODE solver if it is not specified in the setup file
      if isempty(chemistry.odeOptions.RelTol)
        chemistry.odeOptions.RelTol = 1e-3;
      end

      % store configuration about the iterations schemes (pressure, quasineutrality and global) if not pulsed simulation
      if ~chemistry.isPulsed
        chemistry.odeDischargeTime = setup.info.chemistry.timeIntegrationConf.dischargeTime;
        chemistry.odePostDischargeTime = setup.info.chemistry.timeIntegrationConf.postDischargeTime;
        if isfield(setup.info.chemistry.iterationSchemes, 'pressureRelError')
          if setup.info.chemistry.iterationSchemes.pressureRelError < chemistry.odeOptions.RelTol
            warning('Pressure relative error < Solver relative tolerance')
          end
          chemistry.pressureRelError = setup.info.chemistry.iterationSchemes.pressureRelError;
        else
         chemistry.pressureRelError = 2*chemistry.odeOptions.RelTol;
        end  
        chemistry.pressureMaxIterations = setup.info.chemistry.iterationSchemes.pressureMaxIterations;
        if isempty(chemistry.electronKinetics)
          chemistry.neutralityRelError = [];
          chemistry.neutralityMaxIterations = 1;
          chemistry.globalRelError = [];
          chemistry.globalMaxIterations = 1;
          warning(['electronKinetics field is empty. Therefore, neutrality and global cycles have been deactivated\n%s'], ...
              'Relative errors are set to empty and maximum iterations are set to 1.')
        else
          if isfield(setup.info.chemistry.iterationSchemes, 'neutralityRelError')
            if setup.info.chemistry.iterationSchemes.neutralityRelError < chemistry.odeOptions.RelTol
              warning('Neutrality relative error < Solver relative tolerance')
            end
            chemistry.neutralityRelError = setup.info.chemistry.iterationSchemes.neutralityRelError;
          else
            chemistry.neutralityRelError = 2*chemistry.odeOptions.RelTol;
          end
          chemistry.neutralityMaxIterations = setup.info.chemistry.iterationSchemes.neutralityMaxIterations;
          if isfield(setup.info.chemistry.iterationSchemes, 'globalRelError')
            if setup.info.chemistry.iterationSchemes.globalRelError < chemistry.odeOptions.RelTol
              warning('Global relative error < Solver relative tolerance')
            end
            chemistry.globalRelError = setup.info.chemistry.iterationSchemes.globalRelError;
          else
            chemistry.globalRelError = 2*chemistry.odeOptions.RelTol;
          end
          chemistry.globalMaxIterations = setup.info.chemistry.iterationSchemes.globalMaxIterations;
        end
      else
        chemistry.odeDischargeTime = chemistry.pulseFinalTime;
        chemistry.odePostDischargeTime = 0;
      end
    end

    function solve(chemistry)

      if chemistry.isPulsed
        chemistry.solvePulsed();
      else
        chemistry.solveDischargeAndPostDischarge();
      end

    end
    
  end

  methods (Access = private)

    function solveDischargeAndPostDischarge(chemistry)

      % Store initial values of the transport properties for electrons (when the electron kinetics module is active)
      if ~isempty(chemistry.electronKinetics)
        chemistry.electronTransportProperties.reducedDiffCoeff = chemistry.electronKinetics.swarmParam.redDiffCoeff;
        chemistry.electronTransportProperties.reducedMobility = chemistry.electronKinetics.swarmParam.redMobility;
      end

      % Run Chemistry - Kinetics coupling cycles

      [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
        chemistry.globalCycle();

      % evaluate density of parent states in the final solution (time evolution) and final total gas density (final time)
      totalGasDensity = 0;
      for i = 1:chemistry.numberOfSpecies
        if isempty(chemistry.childIDs{i}) && any(i == chemistry.volumePhaseSpeciesIDs)
          totalGasDensity = totalGasDensity+absDensitiesTime(end,i);
        elseif ~isempty(chemistry.childIDs{i})
          absDensitiesTime(:,i) = 0;
          for j = chemistry.childIDs{i}
            if ~isempty(chemistry.childIDs{j})
              absDensitiesTime(:,j) = 0;
              for k = chemistry.childIDs{j}
                absDensitiesTime(:,j) = absDensitiesTime(:,j) + absDensitiesTime(:,k);
              end
            end
            absDensitiesTime(:,i) = absDensitiesTime(:,i) + absDensitiesTime(:,j);
          end
        end
      end

      % evaluate temporal evolution of intermediate temperatures (if thermal model is activated)
      nearWallTemperatureTime = [];
      wallTemperatureTime = [];
      if chemistry.includeThermalModel
        switch chemistry.thermalModelBoundary
          case 'wall'
            nearWallTemperatureTime = zeros(1,length(gasTemperatureTime));
            nearWallTemperatureTime(1) = chemistry.initialNearWallTemperature;
            for i = 1:length(gasTemperatureTime)-1
              chemistry.workCond.gasTemperature = gasTemperatureTime(i);
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
              [~, thermalModel] = kinetics(time(i), [absDensitiesTime(i,:) gasTemperatureTime(i)]', ...
                directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
              nearWallTemperatureTime(i+1) = thermalModel.nearWallTemperature;
            end
          case 'external'
            nearWallTemperatureTime = zeros(1,length(gasTemperatureTime));
            nearWallTemperatureTime(1) = chemistry.initialNearWallTemperature;
            wallTemperatureTime = zeros(1,length(gasTemperatureTime));
            wallTemperatureTime(1) = chemistry.initialWallTemperature;
            for i = 1:length(gasTemperatureTime)-1
              chemistry.workCond.gasTemperature = gasTemperatureTime(i);
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
              chemistry.workCond.wallTemperature = wallTemperatureTime(i);
              [~, thermalModel] = kinetics(time(i), [absDensitiesTime(i,:) gasTemperatureTime(i)]', ...
                directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
              nearWallTemperatureTime(i+1) = thermalModel.nearWallTemperature;
              wallTemperatureTime(i+1) = thermalModel.wallTemperature;
            end
        end
      else
        thermalModel = struct.empty;
      end


      % evaluate final (last time point) rates

      % evaluate time dependent rate coefficients (at final time)
      chemistry.workCond.gasTemperature = gasTemperatureTime(end);
      if chemistry.includeThermalModel
        switch chemistry.thermalModelBoundary
          case 'wall'
            chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(end);
          case 'external'
            chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(end);
            chemistry.workCond.wallTemperature = wallTemperatureTime(end);
        end
      end

      for i = 1:length(timeDependentReactionIDs)
        reactionID = timeDependentReactionIDs(i);
        directRateCoeffs(reactionID) = chemistry.rateCoeffFuncHandles{reactionID}(time(end), ...
          absDensitiesTime(end,:)', totalGasDensity, chemistry.reactionArray(reactionID), chemistry.rateCoeffParams{reactionID}, chemistry);
        if chemistry.reactionArray(reactionID).isReverse
          inverseRateCoeffs(reactionID) = detailedBalance(chemistry.reactionArray(reactionID), ...
            directRateCoeffs(reactionID), Constant.boltzmannInEV*chemistry.workCond.gasTemperature);
        end
      end

      % evaluate reaction rate
      reactionRates = directRateCoeffs.*prod(repmat(absDensitiesTime(end,:)',[1 chemistry.numberOfReactions]).^ ...
        (chemistry.reactantStoiCoeffs+chemistry.catalystStoiCoeffs),1).* ...
        chemistry.workCond.electronDensity.^chemistry.reactantElectrons- ...
        inverseRateCoeffs.*prod(repmat(absDensitiesTime(end,:)',[1 chemistry.numberOfReactions]).^ ...
        (chemistry.productStoiCoeffs+chemistry.catalystStoiCoeffs),1).* ...
        chemistry.workCond.electronDensity.^chemistry.productElectrons;
      for i = chemistry.gasStabilisedReactionIDs
        reactionRates(i) = reactionRates(i)*totalGasDensity;
      end

      % store final discharge solution of the heavy species kinetics in the chemistry properties
      finalDischargeDensity = absDensitiesTime(end,:);
      finalDischargeDensity(chemistry.surfacePhaseSpeciesIDs) = ...
        finalDischargeDensity(chemistry.surfacePhaseSpeciesIDs)./chemistry.workCond.areaOverVolume;
      chemistry.solution.finalDischargeDensity = finalDischargeDensity;
      chemistry.solution.reactionsInfo = struct.empty;
      for id = [chemistry.reactionArray.ID]
        chemistry.solution.reactionsInfo(end+1).reactID = chemistry.reactionArray(id).ID;
        if chemistry.reactionArray(id).isReverse
          chemistry.solution.reactionsInfo(end).rateCoeff = [directRateCoeffs(id) inverseRateCoeffs(id)];
        else
          chemistry.solution.reactionsInfo(end).rateCoeff = directRateCoeffs(id);
        end
        chemistry.solution.reactionsInfo(end).netRate = reactionRates(id);
        chemistry.solution.reactionsInfo(end).energy = chemistry.reactionArray(id).enthalpy;
        chemistry.solution.reactionsInfo(end).description = erase(chemistry.reactionArray(id).descriptionExtended, ' ');
      end
      chemistry.solution.thermalModel = thermalModel;

      %%%%%%%%%%%%%%%%%%%%%%%%% START OF POST-DISCHARGE CODE (WIP) %%%%%%%%%%%%%%%%%%%%%%%%%
      if chemistry.odePostDischargeTime > 0
        % in post-discharge the electron kinetics has a special dependence (checkout kinetics function for details)
        chemistry.electronKineticsDependence = 'postDischarge';
        % in post-discharge quasi-neutrality is imposed during the integration of the rate balance equations
        chemistry.imposeQuasiNeutrality = true;
        % set electron temperature equal to gas temperature
        % (this is overwritten later if the evolution of the eedf is solved or the thermal module is activated)
        chemistry.workCond.electronTemperature = gasTemperatureTime(end)*Constant.boltzmann/Constant.electronCharge;
        % evaluate initial values for the rate coefficients saving IDs of time dependent reactions (post-discharge)
        directRateCoeffs = zeros(1, chemistry.numberOfReactions);
        inverseRateCoeffs = zeros(1, chemistry.numberOfReactions);
        timeDependentReactionIDs = [];
        for reactionID = 1:length(directRateCoeffs)
          % set to zero all "eedf" rate coefficients
          % (this is overwritten later if the evolution of the eedf is solved)
          if strcmp(chemistry.reactionArray(reactionID).type, 'eedf')
            chemistry.reactionArray(reactionID).eedfEquivalent.ineRateCoeff = 0.0;
            chemistry.reactionArray(reactionID).eedfEquivalent.supRateCoeff = 0.0;
            directRateCoeffs(reactionID) = 0.0;
            inverseRateCoeffs(reactionID) = 0.0;
          end
          [directRateCoeffs(reactionID), dependent] = chemistry.rateCoeffFuncHandles{reactionID}(time(end), ...
            absDensitiesTime(end,:)', totalGasDensity, chemistry.reactionArray(reactionID), chemistry.rateCoeffParams{reactionID}, ...
            chemistry);
          if dependent.onTime || dependent.onDensities || (chemistry.includeThermalModel && dependent.onGasTemperature) || ...
            (dependent.onElectronKinetics && chemistry.solveEedf)
            timeDependentReactionIDs(end+1) = reactionID;
          end
          if chemistry.reactionArray(reactionID).isReverse
            inverseRateCoeffs(reactionID) = detailedBalance(chemistry.reactionArray(reactionID), ...
              directRateCoeffs(reactionID), Constant.boltzmannInEV*chemistry.workCond.gasTemperature);
          end
        end
        chemistry.workCond.update('reducedField', 0);
        % flush persistent memory of kinetics function
        kinetics(1, 1, 1, 1, 1, 1, true);
        %chemistry.odeOptions.OutputFcn = @odeProgressBar;          % activate only for debugging purposes
        if chemistry.solveEedf
          chemistry.odeOptions.NonNegative = ...
            1:chemistry.numberOfSpecies+1+chemistry.electronKinetics.energyGrid.cellNumber;
          [timePostDschrg, timeSolutionPostDschrg] = chemistry.odeSolver(@kinetics, ...
            [time(end) time(end)+chemistry.odePostDischargeTime], ...
            [absDensitiesTime(end,:) gasTemperatureTime(end) chemistry.electronKinetics.eedf], chemistry.odeOptions, ...
            directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
          % separate eedf solution
          eedfTimePostDschrg = timeSolutionPostDschrg(:,chemistry.numberOfSpecies+2:end);
        else
          [timePostDschrg, timeSolutionPostDschrg] = chemistry.odeSolver(@kinetics, ...
            [time(end) time(end)+chemistry.odePostDischargeTime], [absDensitiesTime(end,:) gasTemperatureTime(end)], ...
            chemistry.odeOptions, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
        end
        % separate time solutions into its different components
        absDensitiesTimePostDschrg = timeSolutionPostDschrg(:,1:chemistry.numberOfSpecies);
        gasTemperatureTimePostDschrg = timeSolutionPostDschrg(:,chemistry.numberOfSpecies+1);
        % evaluate density of parent states in the post-discharge solution (time evolution) and final total gas density
        totalGasDensity = 0;
        for i = 1:chemistry.numberOfSpecies
          if isempty(chemistry.childIDs{i}) && any(i == chemistry.volumePhaseSpeciesIDs)
            totalGasDensity = totalGasDensity+absDensitiesTimePostDschrg(end,i);
          elseif ~isempty(chemistry.childIDs{i})
            absDensitiesTimePostDschrg(:,i) = 0;
            for j = chemistry.childIDs{i}
              if ~isempty(chemistry.childIDs{j})
                absDensitiesTimePostDschrg(:,j) = 0;
                for k = chemistry.childIDs{j}
                  absDensitiesTimePostDschrg(:,j) = absDensitiesTimePostDschrg(:,j) + absDensitiesTimePostDschrg(:,k);
                end
              end
              absDensitiesTimePostDschrg(:,i) = absDensitiesTimePostDschrg(:,i) + absDensitiesTimePostDschrg(:,j);
            end
          end
        end
        absDensitiesTime = cat(1, absDensitiesTime, absDensitiesTimePostDschrg);
        gasTemperatureTime = [gasTemperatureTime; gasTemperatureTimePostDschrg];
        nearWallTemperatureTimePostDschrg = [];
        wallTemperatureTimePostDschrg = [];
        if chemistry.includeThermalModel
          switch chemistry.thermalModelBoundary
            case 'wall'
              nearWallTemperatureTimePostDschrg = zeros(1,length(gasTemperatureTimePostDschrg));
              nearWallTemperatureTimePostDschrg(1) = nearWallTemperatureTime(end);
              for i = 1:length(gasTemperatureTimePostDschrg)-1
                chemistry.workCond.gasTemperature = gasTemperatureTimePostDschrg(i);
                chemistry.workCond.nearWallTemperature = nearWallTemperatureTimePostDschrg(i);
                [~, thermalModel] = kinetics(timePostDschrg(i), ...
                  [absDensitiesTimePostDschrg(i,:) gasTemperatureTimePostDschrg(i)]', directRateCoeffs, ...
                  inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
                nearWallTemperatureTimePostDschrg(i+1) = thermalModel.nearWallTemperature;
              end
            case 'external'
              nearWallTemperatureTimePostDschrg = zeros(1,length(gasTemperatureTimePostDschrg));
              nearWallTemperatureTimePostDschrg(1) = nearWallTemperatureTime(end);
              wallTemperatureTimePostDschrg = zeros(1,length(gasTemperatureTimePostDschrg));
              wallTemperatureTimePostDschrg(1) = wallTemperatureTime(end);
              for i = 1:length(gasTemperatureTimePostDschrg)-1
                chemistry.workCond.gasTemperature = gasTemperatureTimePostDschrg(i);
                chemistry.workCond.nearWallTemperature = nearWallTemperatureTimePostDschrg(i);
                chemistry.workCond.wallTemperature = wallTemperatureTimePostDschrg(i);
                [~, thermalModel] = kinetics(timePostDschrg(i), ...
                  [absDensitiesTimePostDschrg(i,:) gasTemperatureTimePostDschrg(i)]', directRateCoeffs, ...
                  inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
                nearWallTemperatureTimePostDschrg(i+1) = thermalModel.nearWallTemperature;
                wallTemperatureTimePostDschrg(i+1) = thermalModel.wallTemperature;
              end
          end
        end
        nearWallTemperatureTime = [nearWallTemperatureTime nearWallTemperatureTimePostDschrg];
        wallTemperatureTime = [wallTemperatureTime wallTemperatureTimePostDschrg];
        time = cat(1,time,timePostDschrg);
      end
      %%%%%%%%%%%%%%%%%%%%%%%%% END OF POST-DISCHARGE CODE (WIP) %%%%%%%%%%%%%%%%%%%%%%%%%

      % store time dependent solution of the heavy species kinetics in the chemistry properties
      absDensitiesTime(:,chemistry.surfacePhaseSpeciesIDs) = ...
        absDensitiesTime(:,chemistry.surfacePhaseSpeciesIDs)./chemistry.workCond.areaOverVolume;
      chemistry.solution.time = time;
      chemistry.solution.densitiesTime = absDensitiesTime;
      chemistry.solution.gasTemperatureTime = gasTemperatureTime;
      chemistry.solution.nearWallTemperatureTime = nearWallTemperatureTime;
      chemistry.solution.wallTemperatureTime = wallTemperatureTime;

      % broadcast obtention of a solution for the chemistry equation
      notify(chemistry, 'obtainedNewChemistrySolution');

    end

    function solvePulsed(chemistry)
      
      % Temporal integration of rate balance equations for pulsed simulations
      [time, absDensitiesTime, gasTemperatureTime, electronTemperatureTime, directRateCoeffs, inverseRateCoeffs, ...
        timeDependentReactionIDs] = chemistry.integrateRateBalanceEquations();

      % evaluate density of parent states in and total gas density at each time point
      totalGasDensity = zeros(size(time));
      electronDensity = zeros(size(time));
      for i = 1:chemistry.numberOfSpecies
        if isempty(chemistry.childIDs{i}) && any(i == chemistry.volumePhaseSpeciesIDs)
          totalGasDensity = totalGasDensity+absDensitiesTime(:,i);
          if strcmp(chemistry.stateArray(i).ionCharg, '+')
            electronDensity = electronDensity+absDensitiesTime(:,i);
          elseif strcmp(chemistry.stateArray(i).ionCharg, '-')
            electronDensity = electronDensity-absDensitiesTime(:,i);
          end
        elseif ~isempty(chemistry.childIDs{i})
          absDensitiesTime(:,i) = 0;
          for j = chemistry.childIDs{i}
            if ~isempty(chemistry.childIDs{j})
              absDensitiesTime(:,j) = 0;
              for k = chemistry.childIDs{j}
                absDensitiesTime(:,j) = absDensitiesTime(:,j) + absDensitiesTime(:,k);
              end
            end
            absDensitiesTime(:,i) = absDensitiesTime(:,i) + absDensitiesTime(:,j);
          end
        end
      end

      % evaluate intermediate temperatures at each time point (if thermal model is activated)
      nearWallTemperatureTime = [];
      wallTemperatureTime = [];
      if chemistry.includeThermalModel
        switch chemistry.thermalModelBoundary
          case 'wall'
            nearWallTemperatureTime = zeros(1,length(gasTemperatureTime));
            nearWallTemperatureTime(1) = chemistry.initialNearWallTemperature;
            for i = 1:length(gasTemperatureTime)-1
              chemistry.workCond.gasTemperature = gasTemperatureTime(i);
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
              [~, thermalModel] = kinetics(time(i), [absDensitiesTime(i,:) gasTemperatureTime(i)]', ...
                directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
              nearWallTemperatureTime(i+1) = thermalModel.nearWallTemperature;
            end
          case 'external'
            nearWallTemperatureTime = zeros(1,length(gasTemperatureTime));
            nearWallTemperatureTime(1) = chemistry.initialNearWallTemperature;
            wallTemperatureTime = zeros(1,length(gasTemperatureTime));
            wallTemperatureTime(1) = chemistry.initialWallTemperature;
            for i = 1:length(gasTemperatureTime)-1
              chemistry.workCond.gasTemperature = gasTemperatureTime(i);
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
              chemistry.workCond.wallTemperature = wallTemperatureTime(i);
              [~, thermalModel] = kinetics(time(i), [absDensitiesTime(i,:) gasTemperatureTime(i)]', ...
                directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
              nearWallTemperatureTime(i+1) = thermalModel.nearWallTemperature;
              wallTemperatureTime(i+1) = thermalModel.wallTemperature;
            end
        end
      else
        thermalModel = struct.empty;
      end


      % evaluate rates at each time point
      for i = 1:length(time)
        % evaluate working conditions at each time point (for the evaluation of the rates)
        chemistry.workCond.electronDensity = electronDensity(i);
        if strcmp(chemistry.electronKineticsDependence, 'quasiStationary') && strcmp(chemistry.lookUpMethod, 'localEnergy')
          chemistry.workCond.electronTemperature = electronTemperatureTime(i);
        end
        chemistry.workCond.gasTemperature = gasTemperatureTime(i);
        chemistry.workCond.update('reducedField', chemistry.pulseFunction(time(i), chemistry.pulseFunctionParameters));
        if chemistry.includeThermalModel
          switch chemistry.thermalModelBoundary
            case 'wall'
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
            case 'external'
              chemistry.workCond.nearWallTemperature = nearWallTemperatureTime(i);
              chemistry.workCond.wallTemperature = wallTemperatureTime(i);
          end
        end
        % evaluate time dependent rate coefficients at each time point
        for j = 1:length(timeDependentReactionIDs)
          reactionID = timeDependentReactionIDs(j);
          directRateCoeffs(reactionID) = chemistry.rateCoeffFuncHandles{reactionID}(time(i), ...
            absDensitiesTime(i,:)', totalGasDensity(i), chemistry.reactionArray(reactionID), chemistry.rateCoeffParams{reactionID}, chemistry);
          if chemistry.reactionArray(reactionID).isReverse
            inverseRateCoeffs(reactionID) = detailedBalance(chemistry.reactionArray(reactionID), ...
              directRateCoeffs(reactionID), Constant.boltzmannInEV*chemistry.workCond.gasTemperature);
          end
        end
        % evaluate reaction rate at each time point
        reactionRates = directRateCoeffs.*prod(repmat(absDensitiesTime(i,:)',[1 chemistry.numberOfReactions]).^ ...
          (chemistry.reactantStoiCoeffs+chemistry.catalystStoiCoeffs),1).* ...
          chemistry.workCond.electronDensity.^chemistry.reactantElectrons- ...
          inverseRateCoeffs.*prod(repmat(absDensitiesTime(i,:)',[1 chemistry.numberOfReactions]).^ ...
          (chemistry.productStoiCoeffs+chemistry.catalystStoiCoeffs),1).* ...
          chemistry.workCond.electronDensity.^chemistry.productElectrons;
        for j = chemistry.gasStabilisedReactionIDs
          reactionRates(j) = reactionRates(j)*totalGasDensity(i);
        end
        % store discharge solution of the heavy species kinetics in the chemistry properties at each time point
        chemistry.workCond.currentTime = time(i);
        finalDischargeDensity = absDensitiesTime(i,:);
        finalDischargeDensity(chemistry.surfacePhaseSpeciesIDs) = ...
          finalDischargeDensity(chemistry.surfacePhaseSpeciesIDs)./chemistry.workCond.areaOverVolume;
        chemistry.solution.finalDischargeDensity = finalDischargeDensity;
        chemistry.solution.reactionsInfo = struct.empty;
        for id = [chemistry.reactionArray.ID]
          chemistry.solution.reactionsInfo(end+1).reactID = chemistry.reactionArray(id).ID;
          if chemistry.reactionArray(id).isReverse
            chemistry.solution.reactionsInfo(end).rateCoeff = [directRateCoeffs(id) inverseRateCoeffs(id)];
          else
            chemistry.solution.reactionsInfo(end).rateCoeff = directRateCoeffs(id);
          end
          chemistry.solution.reactionsInfo(end).netRate = reactionRates(id);
          chemistry.solution.reactionsInfo(end).energy = chemistry.reactionArray(id).enthalpy;
          chemistry.solution.reactionsInfo(end).description = erase(chemistry.reactionArray(id).descriptionExtended, ' ');
        end
        chemistry.solution.thermalModel = thermalModel;
        
        % store time dependent solution of the heavy species kinetics in the chemistry properties
        absDensitiesTime(:,chemistry.surfacePhaseSpeciesIDs) = ...
          absDensitiesTime(:,chemistry.surfacePhaseSpeciesIDs)./chemistry.workCond.areaOverVolume;
        chemistry.solution.time = time;
        chemistry.solution.densitiesTime = absDensitiesTime;
        chemistry.solution.gasTemperatureTime = gasTemperatureTime;
        chemistry.solution.nearWallTemperatureTime = nearWallTemperatureTime;
        chemistry.solution.wallTemperatureTime = wallTemperatureTime;
        chemistry.neutralityRelErrorCurrent = 0;
  
        % broadcast obtention of a solution for the chemistry model
        notify(chemistry, 'obtainedNewChemistrySolution');
      end

    end

    function [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
        globalCycle(chemistry)

      % thresholds and other data for cycles
      maxGlobalRelError = chemistry.globalRelError;
      maxGlobalIterations = chemistry.globalMaxIterations;
      currentGlobalIteration = 1;
      if ~isempty(chemistry.globalIterationCurrent)
        currentGlobalIteration = chemistry.globalIterationCurrent;
      end
      currentGlobalRelError = [];

      while currentGlobalIteration<=maxGlobalIterations

        [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
          chemistry.neutralityCycle();

        if ~isempty(chemistry.electronKinetics) && currentGlobalIteration<maxGlobalIterations

          % update gas temperature dependencies of the electron kinetics (in case thermal model is active)
          if chemistry.includeThermalModel
            % update gas temperature used in the electron kinetics
            chemistry.electronKineticsGasTemperature = gasTemperatureTime(end);
            chemistry.workCond.gasTemperature = gasTemperatureTime(end);
            % update possible populations depending on the gas temperature
            for gas = chemistry.electronKinetics.gasArray
              for state = gas.stateArray
                if ~isempty(state.populationFunc) && isempty(state.chemEquivalent)
                  for parameter = state.populationParams
                    if strcmp(parameter{1}, 'gasTemperature')
                      state.evaluatePopulation(chemistry.workCond);
                      break;
                    end
                  end
                end
              end
            end
            % update the densities of states accordingly
            for gas = chemistry.electronKinetics.gasArray
              for state = gas.stateArray
                state.evaluateDensity();
              end
            end
          end

          % evaluate electron kinetics (eedf) densities (final time)
          eedfAbsDensities = absDensitiesTime(end,:);
          eedfTotalGasDensity = 0;
          for i = chemistry.statesIDsToUpdateInElectronKinetics
            if isempty(chemistry.childIDs{i})
              eedfTotalGasDensity = eedfTotalGasDensity + eedfAbsDensities(i);
            else
              eedfAbsDensities(i) = 0;
              for j = intersect(chemistry.childIDs{i}, chemistry.statesIDsToUpdateInElectronKinetics)
                if ~isempty(chemistry.childIDs{j})
                  eedfAbsDensities(j) = 0;
                  for k = intersect(chemistry.childIDs{j}, chemistry.statesIDsToUpdateInElectronKinetics)
                    eedfAbsDensities(j) = eedfAbsDensities(j) + eedfAbsDensities(k);
                  end
                end
                eedfAbsDensities(i) = eedfAbsDensities(i) + eedfAbsDensities(j);
              end
            end
          end

          % NEW - Update chem densities for next iteration
          %for i = 1:chemistry.numberOfSpecies
          %  chemistry.stateArray(i).density = absDensitiesTime(end,i)/eedfTotalGasDensity;
          %end

          % update densities of states in the the electron kinetics gas mixture (normalized to gas density)
          for i = chemistry.statesIDsToUpdateInElectronKinetics
            chemistry.stateArray(i).eedfEquivalent.density = eedfAbsDensities(i)/eedfTotalGasDensity;
          end

          % renormalize electron kinetics gas mixture
          for i = chemistry.gasesIDsToUpdateInElectronKinetics
            chemistry.gasArray(i).eedfEquivalent.renormalizeWithDensities();
            % check for the distribution of states to be properly normalised
            chemistry.gasArray(i).eedfEquivalent.checkPopulationNorms();
          end

          % update the density dependencies of the electron kinetics (this also updates gas temperature dependencies)
          chemistry.electronKinetics.updateDensityDependencies();

          % solve the electron kinetics with new densities
          chemistry.electronKinetics.solve();

          % evaluate global relative error (mean of relative errors of the reduced diffusion coefficient and mobility)
          currentGlobalRelError = ...
            (chemistry.electronTransportProperties.reducedDiffCoeff/chemistry.electronKinetics.swarmParam.redDiffCoeff + ...
            chemistry.electronTransportProperties.reducedMobility/chemistry.electronKinetics.swarmParam.redMobility)/2 - 1;

          % broadcast results of this iteration
          chemistry.globalIterationCurrent = currentGlobalIteration;
          chemistry.globalRelErrorCurrent = currentGlobalRelError;
          notify(chemistry, 'newGlobalCycleIteration');

          % prepare next iteration (new solution of the boltzmann equation)
          chemistry.electronTransportProperties.reducedDiffCoeff = chemistry.electronKinetics.swarmParam.redDiffCoeff;
          chemistry.electronTransportProperties.reducedMobility = chemistry.electronKinetics.swarmParam.redMobility;
          if abs(currentGlobalRelError)>maxGlobalRelError
            currentGlobalIteration = currentGlobalIteration+1;
            if currentGlobalIteration==maxGlobalIterations
              error('Maximum number of iterations reached for the global cycle without convergence.')
            end
          else
            currentGlobalIteration = maxGlobalIterations;
          end
        else
          break;
        end

      end

    end

    function [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
        neutralityCycle(chemistry)

      % thresholds and other data for cycles
      maxNeutralityRelError = chemistry.neutralityRelError;
      maxNeutralityIterations = chemistry.neutralityMaxIterations;
      currentNeutralityIteration = 1;
      if ~isempty(chemistry.neutralityIterationCurrent)
        currentNeutralityIteration = chemistry.neutralityIterationCurrent;
      end

      currentNeutralityRelError = [];
      neutralityRelErrorAll = [];
      neutralityIterationAll = [];
      excitationParameterAll = [];
      while currentNeutralityIteration<=maxNeutralityIterations

        [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
          chemistry.pressureCycle();

        if ~isempty(chemistry.electronKinetics)
          % evaluate the neutrality relative error
          electronDensityObtained = 0;
          for i = 1:chemistry.numberOfSpecies
            chargStr = chemistry.stateArray(i).ionCharg;
            if ~isempty(chargStr)
              electronDensityObtained = electronDensityObtained + ...
                length(regexp(chargStr,'+'))*absDensitiesTime(end,i);
              electronDensityObtained = electronDensityObtained - ...
                length(regexp(chargStr,'-'))*absDensitiesTime(end,i);
            end
          end
          switch chemistry.convergenceParameter
            case 'electronDensity'
              electronDensityPrescribed = chemistry.workCond.electronDensity;
            case 'dischargeCurrent'
              electronDensityPrescribed = chemistry.workCond.dischargeCurrent / ( Constant.electronCharge* ...
                chemistry.electronKinetics.swarmParam.driftVelocity*pi*chemistry.workCond.chamberRadius^2);
            case 'dischargePowerDensity'
              electronDensityPrescribed = chemistry.workCond.dischargePowerDensity / (chemistry.workCond.gasDensity * ...
                chemistry.electronKinetics.power.field * Constant.electronCharge);
          end
          currentNeutralityRelError = (electronDensityObtained - electronDensityPrescribed)/electronDensityPrescribed;
          neutralityIterationAll(end+1) = currentNeutralityIteration;
          neutralityRelErrorAll(end+1) = currentNeutralityRelError;

          % broadcast results of this iteration
          chemistry.neutralityIterationCurrent = currentNeutralityIteration;
          chemistry.neutralityRelErrorCurrent = currentNeutralityRelError;
          notify(chemistry, 'newNeutralityCycleIteration');

          % prepare next iteration (new solution of the boltzmann equation)
          if abs(currentNeutralityRelError)>maxNeutralityRelError
            if strcmp(chemistry.convergenceParameter, 'dischargeCurrent') || ...
                strcmp(chemistry.convergenceParameter, 'dischargePowerDensity')
              chemistry.workCond.update('electronDensity', electronDensityPrescribed);
              % flush persistent variables of kinetics function (to ensure new value of ne^stoichiometric coeff)
              kinetics(1, 1, 1, 1, 1, 1, true);
            end
            switch class(chemistry.electronKinetics)
              case 'Boltzmann'
                excitationParameterAll(end+1) = chemistry.workCond.reducedField;
                chemistry.workCond.update('reducedField', iterateOverParameter(neutralityIterationAll, ...
                  neutralityRelErrorAll, excitationParameterAll, true, true));
              case 'PrescribedEedf'
                excitationParameterAll(end+1) = chemistry.workCond.electronTemperature;
                chemistry.workCond.update('electronTemperature', iterateOverParameter(neutralityIterationAll, ...
                  neutralityRelErrorAll, excitationParameterAll, true, true));
            end
            chemistry.workCond.update('gasTemperature', chemistry.electronKineticsGasTemperature);
            chemistry.electronKinetics.solve();
            chemistry.electronTransportProperties.reducedDiffCoeff = chemistry.electronKinetics.swarmParam.redDiffCoeff;
            chemistry.electronTransportProperties.reducedMobility = chemistry.electronKinetics.swarmParam.redMobility;
            currentNeutralityIteration = currentNeutralityIteration+1;
            if currentNeutralityIteration>maxNeutralityIterations
              error('Maximum number of iterations reached for the neutrality cycle without convergence.')
            end
          else
            break;
          end
        else
          break;
        end

      end

    end

    function [time, absDensitiesTime, gasTemperatureTime, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
        pressureCycle(chemistry)

      % thresholds and other data for cycles
      maxPressureRelError = chemistry.pressureRelError;
      maxPressureIterations = chemistry.pressureMaxIterations;
      currentPressureIteration = 1;
      if ~isempty(chemistry.pressureIterationCurrent)
        currentPressureIteration = chemistry.pressureIterationCurrent;
      end

      currentPressureRelError = [];
      pressureRelErrorAll = [];
      pressureIterationAll = [];
      initialGasPressureAll = [];
      while currentPressureIteration<=maxPressureIterations

        [time, absDensitiesTime, gasTemperatureTime, ~, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs] = ...
          chemistry.integrateRateBalanceEquations();

        % evaluate the pressure relative error
        finalGasDensity = 0;
        for i = chemistry.volumePhaseSpeciesIDs
          if isempty(chemistry.childIDs{i})
            finalGasDensity = finalGasDensity+absDensitiesTime(end,i);
          end
        end
        finalGasPressure = finalGasDensity*Constant.boltzmann*gasTemperatureTime(end);
        currentPressureRelError = (finalGasPressure - chemistry.targetGasPressure)/chemistry.targetGasPressure;
        pressureIterationAll(end+1) = currentPressureIteration;
        pressureRelErrorAll(end+1) = currentPressureRelError;
        initialGasPressureAll(end+1) = chemistry.initialGasPressure;
        
        % broadcast results of this iteration
        chemistry.pressureIterationCurrent = currentPressureIteration;
        chemistry.pressureRelErrorCurrent = currentPressureRelError;
        notify(chemistry, 'newPressureCycleIteration');

        % prepare next iteration (new value of the initial gas pressure + update iteration index)
        if abs(currentPressureRelError)>maxPressureRelError
          chemistry.initialGasPressure = iterateOverParameter(pressureIterationAll, pressureRelErrorAll, ...
            initialGasPressureAll, true, true);
          currentPressureIteration = currentPressureIteration+1;
          if currentPressureIteration>maxPressureIterations
            error('Maximum number of iterations reached for the pressure cycle without convergence.')
          end
        else
          % update final gas pressure and temperature in the workCond object with final discharge values for neutrality cycle
          chemistry.workCond.update('gasPressure', finalGasPressure);
          chemistry.workCond.update('gasTemperature', gasTemperatureTime(end));
          break;
        end

      end

    end

    function [time, absDensitiesTime, gasTemperatureTime, electronTemperatureTime, directRateCoeffs, ...
        inverseRateCoeffs, timeDependentReactionIDs] = integrateRateBalanceEquations(chemistry)
      
      % evaluate initial densities of species (volume phase and surface phase)
      initialGasDensity = chemistry.initialGasPressure/(Constant.boltzmann*chemistry.initialGasTemperature);
      initialAbsDensities = zeros(1, chemistry.numberOfSpecies);
      %%%%%%%%%%%%%%%%%%%% OLD METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      initialAbsDensities(chemistry.volumePhaseSpeciesIDs) = initialGasDensity.*...
        chemistry.initialDensities(chemistry.volumePhaseSpeciesIDs);
      initialAbsDensities(chemistry.surfacePhaseSpeciesIDs) = (chemistry.workCond.areaOverVolume*...
        chemistry.workCond.surfaceSiteDensity).*chemistry.initialDensities(chemistry.surfacePhaseSpeciesIDs);
      %%%%% NEW METHOD - Start heavy-species kinetics at previous global cycle (final discharge solution) %%%%%%%%%%%
      % initialAbsDensities(chemistry.volumePhaseSpeciesIDs) = initialGasDensity.*...
      %   [chemistry.stateArray(chemistry.volumePhaseSpeciesIDs).density];
      % initialAbsDensities(chemistry.surfacePhaseSpeciesIDs) = (chemistry.workCond.areaOverVolume*...
      %   chemistry.workCond.surfaceSiteDensity).*[chemistry.stateArray(chemistry.surfacePhaseSpeciesIDs).density];

      % set initial value for the temperatures
      chemistry.workCond.gasTemperature = chemistry.initialGasTemperature;
      if chemistry.includeThermalModel
        chemistry.workCond.nearWallTemperature = chemistry.initialNearWallTemperature;
        if strcmp(chemistry.thermalModelBoundary, 'external')
          chemistry.workCond.wallTemperature = chemistry.initialWallTemperature;
        end
      end

      % evaluate initial values for the rate coefficients saving IDs of time dependent reactions
      directRateCoeffs = zeros(1, chemistry.numberOfReactions);
        inverseRateCoeffs = zeros(1, chemistry.numberOfReactions);
        timeDependentReactionIDs = [];
      for reactionID = 1:length(directRateCoeffs)
        [directRateCoeffs(reactionID), dependent] = chemistry.rateCoeffFuncHandles{reactionID}(0, ...
          initialAbsDensities', initialGasDensity, chemistry.reactionArray(reactionID), chemistry.rateCoeffParams{reactionID}, chemistry);
        if dependent.onTime || dependent.onDensities || (chemistry.includeThermalModel && dependent.onGasTemperature) || ...
            (chemistry.isPulsed && dependent.onElectronKinetics)
          timeDependentReactionIDs(end+1) = reactionID;
        end
        if chemistry.reactionArray(reactionID).isReverse
          inverseRateCoeffs(reactionID) = detailedBalance(chemistry.reactionArray(reactionID), ...
            directRateCoeffs(reactionID), Constant.boltzmannInEV*chemistry.workCond.gasTemperature);
        end
      end

      % select variables to the solved by the ODE solver (densities of all species + gas temperature always)
      initialValues = [initialAbsDensities chemistry.initialGasTemperature];
      if chemistry.isPulsed && strcmp(chemistry.electronKineticsDependence, 'quasiStationary') && ...
        strcmp(chemistry.lookUpMethod, 'localEnergy')
        initialValues = [initialValues chemistry.workCond.electronTemperature];
      end

      % select time span for the ODE solver
      if chemistry.isPulsed
        % evaluate sampling points prescribed by the user
        if strcmp(chemistry.pulseSamplingType, 'linspace')
          tspan = [0 linspace(chemistry.pulseFirstStep, chemistry.pulseFinalTime, chemistry.pulseSamplingPoints)];
        elseif strcmp(chemistry.pulseSamplingType, 'logspace')
          tspan = [0 logspace(log10(chemistry.pulseFirstStep), log10(chemistry.pulseFinalTime), ...
            chemistry.pulseSamplingPoints)];
        end
      else 
        tspan = [0 chemistry.odeDischargeTime];
      end

      % call the ODE solver
      start = tic;
      notify(chemistry, 'genericStatusMessage', ...
        StatusEventData('\t- Integrating particle rate-balance equations ...\n', 'status'));
      % uncomment the next line for testing event detection functionalities (NOT FOR PRODUCTION!)
      % chemistry.odeOptions.Events = @odeEventFunction;
      [time, timeSolution] = chemistry.odeSolver(@kinetics, tspan, initialValues, ...
        chemistry.odeOptions, directRateCoeffs, inverseRateCoeffs, timeDependentReactionIDs, chemistry, false);
      str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
      notify(chemistry, 'genericStatusMessage', StatusEventData(str, 'status'));

      % separate time solutions into its different components
      absDensitiesTime = timeSolution(:,1:chemistry.numberOfSpecies);
      gasTemperatureTime = timeSolution(:,chemistry.numberOfSpecies+1);
      if chemistry.isPulsed && strcmp(chemistry.electronKineticsDependence, 'quasiStationary') && ...
        strcmp(chemistry.lookUpMethod, 'localEnergy')
        electronTemperatureTime = timeSolution(:,chemistry.numberOfSpecies+2);
      else
        electronTemperatureTime = [];
      end

    end

  end

end

function [derivatives, thermalModel] = kinetics(time, variables, directRateCoeffs, inverseRateCoeffs, ...
  timeDependentReactionIDs, chemistry, clearPersistentVars)
  
  % define persistent variables (do not change during the simulation)
  persistent workCond;
  persistent reactionArray;
  persistent gasArray;
  persistent stateArray;
  persistent electronKinetics;
  persistent rateCoeffFuncHandles;
  persistent rateCoeffParams;
  persistent chamberVolume;
  persistent electronDensitiesReactant;
  persistent electronDensitiesProduct;
  persistent totalReactantStoiCoeffs;
  persistent totalProductStoiCoeffs;
  persistent productMinusReactantStoiCoeffs;
  persistent isTimeDependentReactionReverse;
  persistent reactionEnthalpies;
  persistent numberOfGases;
  persistent numberOfSpecies;
  persistent volumePhaseSpeciesIDs;
  persistent childlessVolumePhaseSpeciesIDs;
  persistent surfacePhaseSpeciesIDs;
  persistent numberOfReactions;
  persistent massArray;                       % array with masses of all gases (needed for the thermal model)
  persistent charLengthThermalModelSquared;   % characteristic length for the thermal model (see documentation)
  persistent lstdsteq1;  
  persistent lstdstgt1;
  persistent lstisteq1;
  persistent lstistgt1;
  persistent drateaux;
  persistent irateaux;
  persistent firstReactionFlowID;
  persistent includeJouleHeating;
  persistent includeGasRotHeating;
  persistent includeGasVibHeating;
  persistent rotStatesPresent;
  persistent vibStatesPresent;
  
  % flush persistent memory for a new simulation (needed before post-discharge integration)
  if clearPersistentVars
    vars = whos;
    vars = vars([vars.persistent]);
    varName = {vars.name};
    clear(varName{:});
    derivatives = [];
    return
  end

  % initialize values of persistent variables
  if isempty(electronDensitiesReactant)
    workCond = chemistry.workCond;
    reactionArray = chemistry.reactionArray;
    gasArray = chemistry.gasArray;
    stateArray = chemistry.stateArray;
    electronKinetics = chemistry.electronKinetics;
    rateCoeffFuncHandles = chemistry.rateCoeffFuncHandles;   
    rateCoeffParams = chemistry.rateCoeffParams;
    chamberVolume = pi*(workCond.chamberRadius^2)*workCond.chamberLength;
    electronDensitiesReactant = workCond.electronDensity.^chemistry.reactantElectrons';
    electronDensitiesProduct = workCond.electronDensity.^chemistry.productElectrons';
    totalReactantStoiCoeffs = (chemistry.reactantStoiCoeffs+chemistry.catalystStoiCoeffs)';
    totalProductStoiCoeffs = (chemistry.productStoiCoeffs+chemistry.catalystStoiCoeffs)';
    productMinusReactantStoiCoeffs = chemistry.productStoiCoeffs - chemistry.reactantStoiCoeffs;
    numberOfGases = length(chemistry.gasArray);
    numberOfSpecies = chemistry.numberOfSpecies;
    volumePhaseSpeciesIDs = chemistry.volumePhaseSpeciesIDs;
    childlessVolumePhaseSpeciesIDs = chemistry.childlessVolumePhaseSpeciesIDs;
    surfacePhaseSpeciesIDs = chemistry.surfacePhaseSpeciesIDs;
    numberOfReactions = chemistry.numberOfReactions;
    lstdsteq1=find(totalReactantStoiCoeffs==1);   %List of positions with stoichiometric coeff == 1 in direct reactions
    lstdstgt1=find(totalReactantStoiCoeffs>1);    %List of positions with stoichiometric coeff > 1 in direct reactions
    lstisteq1=find(totalProductStoiCoeffs==1);
    lstistgt1=find(totalProductStoiCoeffs>1);
    drateaux=ones(size(totalReactantStoiCoeffs));  %Initialize density matrix with 1 for direct rates calculation
    irateaux=ones(size(totalProductStoiCoeffs));
    isTimeDependentReactionReverse = false(size(timeDependentReactionIDs));
    for i = 1:length(timeDependentReactionIDs)
      reactionID = timeDependentReactionIDs(i);
      if reactionArray(reactionID).isReverse
        isTimeDependentReactionReverse(i) = true;
      end
    end
    for i = 1:numberOfReactions
      if strcmp(reactionArray(i).type,'outFlow')
        if isempty(firstReactionFlowID)
          firstReactionFlowID = i;
        end
      end
    end
    if chemistry.includeThermalModel
      % saving masses of volume species
      massArray = zeros(1, numberOfGases);
      for i = 1: numberOfGases
        if gasArray(i).isVolumeSpecies
          massArray(i) = gasArray(i).mass;
        end
      end
      % saving reaction enthalpies in a simple and fast array of doubles
      reactionEnthalpies = zeros(size(reactionArray));
      for i = 1:numberOfReactions
        reactionEnthalpies(i) = reactionArray(i).enthalpy;
      end
      % evaluate characteristic length for the thermal model (radial profile for gas temperature)
      charLengthThermalModelSquared = workCond.chamberRadius^2/8;
      % saving booleans deciding the inclusion of heating terms due to Joule, rotational and vibrational mechanisms
      includeJouleHeating = chemistry.includeJouleHeating;
      includeGasRotHeating = chemistry.includeGasRotHeating;
      includeGasVibHeating = chemistry.includeGasVibHeating;
      rotStatesPresent = any(strcmp({chemistry.stateArray.type}, 'rot'));
      vibStatesPresent = any(strcmp({chemistry.stateArray.type}, 'vib'));
    end
  end

  % separate variables
  densities = variables(1:numberOfSpecies);
  gasTemperature = variables(numberOfSpecies+1);
  if strcmp(chemistry.electronKineticsDependence, 'postDischarge') && chemistry.solveEedf
    eedf = variables(numberOfSpecies+2:end);
  elseif strcmp(chemistry.electronKineticsDependence, 'quasiStationary') && strcmp(chemistry.lookUpMethod, 'localEnergy')
    workCond.electronTemperature = variables(numberOfSpecies+2);
  end

  % evaluate density of parent states, individual gases and total gas density
  totalGasDensity = 0;
  individualGasDensity = zeros(1, numberOfGases);
  for i = 1:numberOfSpecies
    if isempty(stateArray(i).childArray) && stateArray(i).isVolumeSpecies
      totalGasDensity = totalGasDensity+densities(i);
      individualGasDensity(chemistry.gasIDs(i)) = individualGasDensity(chemistry.gasIDs(i))+densities(i);
    elseif ~isempty(stateArray(i).childArray)
      densities(i) = 0;
      IDs1 = chemistry.childIDs{i};
      for j = 1:length(IDs1)
        if ~isempty(stateArray(IDs1(j)).childArray)
          densities(IDs1(j)) = 0;
          IDs2 = chemistry.childIDs{IDs1(j)};
          for k = 1:length(IDs2)
            densities(IDs1(j)) = densities(IDs1(j)) + densities(IDs2(k));
          end
        end
        densities(i) = densities(i) + densities(IDs1(j));
      end
    end
  end

  % evaluate thermal balance equation parameters (if thermal model is activated) --- 1st part of the calculations
  if chemistry.includeThermalModel
    % update working conditions object with current gas temperature
    workCond.gasTemperature = gasTemperature;
    % reevaluate thermal properties of gases that might depend on the gas temperature
    heatCapacityArray = zeros(1,numberOfGases);
    thermalConductivityArray = zeros(1,numberOfGases);
    for i = 1:numberOfGases
      if gasArray(i).isVolumeSpecies
        heatCapacityArray(i) = gasArray(i).evaluateHeatCapacity(workCond);
        thermalConductivityArray(i) = gasArray(i).evaluateThermalConductivity(workCond);
      end
    end
    % evaluate total system heat capacity 
    heatCapacity = dot(individualGasDensity,heatCapacityArray/Constant.avogadro);
    % evaluate total system heat conductivity (can not be vectorized because of exceptions in the loops)
    euckenFactorArray = 0.115+0.354*heatCapacityArray/Constant.idealGas*Constant.electronCharge;
    thermalConductivity = 0;
    for i = 1:numberOfGases
      if individualGasDensity(i) == 0
        continue
      end
      aux = 1;
      for j = 1:numberOfGases
        if j == i || individualGasDensity(j) == 0
          continue
        end
        aux = aux + 1.065/(2*sqrt(2*(1+massArray(i)/massArray(j))))*(1+sqrt(thermalConductivityArray(i)* ...
          euckenFactorArray(j)/(thermalConductivityArray(j)*euckenFactorArray(i))*sqrt(massArray(i)/massArray(j))))^2* ...
          individualGasDensity(j)/individualGasDensity(i);
      end
      thermalConductivity = thermalConductivity + thermalConductivityArray(i) / aux;
    end
    % update near wall temperature in case of 'wall' boundary condition
    if strcmp(chemistry.thermalModelBoundary, 'wall')
      aux = 4*thermalConductivity/workCond.chamberRadius;
      hint = chemistry.intConvCoeff;
      nearWallTemperature = (aux*gasTemperature/hint+workCond.wallTemperature) / (aux/hint+1);
      workCond.nearWallTemperature = nearWallTemperature;
      wallTemperature = workCond.wallTemperature;
    end
  end
  
  % solve evolution of electron density by imposing quasi-neutrality
  if chemistry.imposeQuasiNeutrality
    % evaluate current electron density (imposing quasi-neutrality)
    electronDensity = 0;
    for i = 1:numberOfSpecies
      if strcmp(stateArray(i).ionCharg, '+') && stateArray(i).isVolumeSpecies && isempty(stateArray(i).childArray)
        electronDensity = electronDensity + densities(i);
      elseif strcmp(stateArray(i).ionCharg, '-') && stateArray(i).isVolumeSpecies && isempty(stateArray(i).childArray)
        electronDensity = electronDensity - densities(i);
      end
    end
    workCond.electronDensity = electronDensity;
    electronDensitiesReactant = electronDensity.^chemistry.reactantElectrons';
    electronDensitiesProduct = electronDensity.^chemistry.productElectrons';
  end

  % HERE THE DIFFERENT MODELS FOR THE ELECTRON KINETICS SHOULD BE IMPLEMENTED
  switch chemistry.electronKineticsDependence
    case 'steadyState'
      % electron kinetics is frozen during the integration of the rate balance equations. Nothing to do here. 
      % (placeholder for possible future implementations)
      
    case 'quasiStationary'
      % during pulsed simulations, electron kinetics related info is obtained from the look-up tables
      workCond.update('reducedField', chemistry.pulseFunction(time, chemistry.pulseFunctionParameters));
      switch chemistry.lookUpMethod
        case 'localField'
          lookUpXValues = chemistry.lookUpTableRedFieldValues;
          lookUpQueryPoint = workCond.reducedField;
        case 'localEnergy'
          lookUpXValues = chemistry.lookUpTableEleTempValues;
          lookUpQueryPoint = workCond.electronTemperature;
      end
      DN = interp1(lookUpXValues, chemistry.lookUpTableRedDiffValues, lookUpQueryPoint);
      muN = interp1(lookUpXValues, chemistry.lookUpTableRedMobValues, lookUpQueryPoint);
      if isnan(DN) || isnan(muN)
        error(['Reduced transport parameters out of look-up table range.\n'
          'Please check the look-up table and the range of the reduced field during the simulation.'])
      end
      chemistry.electronTransportProperties.reducedDiffCoeff = DN;
      chemistry.electronTransportProperties.reducedMobility = muN;
    case 'postDischarge'
      % during post-discharge phase, electron kinetics is solved at zero reduced electric field and kept frozen during the
      % integration of the rate balance equations (except chemisty.solveEedf == true, EXPERIMENTAL FEATURE). If the 
      % thermal model is activated the electron temperature is thermalized with the gas temperature at each time step.
      if chemistry.solveEedf
        % update gas temperature dependencies of the electron kinetics (in case thermal model is active)
        if chemistry.includeThermalModel
          % update possible populations depending on the gas temperature
          for gas = electronKinetics.gasArray
            for state = gas.stateArray
              if ~isempty(state.populationFunc) && isempty(state.chemEquivalent)
                for parameter = state.populationParams
                  if strcmp(parameter{1}, 'gasTemperature')
                    state.evaluatePopulation(workCond);
                    break;
                  end
                end
              end
            end
          end
          % update the densities of states accordingly
          for gas = electronKinetics.gasArray
            for state = gas.stateArray
              state.evaluateDensity();
            end
          end
        end
        % evaluate electron kinetics (eedf) densities
        eedfAbsDensities = densities;
        eedfTotalGasDensity = 0;
        for i = chemistry.statesIDsToUpdateInElectronKinetics
          if isempty(chemistry.childIDs{i})
            eedfTotalGasDensity = eedfTotalGasDensity + eedfAbsDensities(i);
          else
            eedfAbsDensities(i) = 0;
            for j = intersect(chemistry.childIDs{i}, chemistry.statesIDsToUpdateInElectronKinetics)
              if ~isempty(chemistry.childIDs{j})
                eedfAbsDensities(j) = 0;
                for k = intersect(chemistry.childIDs{j}, chemistry.statesIDsToUpdateInElectronKinetics)
                  eedfAbsDensities(j) = eedfAbsDensities(j) + eedfAbsDensities(k);
                end
              end
              eedfAbsDensities(i) = eedfAbsDensities(i) + eedfAbsDensities(j);
            end
          end
        end
        % update densities of states in the the electron kinetics gas mixture (normalized to gas density)
        for i = chemistry.statesIDsToUpdateInElectronKinetics
          chemistry.stateArray(i).eedfEquivalent.density = eedfAbsDensities(i)/eedfTotalGasDensity;
        end
        % renormalize electron kinetics gas mixture
        for i = chemistry.gasesIDsToUpdateInElectronKinetics
          chemistry.gasArray(i).eedfEquivalent.renormalizeWithDensities();
          % check for the distribution of states to be properly normalised
          chemistry.gasArray(i).eedfEquivalent.checkPopulationNorms();
        end
        % update the density dependencies of the electron kinetics (this also updates gas temperature dependencies)
        electronKinetics.updateDensityDependencies();
        % evaluate eedf temporal derivatives with new densities
        %       chemistry.electronKinetics.solve();
        eedfDerivatives = eedfTimeDerivative(time, [eedf; electronDensity], electronKinetics, false, false);
        eedfDerivatives = eedfDerivatives(1:end-1);
        % evaluate new electron impact rate coefficients and transport parameters
        electronKinetics.eedf = eedf';
        electronKinetics.evaluateMacroscopicParameters;
        chemistry.electronTransportProperties.reducedDiffCoeff = electronKinetics.swarmParam.redDiffCoeff;
        chemistry.electronTransportProperties.reducedMobility = electronKinetics.swarmParam.redMobility;
      elseif chemistry.includeThermalModel
        workCond.electronTemperature = gasTemperature*Constant.boltzmann/Constant.electronCharge;
      end
  end
  
  % evaluate time dependent rate coefficients
  KbTgInEV = Constant.boltzmannInEV*gasTemperature;  % Copy thermal energy in eV
  for i = 1:length(timeDependentReactionIDs)
    reactionID = timeDependentReactionIDs(i);
    % for flow model 'ensureIsobaric' the rate coefficients are evaluated after calculating the temperature derivative
    if (chemistry.ensureIsobaric && ~strcmp(reactionArray(reactionID).type,'outFlow')) || ~chemistry.ensureIsobaric 
      directRateCoeffs(reactionID) = rateCoeffFuncHandles{reactionID}(time, densities, totalGasDensity, ...
        reactionArray(reactionID), rateCoeffParams{reactionID}, chemistry);
      if isTimeDependentReactionReverse(i)
        inverseRateCoeffs(reactionID) = detailedBalance(reactionArray(reactionID), directRateCoeffs(reactionID), KbTgInEV);
      end  
    end
  end

  % evaluate rates of the reactions
  densitiesMatrix = repmat(densities',[numberOfReactions 1]);
  drateaux(lstdsteq1)=densitiesMatrix(lstdsteq1);
  drateaux(lstdstgt1)=densitiesMatrix(lstdstgt1).^totalReactantStoiCoeffs(lstdstgt1);
  irateaux(lstisteq1)=densitiesMatrix(lstisteq1);
  irateaux(lstistgt1)=densitiesMatrix(lstistgt1).^totalProductStoiCoeffs(lstistgt1);
  reactionRates = directRateCoeffs'.*(electronDensitiesReactant.*prod(drateaux,2)) - ...
    inverseRateCoeffs'.*(electronDensitiesProduct.*prod(irateaux,2));
  for i = chemistry.gasStabilisedReactionIDs
    reactionRates(i) = reactionRates(i)*totalGasDensity;
  end
  
  % evaluate temporal derivatives for the densities
  densityDerivatives = productMinusReactantStoiCoeffs*reactionRates;
  
  % evaluate temporal derivative of the gas temperature (if thermal model is activated) --- 2nd part of the calculations
  if chemistry.includeThermalModel
    volumeReactionIDs = chemistry.volumeReactionIDs;
    surfaceReactionIDs = chemistry.surfaceReactionIDs;
    transportReactionIDs = chemistry.transportReactionIDs;
    fw = chemistry.thermalModelWallFraction;
    % evaluation of source terms
    if includeJouleHeating
      elasticCollisions = 0;
      wallSource = 0;     
      if ~isempty(electronKinetics) 
        volumeSource = workCond.electronDensity*totalGasDensity*electronKinetics.power.field;
      elseif chemistry.isPulsed
        volumeSource = workCond.electronDensity*chemistry.electronTransportProperties.reducedMobility* ...
          workCond.reducedFieldSI^2*totalGasDensity;
      end
    else  
      volumeSource = dot(reactionRates(volumeReactionIDs), reactionEnthalpies(volumeReactionIDs)); 
      wallSource = dot(reactionRates([transportReactionIDs surfaceReactionIDs]), ...
        reactionEnthalpies([transportReactionIDs surfaceReactionIDs]));
      if ~isempty(electronKinetics)
        elasticCollisions = -workCond.electronDensity*totalGasDensity*electronKinetics.power.elasticNet;      
        % consider additional sources of heating directly from electron kinetics (if needed)
        for i = 1:numberOfGases
          if includeGasRotHeating(i) || includeGasVibHeating(i)  
            gasName = chemistry.gasArray(i).name;
            gasPower = electronKinetics.power.gases.(gasName);
          end  
          if includeGasRotHeating(i)
            volumeSource = volumeSource - workCond.electronDensity*totalGasDensity*gasPower.rotationalNet;
          end
          if includeGasVibHeating(i)
            volumeSource = volumeSource - workCond.electronDensity*totalGasDensity*gasPower.vibrationalNet;
          end
        end
        % include rotational heating from the CAR operator of the electron kinetics (if applicable)
        volumeSource = volumeSource - workCond.electronDensity*totalGasDensity*electronKinetics.power.carNet;
      elseif chemistry.isPulsed
        elasticCollisions = workCond.electronDensity*totalGasDensity* ...
          interp1(lookUpXValues, chemistry.lookUpTablePowerElasticNetValues, lookUpQueryPoint);
        % consider additional sources of heating directly from power look-up tables (if needed)
        if ~rotStatesPresent
          volumeSource = volumeSource - workCond.electronDensity*totalGasDensity* ...
            (interp1(lookUpXValues, chemistry.lookUpTablePowerRotationalNetValues, lookUpQueryPoint) + ...
            interp1(lookUpXValues, chemistry.lookUpTablePowerCARNetValues, lookUpQueryPoint));
        end
        if ~vibStatesPresent
          volumeSource = volumeSource - workCond.electronDensity*totalGasDensity* ...
            interp1(lookUpXValues, chemistry.lookUpTablePowerVibrationalNetValues, lookUpQueryPoint);
        end
      else
        elasticCollisions = 0;
      end
    end   
    % update near wall temperature and wall temperature in case of 'external' boundary condition
    if strcmp(chemistry.thermalModelBoundary, 'external')
        aux = 4*thermalConductivity/workCond.chamberRadius;
        hint = chemistry.intConvCoeff;
        hext = chemistry.extConvCoeff;
        wallTemperature = (aux*gasTemperature/hext + (aux/hint+1)*workCond.extTemperature + ...
          (aux/hint/hext+1/hext)*wallSource*workCond.chamberRadius*0.5) / (aux/hint+aux/hext+1);
        nearWallTemperature = ((aux/hint+aux/hext)*gasTemperature + workCond.extTemperature + ...
          wallSource*workCond.chamberRadius*0.5/hext) / (aux/hint+aux/hext+1);
        workCond.wallTemperature = wallTemperature;
        workCond.nearWallTemperature = nearWallTemperature;
    end
    conduction = -thermalConductivity*(gasTemperature-nearWallTemperature)/charLengthThermalModelSquared;
    temperatureDerivative = (conduction + elasticCollisions + volumeSource + fw*wallSource)/heatCapacity;
    thermalModel = struct('nearWallTemperature', nearWallTemperature, 'wallTemperature', wallTemperature, ...
      'conduction', conduction, 'elasticCollisions', elasticCollisions, 'volumeSource', volumeSource, ...
      'wallSource', wallSource);
  else
    temperatureDerivative = 0;
    thermalModel = struct.empty;
  end
  
  % calculate the ouflow for the 'ensureIsobaric' model (constant pressure)
  if chemistry.ensureIsobaric

    % evaluate the current value of the outflow (in sccm) and update it in the working conditions
    totalDensDerivativeExceptOutFlow = sum(densityDerivatives(childlessVolumePhaseSpeciesIDs),'all');
    temperatureLogarithmicDerivative = temperatureDerivative/gasTemperature;
    workCond.totalSccmOutFlow = ...
        (totalDensDerivativeExceptOutFlow + temperatureLogarithmicDerivative*totalGasDensity) * ...
        chamberVolume/Constant.sccmToParticleRateCoeff;

    % evaluate the current value of the outflow rate coefficient (in s-1)
    directRateCoeffFlow = rateCoeffFuncHandles{firstReactionFlowID}(time, densities, totalGasDensity, ...
        reactionArray(firstReactionFlowID), rateCoeffParams{firstReactionFlowID}, chemistry);

    % evaluate the time dependent outflow rate (in m-3 s-1)
    densityDerivatives(childlessVolumePhaseSpeciesIDs) = densityDerivatives(childlessVolumePhaseSpeciesIDs) - ...
      directRateCoeffFlow*densities(childlessVolumePhaseSpeciesIDs);  
  end

  % combine derivatives of the diferent variables
  derivatives = [densityDerivatives; temperatureDerivative];
  if strcmp(chemistry.electronKineticsDependence, 'postDischarge') && chemistry.solveEedf
    derivatives = [derivatives; eedfDerivatives];
  elseif strcmp(chemistry.electronKineticsDependence, 'quasiStationary') && strcmp(chemistry.lookUpMethod, 'localEnergy')
    powerField = chemistry.electronTransportProperties.reducedMobility*totalGasDensity*workCond.electricFieldSI^2;
    powerCollisions = totalGasDensity * interp1(lookUpXValues, chemistry.lookUpTablePowerFieldValues, lookUpQueryPoint);
    TeDerivative = 2/3*(powerField + powerCollisions);
    derivatives = [derivatives; TeDerivative];
  end
  
end

function inverseRateCoeff = detailedBalance(reaction, directRateCoeff, KbTgInEV)
% detailedBalance evaluate the rate coefficient of the inverse reaction by taking into account the detailed
% balance. For the moment, the detailed balance is only implemented for reactions of type "eedf", i.e. rate
% coefficients calculated through the integration of a cross section with an EEDF, and for two body reactions
% between heavy species (not electrons).
  
  persistent statWeightRatio;
  persistent exponent;
  
  % in case the inverse rate coefficient is already evaluated by the rate coefficient function return that value
  if ~isempty(reaction.backRateCoeff)
    inverseRateCoeff = reaction.backRateCoeff;
    return;
  end

  % in case the inverse rate coefficient has not been calculated check if it falls into one of the "standard" categories
  switch reaction.type
    case 'eedf'
      inverseRateCoeff = reaction.eedfEquivalent.supRateCoeff;
    otherwise
      reactionID = reaction.ID;
      if length(statWeightRatio)<reactionID || statWeightRatio(reactionID) == 0
        if reaction.reactantElectrons == 0 && reaction.productElectrons == 0 && ...
            sum(reaction.reactantStoiCoeff)+sum(reaction.catalystStoiCoeff) == 2 && ...
            sum(reaction.productStoiCoeff)+sum(reaction.catalystStoiCoeff) == 2
          statWeightRatio(reactionID) = 1;
          exponent(reactionID) = 0;
          for i = 1:length(reaction.reactantArray)
            state = reaction.reactantArray(i);
            if isempty(state.statisticalWeight)
              error('Unable to find %s statatistical weight to evaluate the detail balance of reaction: \n%s', ...
                state.name, reaction.description);
            elseif isempty(state.energy)
              error('Unable to find %s energy to evaluate the detail balance of reaction: \n%s', ...
                state.name, reaction.description);
            end
            statWeightRatio(reactionID) = statWeightRatio(reactionID) * ...
              state.statisticalWeight^reaction.reactantStoiCoeff(i);
            exponent(reactionID) = exponent(reactionID) - ...
              state.energy*reaction.reactantStoiCoeff(i);
          end
          for i = 1:length(reaction.productArray)
            state = reaction.productArray(i);
            if isempty(state.statisticalWeight)
              error('Unable to find %s statatistical weight to evaluate the detail balance of reaction: \n%s', ...
                state.name, reaction.description);
            elseif isempty(state.energy)
              error('Unable to find %s energy to evaluate the detail balance of reaction: \n%s', ...
                state.name, reaction.description);
            end
            statWeightRatio(reactionID) = statWeightRatio(reactionID) / ...
              state.statisticalWeight^reaction.productStoiCoeff(i);
            exponent(reactionID) = exponent(reactionID) + ...
              state.energy*reaction.productStoiCoeff(i);
          end        
        else
          error(['Error found when evaluating the inverse rate coefficient for reaction:\n%s\n' ...
            'Detailed balance is not implemented for this reaction.'], reaction.description);
        end
      end
      inverseRateCoeff = directRateCoeff*statWeightRatio(reactionID)*exp(exponent(reactionID)/KbTgInEV);
  end
  
end

function newParameter = iterateOverParameter(iterationIDs, relErrors, parameters, limitLinearExtrapolation, ...
  bissectionToBeDone)
% iterateOverParameter evaluates the new value of a certain parameter that it is iterated over based on data from 
% previous iterations (the IDs of the iterations, the errors obtained and the values of the parameter for those
% iterations)
% The following boleean parameters are intended to optimize the iteration algorithm:
%  - limitLinearExtrapolation: by default, it is set 'true' in all iterative cycles. One can consider setting this 
%                              parameter 'false' in the neutralityCycle, but this requires discretion
%  - bissectionToBeDone: by default, it is set 'true' in the pressureCycle and the neutralityCycle, and 'false' in the 
%                        other cycles.

  if length(parameters)>5
    iterationIDs = iterationIDs(end-4:end);
    relErrors = relErrors(end-4:end);
    parameters = parameters(end-4:end);
  end

  
  if  length(iterationIDs)==1 || ~( any(relErrors>0) && any(relErrors<0))
    if relErrors(end) <= -1 % important for neutralityCycle in electronegative gases (net charge density can be negative for low E/N)
      newParameter = 2*parameters(end);
    elseif limitLinearExtrapolation
      limitedRelError = max([min([relErrors(end) 1]) -0.5]);
      newParameter = parameters(end)/(1+limitedRelError);  
    else
      newParameter = parameters(end)/(1+relErrors(end));  
    end
  else
    % order arrays of parameters and relative errors in ascending order of parameter
    [auxParameters, indecesParameter] = sort(parameters, 'ascend');
    auxRelErrors = relErrors(indecesParameter); 
    
    % find smallest negative and positive relative errors
    posErr = inf;
    negErr = -inf;
    for i = 1:length(auxRelErrors)
      if auxRelErrors(i) > 0 && auxRelErrors(i) < posErr
        posErr = auxRelErrors(i);
        positiveErrorID = i;
      elseif auxRelErrors(i) < 0 && auxRelErrors(i) > negErr
        negErr = auxRelErrors(i);
        negativeErrorID = i;
      end
    end
    
    % bisection method when the errors are too big
    if (posErr > 1 || negErr < -1) && bissectionToBeDone
      newParameter = (auxParameters(positiveErrorID)+auxParameters(negativeErrorID))/2;
      % avoid any parameter already used
      if any(auxParameters==newParameter)
        newParameter = newParameter*(rand*1.5+0.5);
      end
      return;
    end
    
    % select range of parameters and relative errors to be used in the interpolation function
    switch positiveErrorID
      case negativeErrorID-1
        minID = positiveErrorID;
        error = posErr;
        for i = positiveErrorID-1:-1:1
          if auxRelErrors(i)>error
            minID = i;
            error = auxRelErrors(i);
          else
            break;
          end
        end
        maxID = negativeErrorID;
        error = negErr;
        for i = negativeErrorID+1:length(auxRelErrors)
          if auxRelErrors(i)<error
            maxID = i;
            error = auxRelErrors(i);
          else
            break;
          end
        end
      case negativeErrorID+1
        maxID = positiveErrorID;
        error = posErr;
        for i = positiveErrorID+1:length(auxRelErrors)
          if auxRelErrors(i)>error
            maxID = i;
            error = auxRelErrors(i);
          else
            break;
          end
        end
        minID = negativeErrorID;
        error = negErr;
        for i = negativeErrorID-1:-1:1
          if auxRelErrors(i)<error
            minID = i;
            error = auxRelErrors(i);
          else
            break;
          end
        end
      otherwise
        minID = min([positiveErrorID negativeErrorID]);
        maxID = max([positiveErrorID negativeErrorID]);
        while maxID-minID > 1
          if sign(auxRelErrors(minID)) == sign(auxRelErrors(minID+1))
            minID = minID+1;
          end
          if sign(auxRelErrors(maxID)) == sign(auxRelErrors(maxID-1))
            maxID = maxID-1;
          end
        end
    end
    
    % interpolate new value for the parameter to ensure zero relative error
    if maxID-minID > 2
      newParameter = interp1(auxRelErrors(minID:maxID), auxParameters(minID:maxID), 0, 'spline');
      if newParameter > max([auxParameters(positiveErrorID) auxParameters(negativeErrorID)]) || ...
          newParameter < min([auxParameters(positiveErrorID) auxParameters(negativeErrorID)])
        newParameter = interp1(auxRelErrors(minID:maxID), auxParameters(minID:maxID), 0, 'linear');
      end
    else
      newParameter = interp1(auxRelErrors(minID:maxID), auxParameters(minID:maxID), 0, 'linear');
    end
    
    % avoid any parameter already used
    if any(auxParameters==newParameter)
      newParameter = newParameter*(rand*1.5+0.5);
    end
  end
  
end

function [value, isTerminal, direction] = odeEventFunction(t, variables, varargin)
% odeEventFunction is a function that evaluates the event function for the ODE solver. It is used to stop the
% integration when a certain condition is met

  persistent dischargeTime;
  
  if isempty(dischargeTime)
    dischargeTime = varargin{5}.odeDischargeTime;
  end

  value = t-dischargeTime; % stops when the integration time reaches the discharge time defined by the user
  isTerminal = 1; % stop the integration
  direction = 0; % any direction

end

function status = odeProgressBar(t,variables,flag, varargin)

  persistent progressFigure1;
  persistent progressFigure2;
  persistent progressGraph1;
  persistent progressGraph2;
  persistent integrationTimeStr1;
  persistent integrationTimeStr2;
  persistent initialClock;
  persistent initialTime;
  persistent finalTime;
  persistent allt;
  persistent allne;
  
  switch(flag)
    case 'init'
      initialTime = t(1)-1e3;
      finalTime = t(2)-1e3;
      progressFigure1 = figure('Name', 'Chemistry solver debugging window', 'NumberTitle', 'off', 'MenuBar', 'none');
      progressFigure2 = figure('Name', 'Chemistry solver debugging window', 'NumberTitle', 'off', 'MenuBar', 'none');
      integrationTimeStr1 = uicontrol('Parent', progressFigure1, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.1 0.9 0.5 0.05], 'HorizontalAlignment', 'left', 'String', 'Computational time: 0 s');
      integrationTimeStr2 = uicontrol('Parent', progressFigure2, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.1 0.9 0.5 0.05], 'HorizontalAlignment', 'left', 'String', 'Computational time: 0 s');
      progressGraph1 = axes('Parent', progressFigure1, 'Units', 'normalized', 'OuterPosition', [0 0 1 0.9], ...
        'Box', 'on', 'XScale', 'log');
      progressGraph2 = axes('Parent', progressFigure2, 'Units', 'normalized', 'OuterPosition', [0 0 1 0.9], ...
        'Box', 'on');
      initialClock = clock;
      allt = initialTime;
      allne = variables(8)+variables(9)+variables(23)+variables(24)+variables(25);
    case 'done'
      close(progressFigure1)
      vars = whos;
      vars = vars([vars.persistent]);
      varName = {vars.name};
      clear(varName{:});
    otherwise
      allt(end+1) = t-1e3;
      allne(end+1) = variables(8)+variables(9)+variables(23)+variables(24)+variables(25);
      loglog(progressGraph1, allt, allne)
      legend({'n_e'});
      semilogy(progressGraph2, 1:1000, variables(end-999:end))
      legend({'eedf'});
      integrationTimeStr1.String = sprintf('Computational time: %.1f s\nPhysical time: %e s', ...
        etime(clock, initialClock), t-1e3);
      integrationTimeStr2.String = sprintf('Computational time: %.1f s\nPhysical time: %e s', ...
        etime(clock, initialClock), t-1e3);
  end
  status = 0;
  drawnow;

end
