% LoKI-B solves a time and space independent form of the two-term 
% electron Boltzmann equation (EBE), for non-magnetised non-equilibrium 
% low-temperature plasmas excited by DC/HF electric fields from 
% different gases or gas mixtures.
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
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

classdef Boltzmann < handle
  %Boltzmann Class that solves the boltzmann equation under certain conditions to
  %obtain an EEDF
  %
  %   Boltzmann creates the different boltzmann matrices corresponding
  %   to the different elements of the boltzmann equation. In order to
  %   create a BoltzmannEq object, a gasArray with the gas mixture
  %   information and an energyGrid must be sent to the constructor.
  
  properties
    
    gasArray = Gas.empty;                   % handle array to the gas mixture
    energyGrid = Grid.empty;                % handle to the energy grid in which the Boltzmann equation is to be solved
    workCond = WorkingConditions.empty;     % handle to the working conditions of the simulation
    CARgases = [];                          % gases described by the continuous approximation for rotations (CAR)
    
    ionCollOpType;            % describes the type of ionization collision operator
    eDensGrowthModel;         % describes the electron density growth model
    includeEECollisions;      % indicates if electron-electrons collisions are accounted for
    
    nonLinearAlgorithm;       % algorithm used to solve the non-linear operators (ionization, attachment or e-e)
    mixingParameter;          % solution mixing fraction used on the iterative scheme of the growth routines
    maxEedfRelError;          % maximum relative difference for the eedf between two consecutive iterations (non-linear routines)
    maxPowerBalanceRelError;  % maximum value for the relative power balance (threshold for the warning message)
    odeOptions;               % user defined options for the ODE solver (used when nonLinearAlgorithm is temporalIntegration)
    
    isTimeDependent = false;        % false when solving quasi-stationary Boltzmann equation, true for time-dependent pulsed Boltzmann equation
    pulseFunction;                  % function handle to the reduced electric field pulse function
    pulseFirstStep = [];            % first time step for the pulse (only for the sampling of the solution)
    pulseFinalTime = [];            % final time for the pulse
    pulseSamplingType = [];         % type of sampling for time (either 'linspace' or 'logspace')
    pulseSamplingPoints = [];       % number of sampling points for the tinal solution
    pulseFunctionParameters = {};   % additional parameters to be sent to the reduced electric field pulse function
    eDensIsTimeDependent = false;   % true if electron density should evolve in pulsed simulations

    totalCrossSection = [];   % total momentum transfer cross section
    elasticCrossSection = []; % total elastic cross section
    
    g_E = [];                   % elements of the field operator
    g_c = [];                   % elements of the elastic operator
    g_car = [];                 % elements of the CAR operator
    g_fieldSpatialGrowth = [];  % elements of the field operator with spatial electron density growth model
    g_fieldTemporalGrowth = []; % elements of the field operator with temporal electron density growth model
    Afinal = [];                % elements of the energy upflux from electron-electron collisions (used in power balance)
    Bfinal = [];                % elements of the energy downflux from electron-electron collisions (used in power balance)
    
    ionizationThresholdIsSmallerThanUmax;     % indicates if at least one ionization threshold is within the energy grid range
    includeNonConservativeIonization = false; % boolean that indicates if growth models are to be used because of ionization
    includeNonConservativeAttachment = false; % boolean that indicates if growth models are to be used because of attachment
    attachmentThresholdIsSmallerThanUmax;     % indicates if at least one attachment threshold is within the energy grid range
    
    CIEff;                    % effective ionization rate coefficient
    alphaRedEff;              % reduced first Townsend ionization coefficient
    alphaee;                  % alpha constant of electron-electron collisions
    
    fieldMatrix = [];           % matrix of the electric field operator of the boltzmann equation (continuous term)
    fieldMatrixTempGrowth = []; % matrix of the electric field operator of the boltzmann equation with temporal growth model
    fieldMatrixSpatGrowth = []; % matrix of another part of the electric field operator of the boltzmann equation with spatial growth model
    elasticMatrix = [];         % matrix of the elastic collision operator of the boltzmann equation (continuous term)
    CARMatrix = [];             % matrix of the CAR operator of the boltzmann equation (continuous term)
    continuousMatrix = [];      % matrix of the continuous operator of the boltzmann equation (sum of all continuos terms)
    inelasticMatrix = [];       % matrix of the discrete operator of the boltzmann equation (sum of all discrete terms)
    ionizationConservativeMatrix = []; % matrix of the conservative ionization operator of the boltzmann equation (discrete terms)
    ionizationMatrix = [];      % matrix of the non-conservative ionization operator of the boltzmann equation (discrete terms)
    attachmentConservativeMatrix = []; % matrix of the conservative attachment operator of the boltzmann equation (discrete terms)
    attachmentMatrix = [];      % matrix of the non-conservative attachment operator of the boltzmann equation (discrete terms)
    ionTemporalGrowth= [];      % matrix of the electron-density temporal growth term
    ionSpatialGrowthD = [];     % matrix of the electron-density spatial growth diffusion term
    ionSpatialGrowthU = [];     % matrix of the electron-density spatial growth mobility term
    Aee = [];                   % auxiliary matrix to calculate upflux vector of electron-electron collisions
    Bee = [];                   % auxiliary matrix to calculate downflux vector of electron-electron collisions
    
    eedf = [];                    % eedf obtained as solution to the boltzmann equation
    power = struct.empty;         % power balance of the boltzmann equation
    swarmParam = struct.empty;    % swarm parameters obtained with the eedf
    rateCoeffAll = struct.empty;  % rate coefficients obtained with the eedf and collisions in gasArray
    rateCoeffExtra = struct.empty;% extra rate coefficients for collisions not taken into account to obtain the eedf
    firstAnisotropy = [];         % value of the first anisotropy of the EDF (two term approximation)
    
  end
  
  events
    genericStatusMessage;
    obtainedNewEedf;
  end
   
  methods (Access = public)
    
    function boltzmann = Boltzmann(setup)
      
      % store the gas array
      boltzmann.gasArray = setup.electronKineticsGasArray;
      
      % store the energy grid and add corresponding listener
      boltzmann.energyGrid = setup.energyGrid;
      addlistener(boltzmann.energyGrid, 'updatedMaxEnergy2', @boltzmann.evaluateMatrix);
      
      % store working conditions and add corresponding listeners
      boltzmann.workCond = setup.workCond;
      addlistener(boltzmann.workCond, 'updatedGasTemperature', @boltzmann.evaluateMatrix);
      %       addlistener(workCond, 'updatedGasDensity', @boltzmann.);
      addlistener(boltzmann.workCond, 'updatedReducedField', @boltzmann.evaluateFieldOperator);
      addlistener(boltzmann.workCond, 'updatedReducedExcitationFrequency', @boltzmann.evaluateFieldOperator);
      
      % store gases for which the CAR is activated (in case there is any)
      if isfield(setup.info.electronKinetics, 'CARgases')
        boltzmann.CARgases = setup.info.electronKinetics.CARgases;
      end
      
      % store electron-impact ionization operator type
      boltzmann.ionCollOpType = setup.info.electronKinetics.ionizationOperatorType;
      
      % store electron density growth model type
      boltzmann.eDensGrowthModel = setup.info.electronKinetics.growthModelType;
      
      % store electron-electron collisions setups
      boltzmann.includeEECollisions = setup.info.electronKinetics.includeEECollisions;
      
      % store information about the numerical details on how to solve the Boltzmann equation
      boltzmann.nonLinearAlgorithm = setup.info.electronKinetics.numerics.nonLinearRoutines.algorithm;
      boltzmann.maxEedfRelError = setup.info.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError;
      if strcmp(boltzmann.nonLinearAlgorithm, 'mixingDirectSolutions') 
        boltzmann.mixingParameter = setup.info.electronKinetics.numerics.nonLinearRoutines.mixingParameter;
      else
        % store configuration of the ODE solver
        options = odeset();
        for parameter = fields(options)'
          if isfield(setup.info.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters') && ...
              isfield(setup.info.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, parameter{1})
            options.(parameter{1}) = setup.info.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.(parameter{1});
          end
        end
        boltzmann.odeOptions = options;
      end
      boltzmann.maxPowerBalanceRelError = setup.info.electronKinetics.numerics.maxPowerBalanceRelError;
      
      % check if the the simulation is steady-state or pulsed
      if setup.pulsedSimulation
        % store information about pulsed simulation in case it is activated
        boltzmann.isTimeDependent = true;
        boltzmann.pulseFunction = setup.pulseInfo.function;
        boltzmann.pulseFirstStep = setup.pulseInfo.firstStep;
        boltzmann.pulseFinalTime = setup.pulseInfo.finalTime;
        boltzmann.pulseSamplingType = setup.pulseInfo.samplingType;
        boltzmann.pulseSamplingPoints = setup.pulseInfo.samplingPoints;
        boltzmann.pulseFunctionParameters = setup.pulseInfo.functionParameters;
        % set initial value of the reduced electric field in the working conditions object
        boltzmann.workCond.reducedField = boltzmann.pulseFunction(0, boltzmann.pulseFunctionParameters);
        boltzmann.workCond.reducedFieldSI = boltzmann.workCond.reducedField*1e-21;
        
      end
      
      % allocate memory for different properties
%       boltzmann.totalCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.elasticCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.fieldMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.elasticMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
      boltzmann.CARMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.g_car = zeros(1,boltzmann.energyGrid.cellNumber+1);
%       boltzmann.continuousMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.discreteMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.ionizationIneMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
      
      % initialize quantities of the ionization and electron-electron routines
%       boltzmann.ionizationMatrix = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.Aee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.Bee = zeros(boltzmann.energyGrid.cellNumber,boltzmann.energyGrid.cellNumber);
%       boltzmann.ionSpatialGrowthD  = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.ionSpatialGrowthU = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixSpatGrowth = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.fieldMatrixTempGrowth = zeros(boltzmann.energyGrid.cellNumber);
%       boltzmann.ionTemporalGrowth = zeros(boltzmann.energyGrid.cellNumber);
      
      %evaluate boltzmann matrix
      boltzmann.evaluateMatrix();
      
    end
    
    function solve(boltzmann)
      % solve solves the electron Boltzmann equation according to the setup specified by the user
      
      % logging start of the boltzmann calculations
      start = tic;
      notify(boltzmann, 'genericStatusMessage', StatusEventData('\t- Solving Boltzmann ...\n', 'status'));
      
      % select the time-dependent or time-independent solutions
      if boltzmann.isTimeDependent
        boltzmann.obtainTimeDependentSolution();
      else
        boltzmann.obtainTimeIndependentSolution();
      end

      % logging end of the boltzmann calculations
      str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
      notify(boltzmann, 'genericStatusMessage', StatusEventData(str, 'status'));
      
    end

    function evaluateMacroscopicParameters(boltzmann)
    % evaluateMacroscopicParameters evaluate electron-impact rate coefficients, swarm parameters and power channels
    % using Boltzmann operators from the Boltzmann equation 

      boltzmann.evaluatePower(false);
      boltzmann.evaluateSwarmParameters;
      boltzmann.evaluateRateCoeff;

    end
    
    function updateDensityDependencies(boltzmann)
    % updateDensityDependencies is a function that evaluates all the species density dependencies of the boltzmann
    % object.
      
      boltzmann.evaluateMatrix();
      
    end
    
  end
  
  methods (Access = private)
    
    function updateGasTemperatureDependencies(boltzmann, ~, ~)
    % updateGasTemperatureDependencies is a function that evaluates all the dependencies of the boltzmann object (either
    % direct or indirect) on the gasTemperature property of the workCond object
    
      % evaluate population of state dependent on the gas temperature (e.g. Boltzmann distributions)
      populationsDependentOnGasTemperature = false;
      for gas = boltzmann.gasArray
        for state = gas.stateArray
          if ~isempty(state.populationFunc)
            for parameter = state.populationParams
              if strcmp(parameter{1}, 'gasTemperature')
                state.evaluatePupulation(boltzmann.workCond);
                populationsDependentOnGasTemperature = true;
                break;
              end
            end
          end
        end
      end
      
      % evaluate Boltzmann operators dependent on the gas temperature
      if populationsDependentOnGasTemperature 
        boltzmann.evaluateMatrix;
      else
        boltzmann.evaluateContinuousOperators;
      end
      
    end
    
    function evaluateMatrix(boltzmann, ~, ~)
      
      % evaluate continuous operators of the Boltzmann equation (field+elastic+CAR operators)
      boltzmann.evaluateContinuousOperators();
      
      % evaluate discrete operators of the Boltzmann equation (inelastic and superelastic operators without including
      % the contributions from ionization nor attachment)
      boltzmann.evaluateInelasticOperators();
      
      % evaluate ionization operator of the Boltzmann equation
      boltzmann.evaluateIonizationOperator();
      
      % evaluate attachment operator of the Boltzmann equation
      boltzmann.evaluateAttachmentOperator();
      
    end
    
    function evaluateContinuousOperators(boltzmann)
      
      % evaluate total and elastic momentum transfer cross sections for continuous operators
      boltzmann.evaluateTotalAndElasticCrossSections();
      
      % evaluate field operator
      boltzmann.evaluateFieldOperator();
      
      % evaluate elastic operator
      boltzmann.evaluateElasticOperator();
      
      % evaluate CAR operator
      if ~isempty(boltzmann.CARgases)
        boltzmann.evaluateCAROperator();
      end
      
    end
    
    function evaluateTotalAndElasticCrossSections(boltzmann)
      
      % reset values to zero
      boltzmann.totalCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
      boltzmann.elasticCrossSection = zeros(1,boltzmann.energyGrid.cellNumber+1);
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
        % avoid dummy gasses
        if isempty(gas.collisionArray)
          continue;
        end
        % evaluation of the mass ratio
        massRatio = Constant.electronMass/gas.mass;
        % loop over each collision with the gas
        for collision = gas.collisionArray
          % avoid effective collisions
          if strcmp(collision.type, 'Effective')
            continue;
          end
          % add collision cross section to the total momentum transfer cross section (also superelastic)
          if isempty(collision.momentumTransferCrossSection)
            boltzmann.totalCrossSection = boltzmann.totalCrossSection + ...
              collision.target.density*collision.crossSection;
          else
            boltzmann.totalCrossSection = boltzmann.totalCrossSection + ...
              collision.target.density*collision.momentumTransferCrossSection;
          end
          if collision.isReverse
            [superElasticCS, superElasticMTCS] = collision.superElasticCrossSection;
            if isempty(superElasticMTCS)
              boltzmann.totalCrossSection = boltzmann.totalCrossSection + ...
                collision.productArray.density*superElasticCS;
            else
              boltzmann.totalCrossSection = boltzmann.totalCrossSection + ...
                collision.productArray.density*superElasticMTCS;
            end
          end
          % add elastic collision cross section to the total elastic cross section (weighted by the mass ratio)
          if strcmp(collision.type, 'Elastic')
            boltzmann.elasticCrossSection = boltzmann.elasticCrossSection + ...
              massRatio*collision.target.density*collision.crossSection;
            continue;
          end
        end
      end
      
    end
    
    function evaluateFieldOperator(boltzmann, ~, ~)
      % evaluateFieldOperator is in charge of the evaluation of the elements of the field operator. The elements of
      % the operator is given by the following expression: (please note that the factor EoN^2 is not included because
      % of compatibility reasons with the simulation of field pulses, this factor is later added when doing calculations
      % with g_E and/or fieldMatrix)
      %
      %              g_E(u) = -N*sqrt(2*e/me)*EoN^2/(sigmaT(u)*(1+WoN^2/(2*e*u*sigmaT^2(u)/me)))
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      gamma = Constant.gamma;                       % gamma parameter (sqrt(2e/me))
      WoN = boltzmann.workCond.reducedExcFreqSI;    % reduced angular exitation frequency (SI units)
      cellNumber = boltzmann.energyGrid.cellNumber; % number of energy cells in the energy grid
      
      % evaluation of the elements of the operator
      boltzmann.g_E = boltzmann.energyGrid.node./(3*boltzmann.totalCrossSection.*(1+(WoN/gamma)^2./...
        (boltzmann.energyGrid.node.*boltzmann.totalCrossSection.^2)));
      
      % evaluation of the boundary conditions
      boltzmann.g_E(1) = 0;
      boltzmann.g_E(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.fieldMatrix)
        boltzmann.fieldMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      % evaluate diagonal elements
      boltzmann.fieldMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
        - (boltzmann.g_E(1:cellNumber)+boltzmann.g_E(2:cellNumber+1))./boltzmann.energyGrid.step^2;
      % evaluate inferior diagonal elements
      boltzmann.fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) = ...
        boltzmann.g_E(2:cellNumber)./boltzmann.energyGrid.step^2;
      % evaluate superior diagonal elements
      boltzmann.fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = ...
        boltzmann.g_E(2:cellNumber)./boltzmann.energyGrid.step^2;
      
    end
    
    function evaluateElasticOperator(boltzmann)
      % evaluateElasticOperator is in charge of the evaluation of the elements of the elastic collision operator. The
      % elements of the operator are given by the following expression:
      %
      %                               g_c(u) = -N*sqrt(2*e/me)*2*u^2*sigmaC(u)
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      kb = Constant.boltzmannInEV;                  % boltzmann constant in eV
      Tg = boltzmann.workCond.gasTemperature;       % gas temperature in K
      cellNumber = boltzmann.energyGrid.cellNumber; % number of energy cells in the energy grid
      factor1 = (kb*Tg/boltzmann.energyGrid.step+0.5)/boltzmann.energyGrid.step;
      factor2 = (kb*Tg/boltzmann.energyGrid.step-0.5)/boltzmann.energyGrid.step;
      
      % evaluation of the elements of the operator
      boltzmann.g_c = 2*boltzmann.energyGrid.node.^2.*boltzmann.elasticCrossSection;
      
      % evaluation of the boundary conditions
      boltzmann.g_c(1) = 0;
      boltzmann.g_c(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.elasticMatrix)
        boltzmann.elasticMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      % evaluate diagonal elements
      boltzmann.elasticMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
        - (boltzmann.g_c(1:cellNumber)*factor1+boltzmann.g_c(2:cellNumber+1)*factor2);
      % evaluate inferior diagonal elements
      boltzmann.elasticMatrix(2:cellNumber+1:cellNumber*cellNumber) = boltzmann.g_c(2:cellNumber).*factor2;
      % evaluate superior diagonal elements
      boltzmann.elasticMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = boltzmann.g_c(2:cellNumber).*factor1;
      
    end
    
    function evaluateCAROperator(boltzmann)
      % evaluateCAROperator is in charge of the evaluation of the elements of the continuous approximation for rotation
      % (CAR) operator with Chapman-Cowling correction. The operator is developed for gases well described by the
      % rotational cross sections proposed by Gerjuoy and Stein (complete development of the operator can be found in
      % "M A Ridenti, L L Alves, V Guerra and J Amorim, Plasma Sources Sci. Technol. 24 (2015) 035002 (16pp)).
      %
      % The elements of the operator are given by the following expression:
      %
      %                               g_car(u) = -N*sqrt(2*e/me)*4*u*sigma0*B
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % definition of intermediate variables
      kb = Constant.boltzmannInEV;                  % boltzmann constant in eV
      Tg = boltzmann.workCond.gasTemperature;       % gas temperature in K
      cellNumber = boltzmann.energyGrid.cellNumber; % number of energy cells in the energy grid
      factor1 = (kb*Tg/boltzmann.energyGrid.step+0.5)/boltzmann.energyGrid.step;
      factor2 = (kb*Tg/boltzmann.energyGrid.step-0.5)/boltzmann.energyGrid.step;
      
      % evaluation of the elements of the operator
      sigma0B = 0;
      for gasName = boltzmann.CARgases
        gasID = Gas.find(gasName, boltzmann.gasArray);
        gas = boltzmann.gasArray(gasID);
        sigma0B = sigma0B + gas.fraction*gas.electricQuadrupoleMoment^2*gas.rotationalConstant;
      end
      sigma0B = 8.0*pi*sigma0B/(15.0*Constant.electronCharge^2*Constant.bohrRadius^2);
      boltzmann.g_car =4*boltzmann.energyGrid.node.*sigma0B;
      
      % evaluation of the boundary conditions
      boltzmann.g_car(1) = 0;
      boltzmann.g_car(end) = 0;
      
      % loading of the matrix
      if isempty(boltzmann.CARMatrix)
        boltzmann.CARMatrix = zeros(boltzmann.energyGrid.cellNumber);
      end
      
      % evaluate diagonal elements
      boltzmann.CARMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
        - (boltzmann.g_car(1:cellNumber)*factor1+boltzmann.g_car(2:cellNumber+1)*factor2);
      % evaluate inferior diagonal elements
      boltzmann.CARMatrix(2:cellNumber+1:cellNumber*cellNumber) = boltzmann.g_car(2:cellNumber).*factor2;
      % evaluate superior diagonal elements
      boltzmann.CARMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = boltzmann.g_car(2:cellNumber).*factor1;
      
    end
    
    function evaluateInelasticOperators(boltzmann)
      % evaluateInelasticOperators is in charge of the evaluation of the elements of the inelastic and superelastic 
      % collision operator (it does not include the contribution due to ionization nor attachment collisions)
      %
      % For more information see the notes about the discretisation of the boltzmann equation.
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      
      % allocate memory for the matrix
      inelasticMatrixAux = zeros(cellNumber);
      
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          collisionType = collision.type;
          threshold = collision.threshold;
          % avoid Effective, Elastic, Ionization and attachment collisions, also avoid collisions which threshold is 
          % larger than the maximum energy or smaller than the energy step
          if strcmp(collisionType, 'Effective') || strcmp(collisionType, 'Elastic') || ...
              strcmp(collisionType, 'Ionization') || strcmp(collisionType, 'Attachment') || ...
              threshold > energyNode(end) || threshold < energyStep
            continue;
          end
          % evaluate numerical threshold
          numThreshold=floor(threshold/energyStep);
          % evaluate cross section at cell positions
          nodeCrossSection = collision.crossSection;
          cellCrossSection = 0.5*(nodeCrossSection(1:end-1)+nodeCrossSection(2:end));
          % store density of the target of the collision
          targetDensity = collision.target.density;
          % load matrix

          % evaluation of core inelastic/superelastic matrix elements
          energyTimesCrossSection = energyCell.*cellCrossSection;
          % evaluation of inleastic matrix elements
          inelasticMatrixElements = targetDensity*energyTimesCrossSection;
          % fill "exits" in the inelastic matrix (inelastic contribution)
          inelasticMatrixAux(1:cellNumber+1:cellNumber*cellNumber) = ...
            inelasticMatrixAux(1:cellNumber+1:cellNumber*cellNumber) - ...
            inelasticMatrixElements;
          % fill "entrances" in the inelastic matrix (inelastic contribution)
          inelasticMatrixAux(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) = ...
            inelasticMatrixAux(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) + ...
            inelasticMatrixElements(1+numThreshold:cellNumber);
          if collision.isReverse
            statWeightRatio = collision.target.statisticalWeight/collision.productArray.statisticalWeight;
            productDensity = collision.productArray.density;
            % evaluation of superelastic matrix elements
            superelasticMatrixElements = statWeightRatio*productDensity*energyTimesCrossSection;
            % fill "exits" in the superelastic matrix (superelastic contribution)
            inelasticMatrixAux(1:cellNumber+1:cellNumber*(cellNumber-numThreshold)) = ...
              inelasticMatrixAux(1:cellNumber+1:cellNumber*(cellNumber-numThreshold)) - ...
              superelasticMatrixElements(1+numThreshold:cellNumber);
            % fill "entrances" in the superelastic matrix (superelastic contribution)
            inelasticMatrixAux(1+numThreshold:cellNumber+1:cellNumber*(cellNumber-numThreshold)) = ...
              inelasticMatrixAux(1+numThreshold:cellNumber+1:cellNumber*(cellNumber-numThreshold)) + ...
              superelasticMatrixElements(1+numThreshold:cellNumber);
          end

        end
      end
      
      % save discrete matrix in the boltzmann object
      boltzmann.inelasticMatrix = inelasticMatrixAux;
      
    end
    
    function evaluateIonizationOperator(boltzmann)
      % evaluateIonizationOperator is in charge of the evaluation of the elements of the ionization
      % conservative or non-conservative collision operator
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      ionConservativeMatrix = zeros(cellNumber);
      ionNonConservativeMatrix = zeros(cellNumber);
      thresholdIsSmallerThanUmax = false;
      
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          threshold = collision.threshold;
          if strcmp(collision.type, 'Ionization') && threshold <= energyNode(end)
            % set flag for non-conservative algorithms
            thresholdIsSmallerThanUmax = true;
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1) + collision.crossSection(2:end));
            % local copy of target density
            density = collision.target.density;
            % evaluate numerical threhold
            numThreshold=floor(threshold/energyStep);
            % evaluate elements of the operator
            operatorElements = density.*energyCell.*cellCrossSection;
            % evaluation of nonconservative ionization collisional operator (in case it is needed)
            switch boltzmann.ionCollOpType
              case 'oneTakesAll'
                % fill "exits" in the ionNonConservativeMatrix
                ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
                  ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) - operatorElements;
                % fill "entrances" in the ionNonConservativeMatrix (scattered electrons)
                ionNonConservativeMatrix(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) = ...
                  ionNonConservativeMatrix(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) + ...
                  operatorElements(1+numThreshold:cellNumber);
                % fill "entrances" in the ionNonConservativeMatrix (secondary electrons)
                ionNonConservativeMatrix(1,:) = ionNonConservativeMatrix(1,:) + operatorElements;
              case 'equalSharing'
                % fill "exits" in the ionNonConservativeMatrix
                ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
                  ionNonConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) - operatorElements;
                % fill "entrances" in the ionNonConservativeMatrix (both scattered and secondary electrons)
                ionNonConservativeMatrix(cellNumber*(numThreshold+1)+1:2*cellNumber+1:cellNumber*cellNumber) = ...
                  ionNonConservativeMatrix(cellNumber*(numThreshold+1)+1:2*cellNumber+1:cellNumber*cellNumber) + ...
                  4*operatorElements(2+numThreshold:2:2*floor((cellNumber-numThreshold)/2)+numThreshold);
              case 'usingSDCS'
                W = gas.OPBParameter;
                if isempty(W)
                  W = threshold;
                end
                aux1 = 1./(W.*atan((energyCell-threshold)./(2*W)));
                aux2 = 1./(1+(energyCell./W).^2);
                aux3 = energyStep*operatorElements.*aux1;
                aux4 = cumsum(aux2);
                for k=1:cellNumber
                  half = floor((k-numThreshold)/2);
                  final = min([2*k+numThreshold cellNumber]);
                  % fill "exits" in the ionNonConservativeMatrix
                  if half>0
                    ionNonConservativeMatrix(k,k)= ionNonConservativeMatrix(k,k) - aux3(k)*aux4(half);
                  end
                  % fill "entrances" in the ionNonConservativeMatrix (both scattered and secondary electrons)
                  if k+numThreshold+1<=cellNumber
                    ionNonConservativeMatrix(k,k+numThreshold+1:final) = ...
                      ionNonConservativeMatrix(k,k+numThreshold+1:final) + aux3(k+numThreshold+1:final).*...
                      aux2(1:final-k-numThreshold);
                  end
                  ionNonConservativeMatrix(k,(2*k+numThreshold):cellNumber) = ...
                    ionNonConservativeMatrix(k,(2*k+numThreshold):cellNumber) + aux3(2*k+numThreshold:cellNumber).*...
                    aux2(k);
                end
            end
            % avoid writing of conservative collision operator if threshold is smaller than the energy step
            if numThreshold==0
              continue
            end
            % evaluation of conservative ionization collisional operator (always needed)
            % fill "exits" in the ionConservativeMatrix
            ionConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) = ...
              ionConservativeMatrix(1:cellNumber+1:cellNumber*cellNumber) - ...
              operatorElements;
            % fill "entrances" in the ionConservativeMatrix
            ionConservativeMatrix(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) = ...
              ionConservativeMatrix(cellNumber*numThreshold+1:cellNumber+1:cellNumber*cellNumber) + ...
              operatorElements(1+numThreshold:cellNumber);
          end
        end
      end
      
      % save matrixes in as properties of the boltzmann object
      boltzmann.ionizationConservativeMatrix = ionConservativeMatrix;
      boltzmann.ionizationMatrix = ionNonConservativeMatrix;
      boltzmann.ionizationThresholdIsSmallerThanUmax = thresholdIsSmallerThanUmax;
      
      if ~strcmp(boltzmann.ionCollOpType, 'conservative') && thresholdIsSmallerThanUmax
        boltzmann.includeNonConservativeIonization = true;
      end
      
    end
    
    function evaluateAttachmentOperator(boltzmann)
      % evaluateAttachmentOperator is in charge of the evaluation of the elements of the attachment
      % conservative or non-conservative collision operator
      
      % define local copies of variables used multiple times along the function
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      attConservativeMatrix = zeros(cellNumber);
      attNonConservativeMatrix = zeros(cellNumber);
      thresholdIsSmallerThanUmax = false;
      
      for gas = boltzmann.gasArray
        % loop over each collision with the gas
        for collision = gas.collisionArray
          threshold = collision.threshold;
          if strcmp(collision.type, 'Attachment') && threshold <= energyNode(end)
            % set flag for non-conservative algorithms
            thresholdIsSmallerThanUmax = true;
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1) + collision.crossSection(2:end));
            % local copy of target density
            density = collision.target.density;
            % evaluate numerical threshold
            numThreshold=floor(threshold/energyStep);
            
            % evaluation of nonconservative attachment collisional operator (always)
            for k=1:cellNumber
              attNonConservativeMatrix(k,k) = attNonConservativeMatrix(k,k) - ...
                density*energyCell(k)*cellCrossSection(k);
            end
            % avoid writing of conservative operator if threshold is smaller than the energy step
            if numThreshold==0
              continue
            end
            % evaluation of conservative attachment collisional operator (always needed)
            for k=1:cellNumber
              if (k<=cellNumber-numThreshold)
                attConservativeMatrix(k,k+numThreshold) = attConservativeMatrix(k,k+numThreshold) + ...
                  density*energyCell(k+numThreshold)*cellCrossSection(k+numThreshold);
              end
              attConservativeMatrix(k,k) = attConservativeMatrix(k,k) - density*energyCell(k)*cellCrossSection(k);
            end
            
          end
        end
      end
      
      % save matrixes in as properties of the boltzmann object
      boltzmann.attachmentConservativeMatrix = attConservativeMatrix;
      boltzmann.attachmentMatrix = attNonConservativeMatrix;
      boltzmann.attachmentThresholdIsSmallerThanUmax = thresholdIsSmallerThanUmax;
      
      if thresholdIsSmallerThanUmax
        boltzmann.includeNonConservativeAttachment = true;
      end
      
    end

    function obtainTimeIndependentSolution(boltzmann)
      % obtainTimeIndependentSolution solves the time-independent electron boltzmann equation 
      % (to be used for steady-state simulations) 
      
      % save appropiate method for the selected non-linear algorithm
      switch boltzmann.nonLinearAlgorithm
        case 'mixingDirectSolutions'
          nonLinearSolver = str2func('mixingDirectSolutions');
        case 'temporalIntegration'
          nonLinearSolver = str2func('temporalIntegration');
      end

      % check for the presence of non-linear operators in the Boltzmann equation
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
          boltzmann.includeEECollisions
        % iterate non-linear solver until convergence
        auxEedf = nonLinearSolver(boltzmann);
      else
        % invert the boltzmann matrix to obtain an eedf (stationary solution)
        auxEedf = boltzmann.linearSolver;
      end
      
      % when the smart grid is activated the minimum number of decades of decay in the eedf is ensured
      if boltzmann.energyGrid.isSmart
        decades = log10(auxEedf(1))-log10(auxEedf(end));
        while decades < boltzmann.energyGrid.minEedfDecay
          % increase maximum value of the energy grid
          boltzmann.energyGrid.updateMaxValue(boltzmann.energyGrid.node(end)*(1+boltzmann.energyGrid.updateFactor));
          % check for the presence of non-linear operators in the Boltzmann equation
          if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
              boltzmann.includeEECollisions
            % iterate non-linear solver until convergence
            auxEedf = nonLinearSolver(boltzmann);
          else
            % invert the boltzmann matrix to obtain an eedf (stationary solution)
            auxEedf = boltzmann.linearSolver;
          end
          % check decades of decay
          decades = log10(auxEedf(1))-log10(auxEedf(end));
        end
        while decades > boltzmann.energyGrid.maxEedfDecay
          % decrease maximum value of the energy grid
          boltzmann.energyGrid.updateMaxValue(boltzmann.energyGrid.node(end)/(1+boltzmann.energyGrid.updateFactor));
          if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
              boltzmann.includeEECollisions
            % iterate non-linear solver until convergence
            auxEedf = nonLinearSolver(boltzmann);
          else
            % invert the boltzmann matrix to obtain an eedf (stationary solution)
            auxEedf = boltzmann.linearSolver;
          end
          % check decades of decay
          decades = log10(auxEedf(1))-log10(auxEedf(end));
        end
      end
      
      % evaluate power balance
      boltzmann.evaluatePower(true);
      
      % evaluate rate coefficients
      boltzmann.evaluateRateCoeff();
      
      % evaluate transport parameters
      boltzmann.evaluateSwarmParameters();
      
      % evaluate first anisotropy 
      boltzmann.evaluateFirstAnisotropy();
      
      % bradcast obtention of a solution for the boltzmann equation
      notify(boltzmann, 'obtainedNewEedf');
      
    end
    
    function obtainTimeDependentSolution(boltzmann)
      % obtainTimeDependentSolutions solves the time-dependent electron Boltzmann equation
      % (to be used for pulsed simulations)
      
%       % obtain initial solution (steady-state solution for the initial value of the reduced electric field)
%       if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment || ...
%           boltzmann.includeEECollisions
%         % iterate non-linear solver until convergence
%         nonLinearSolver = str2func(boltzmann.nonLinearAlgorithm);
%         boltzmann.isTimeDependent = false;
%         auxEedf = nonLinearSolver(boltzmann);
%         boltzmann.isTimeDependent = true;
%       else
%         % invert the boltzmann matrix to obtain an eedf (linear solution)
%         auxEedf = boltzmann.linearSolver;
%       end
      
      % obtain initial solution (as a Maxwellian EEDF at the gas temperature)
      auxEedf = exp(-boltzmann.energyGrid.cell/(Constant.boltzmannInEV*boltzmann.workCond.gasTemperature))/ ...
        (boltzmann.energyGrid.step*dot(sqrt(boltzmann.energyGrid.cell),exp(-boltzmann.energyGrid.cell/ ...
        (Constant.boltzmannInEV*boltzmann.workCond.gasTemperature))));
      
      % switch for electron density evolution
      boltzmann.eDensIsTimeDependent = false; 
      eDensIsTimeDependentAux = boltzmann.eDensIsTimeDependent;
      
      % save initial value of electron density
      auxNe = boltzmann.workCond.electronDensity;
      
      % evaluate sampling points prescribed by the user
      if strcmp(boltzmann.pulseSamplingType, 'linspace')
        times = [0 linspace(boltzmann.pulseFirstStep, boltzmann.pulseFinalTime, boltzmann.pulseSamplingPoints)];
      elseif strcmp(boltzmann.pulseSamplingType, 'logspace')
        times = [0 logspace(log10(boltzmann.pulseFirstStep), log10(boltzmann.pulseFinalTime), ...
          boltzmann.pulseSamplingPoints)];
      end
      
      % check if the GUI is activated
      if isempty(findobj('name', 'LoKI Simulation Tool'))
        GUIisActive = false;
      else
        GUIisActive = true;
      end
      
      % initialise ODE solver parameters
      eedfTimeDerivative(0, [auxEedf auxNe]', boltzmann, true, eDensIsTimeDependentAux);    % flush persistent variables previously stored in function memory
      odeOptionsLocal = boltzmann.odeOptions;                     % load ode options specified by the user
      odeOptionsLocal.NonNegative = 1:length(auxEedf)+1;          % ensure non negative values for the components of the eedf
      if GUIisActive
        odeOptionsLocal.OutputFcn = @odeProgressBar;
        odeProgressBar(0, 0, 'firstInit', boltzmann.pulseFirstStep, boltzmann.pulseFinalTime);
      end
      
      % perform temporal integration along the pulse (dividing integration into smaller chunks for performance sake)
      initialTimeIdx = 1;
      finalTimeIdx = 1;
      finalTimes = [];
      finalEedfs = [];
      finalNe = [];
      
      while finalTimeIdx < length(times)
        % evaluate final target time of the current chunk
        if initialTimeIdx == 1
          targetFinalTime = times(2)*10;
        else
          targetFinalTime = times(initialTimeIdx)*10;
        end
        % find first element of the times array larger than the targetFinalTime
        for idx = finalTimeIdx+1:length(times)
          if times(idx) > targetFinalTime
            finalTimeIdx = idx;
            break
          elseif idx == length(times)
            finalTimeIdx = idx;
          end
        end
        % perform temporal integration of the current chunk
%         [timesSol,eedfsSol,~,~,~] = ...
        [timesSol,solutions] = ode15s(@eedfTimeDerivative,times(initialTimeIdx:finalTimeIdx), [auxEedf auxNe]', ...
          odeOptionsLocal, boltzmann, false, eDensIsTimeDependentAux);
        % separate solutions
        eedfsSol = solutions(:,1:end-1);
        neSol = solutions(:,end);
        % renormalize eedfs of the current chunck
        for idx = 1:length(timesSol)
          norm = (boltzmann.energyGrid.step*dot(eedfsSol(idx,:), sqrt(boltzmann.energyGrid.cell)));
          eedfsSol(idx,:) = eedfsSol(idx,:)/norm;
        end
        % evaluate initial conditions for the integration of the next chunk
        initialTimeIdx = finalTimeIdx;
        auxEedf = eedfsSol(end,:);
        auxNe = neSol(end);
        odeOptionsLocal.InitialStep = timesSol(end)-timesSol(end-1);
        % acumulate solution of the current chuck with the previous ones
        if isempty(finalTimes)
          finalTimes = timesSol;
          finalEedfs = eedfsSol;
          finalNe = neSol;
        else
          finalTimes = [finalTimes; timesSol(2:end)];
          finalEedfs = [finalEedfs; eedfsSol(2:end,:)];
          finalNe = [finalNe; neSol(2:end)];
        end
      end
      
      % close ode integration progress bar in case it was active
      if GUIisActive
        odeProgressBar(0,0,'lastDone');
      end
      
      % evaluate output and save solutions for each time step
      for idx = 1:length(finalTimes)
        
        % evaluate current working conditions(and save them into the working conditions object)
        boltzmann.workCond.currentTime = finalTimes(idx);
        boltzmann.workCond.reducedField = boltzmann.pulseFunction(finalTimes(idx), boltzmann.pulseFunctionParameters);
        boltzmann.workCond.reducedFieldSI = boltzmann.workCond.reducedField*1e-21;
        boltzmann.workCond.electronDensity = finalNe(idx);
        
        % save current eedf in the properties of the Boltzmann object (to be used in other methods)
        boltzmann.eedf = finalEedfs(idx,:);
        
        % evaluate intermidiate quantities needed for power evaluation
        eedfTimeDerivative(finalTimes(idx), [finalEedfs(idx,:) finalNe(idx)]', boltzmann, false, eDensIsTimeDependentAux);
        
        % evaluate power balance
        boltzmann.evaluatePower(true);
        
        % evaluate rate coefficients
        boltzmann.evaluateRateCoeff();
        
        % evaluate transport parameters
        boltzmann.evaluateSwarmParameters();
        
        % evaluate first anisotropy
        boltzmann.evaluateFirstAnisotropy();
        
        % bradcast obtention of a solution for the boltzmann equation
        notify(boltzmann, 'obtainedNewEedf');
      end
      
    end
    
    function eedf = linearSolver(boltzmann)
      % linearSolver invert the matrix of the discretized Boltzmann equation without considering non-linear
      % operators (non-conservative ionization or attachment, growth models, e-e collisions) in order to obtain an eedf
      
      % sum all contributions to the boltzmann matrix and reescale it to avoid very small numbers
      matrix = 1.e20*(boltzmann.workCond.reducedFieldSI^2*boltzmann.fieldMatrix + boltzmann.elasticMatrix + ...
        boltzmann.CARMatrix + boltzmann.inelasticMatrix + boltzmann.ionizationConservativeMatrix + ...
        boltzmann.attachmentConservativeMatrix);
      
      % invert the matrix
      eedf = matrixInversion(matrix, boltzmann.energyGrid);
      
      % store eedf in the properties of the boltzmann object
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = temporalIntegration(boltzmann)
      % eedfTimeIntegration performs a time integration of the Boltzmann equation (to be completed when the function is
      % completely developed)
      
      % initial guess for the eedf from the linear Boltzmann equation (i.e. without non-linear operators)
      eedf = boltzmann.linearSolver();
      
      % iterative ode solver
      eedfTimeDerivative(0, [eedf boltzmann.workCond.electronDensity]', boltzmann, true, false); % flush persistent variables previously stored in function memory
      odeOptionsLocal = boltzmann.odeOptions;
      odeOptionsLocal.NonNegative = 1:length(eedf);
      odeOptionsLocal.Events = @steadyStateEventFcn;
      [~,~,~,variablesSolution,~] = ode15s(@eedfTimeDerivative, [0 inf], [eedf boltzmann.workCond.electronDensity]', odeOptionsLocal, boltzmann, false, false);
      
      % select final solution if more than one is found
      if length(variablesSolution(:,1))>1
        % select eedf from final solution
        eedf = variablesSolution(end,1:end-1);
      else
        % select eedf from solution
        eedf = variablesSolution(1,1:end-1);
      end
      
      % ensure proper normalization
      norm = sum(eedf.*sqrt(boltzmann.energyGrid.cell))*boltzmann.energyGrid.step;
      eedf = eedf/norm;
      
      % store eedf in the properties of the boltzmann object
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = mixingDirectSolutions(boltzmann)
      % mixingDirectSolutions solves the non-linear Boltzmann equation (i.e. including one of more of the following
      % operators: non-conservative ionization, non-conservative attachment or e-e collision) with an algorithm of
      % mixing direct solutions of a linearized Boltzmann equation.
      
      % local copies of energy grid variables
      cellNumber = boltzmann.energyGrid.cellNumber;
      
      % initial guess for the eedf from the linear Boltzmann equation (i.e. without non-linear operators)
      boltzmann.linearSolver();
      
      % select the appropriate method depending on the growth model (if non-conservative operators are activated)
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'spatial'
            solveEGrowthModel = str2func('solveSpatialGrowthMatrix');
            % reset spatial growth related matrices
            boltzmann.ionSpatialGrowthD  = zeros(cellNumber);
            boltzmann.ionSpatialGrowthU = zeros(cellNumber);
            boltzmann.fieldMatrixSpatGrowth = zeros(cellNumber);
          case 'temporal'
            solveEGrowthModel = str2func('solveTemporalGrowthMatrix');
            % reset temporal growth related matrices
            boltzmann.fieldMatrixTempGrowth = zeros(cellNumber);
            boltzmann.ionTemporalGrowth = zeros(cellNumber);
        end
      else
        solveEGrowthModel = [];
      end
      
      % solve the non-linear Boltzmann system to obtain an eedf (different cases depending on the activated operators)
      if boltzmann.includeEECollisions
        % reset e-e collisions related variables
        boltzmann.alphaee = 0;
        boltzmann.Aee = zeros(cellNumber);
        boltzmann.Bee = zeros(cellNumber);
        boltzmann.Afinal = zeros(1, cellNumber);
        boltzmann.Bfinal = zeros(1, cellNumber);
        if isempty(solveEGrowthModel)
          % solution with the e-e collisions only
          eedf = boltzmann.solveEEColl();
        else
          % solution with e-e collisions and growth model
          globalIteration = 0;        % counter for the number of iterations of the global cycle
          maxGlobalIterations = 20;   % maximum number of iterations for the global cycle
          while globalIteration < maxGlobalIterations
            eedfOld = solveEGrowthModel(boltzmann);   % cycle updating only the growth model operator
            eedf = boltzmann.solveEEColl();           % cycle updating only the e-e collision operator
            % exit global cycle in case of convergence (maximum relative difference of each value of the eedf < 1e-9)
            if max(abs(eedf-eedfOld)./eedfOld) < boltzmann.maxEedfRelError 
              break;
            end
            globalIteration = globalIteration+1;
            % throw a warning message if the global cycle has not converged
            if globalIteration == maxGlobalIterations
              warning('Global cycle (mixing of direct solutions): EEDF did not converge after %d iterations\n', ...
                maxGlobalIterations);
            end
          end
        end
      else
        % solution with growth model only 
        eedf = solveEGrowthModel(boltzmann);
      end
      
    end
    
    function eedf = solveSpatialGrowthMatrix(boltzmann)
      
      gamma = Constant.gamma;                     % gamma parameter (sqrt(2e/me))
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      ionizationMatrixAux = boltzmann.ionizationMatrix;
      attachmentMatrixAux = boltzmann.attachmentMatrix;
      D = zeros(cellNumber);
      U = zeros(cellNumber);
      D1 = zeros(cellNumber);
      MatrixFieldSpatialGrowth = zeros(cellNumber);
      totalCrossSectionAux = boltzmann.totalCrossSection;
      cellTotalCrossSectionAux = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
      
      % writing of the Boltzmann matrix without the growth model and the electron-electron collisional operator
      baseMatrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + EoN^2*boltzmann.fieldMatrix + ...
        boltzmann.inelasticMatrix + boltzmann.ionizationMatrix + boltzmann.attachmentMatrix;
      
      % electron-electron collisions terms
      Mee = zeros(cellNumber);
      if boltzmann.includeEECollisions
        alphaEE = boltzmann.alphaee;
        auxA=boltzmann.Aee;
        auxB=boltzmann.Bee;
        A = (alphaEE/energyStep)*(auxA*eedf');
        B = (alphaEE/energyStep)*(auxB*eedf');
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
      end
      
      mixingParam = boltzmann.mixingParameter;
      
      % evaluation of the effective ionization rate integrand
      integrandCI = gamma*energyStep*sum(ionizationMatrixAux+attachmentMatrixAux)';
      
      % evaluation of the diffusion and mobility components of the spatial growth terms
      D0 = energyCell./(3*cellTotalCrossSectionAux);
      D0node = energyNode./(3*totalCrossSectionAux);
      U0sup = EoN/(2*energyStep)*[0 D0node(2:cellNumber)]; 
      U0inf = - EoN/(2*energyStep)*[D0node(2:cellNumber) 0];
      U0mid = EoN/(2*energyStep)*[-D0node(2) D0node(2:end-2)-D0node(3:end-1) D0node(end-1)];          
      U0 = U0sup+U0mid+U0inf;
          
      % evaluation of the reduced diffusion coefficient and the mobility times the electric field
      ND = gamma*energyStep*dot(D0,eedf);
      muE = -gamma*energyStep*dot(U0,eedf);

      % evaluations of the effective ionization rate
      CI = dot(eedf,integrandCI);
      
      % evaluate reduced Townsend coefficient 
      % NOTE: the initial eedf guess may lead to negative root argument, in which case, it can be calculated 
      %  ... either assuming that there is no electron density gradient: alpha/N = CI/vd = CI/(mu*E)
      %  ... or assuming that the discriminant is zero (and using the expression for the root of the quadratic equation):
      %      alpha/N = mu*E/(2D) = 2*CI/vd =  2*CI/(mu*E)  
      if muE^2-4*CI*ND < 0
        alphaRedEffNew = CI/muE;
%         alphaRedEffNew = 2*CI/muE;
      else
        alphaRedEffNew = (muE - sqrt(muE^2-4*CI*ND))/(2*ND);
      end
      
      % initialize cycle counter
      iteration = 0;

      while true

        % writing of MatrixFieldSpatialGrowth which refers to the additional electric field terms of the spatial 
        % growth model
        g_fieldSpatialGrowthAux = 0.5*alphaRedEffNew.*D0node;
        g_fieldSpatialGrowthAux(1) = 0;
        g_fieldSpatialGrowthAux(end) = 0;
        for k=1:cellNumber
          MatrixFieldSpatialGrowth(k,k) = -EoN*(g_fieldSpatialGrowthAux(k) - g_fieldSpatialGrowthAux(k+1))/energyStep;
          if k>1
            MatrixFieldSpatialGrowth(k,k-1) = -EoN*g_fieldSpatialGrowthAux(k)/energyStep;
          end
          if k<cellNumber
            MatrixFieldSpatialGrowth(k,k+1) = EoN*g_fieldSpatialGrowthAux(k+1)/energyStep;
          end
        end
        
        % calculation of the diffusion matrix of the spatial growth model
        D(1:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*alphaRedEffNew*D0;

        % calculation of the mobility matrices of the spatial growth model
        U(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*U0sup(2:cellNumber);
        U(2:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*U0inf(1:cellNumber-1);
        D1(1:cellNumber+1:cellNumber*cellNumber) = alphaRedEffNew*U0mid;
        
        % writting of the Boltzmann matrix (expansion of the following expression because of performance reasons)
        % matrixAux =  1.e20*(MatrixI + D + D1 + U + baseMatrix + Mee)
        matrixAux = 1.e20*baseMatrix;
        matrixAux(1:cellNumber+1:cellNumber*cellNumber) = matrixAux(1:cellNumber+1:cellNumber*cellNumber) + ...
          1.e20*(MatrixFieldSpatialGrowth(1:cellNumber+1:cellNumber*cellNumber) + ...
          D(1:cellNumber+1:cellNumber*cellNumber) + D1(1:cellNumber+1:cellNumber*cellNumber) + ...
        Mee(1:cellNumber+1:cellNumber*cellNumber));
        matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = ...
          matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldSpatialGrowth(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...
          U(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(2:cellNumber+1:cellNumber*cellNumber) = matrixAux(2:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldSpatialGrowth(2:cellNumber+1:cellNumber*cellNumber) + ...
          U(2:cellNumber+1:cellNumber*cellNumber) + Mee(2:cellNumber+1:cellNumber*cellNumber) );
        
        % save previous solution
        eedfOld = eedf;
        alphaRedEffOld = alphaRedEffNew;
        
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(matrixAux, boltzmann.energyGrid);
                
        % evaluation of the reduced diffusion coefficient and the mobility times the electric field
        ND = gamma*energyStep*dot(D0,eedf);
        muE = -gamma*energyStep*dot(U0,eedf);
      
        %  evaluation of the effective ionization rate
        CI = dot(eedf,integrandCI);
         
        % calculation of the new effective reduced first Townsend coefficient
        if muE^2-4*CI*ND < 0
          alphaRedEffNew = CI/muE;
%           alphaRedEffNew = 2*CI/muE;
        else
          alphaRedEffNew = (muE - sqrt(muE^2-4*CI*ND))/(2*ND);
        end
        
        % evaluate convergence criteria
        if max(abs(eedf-eedfOld)./eedfOld) < boltzmann.maxEedfRelError && ...
            (alphaRedEffNew==0 || (alphaRedEffNew-alphaRedEffOld)/alphaRedEffOld<1e-9)
          break;
        elseif iteration==400 
          if ~boltzmann.includeEECollisions
            warning('Spatial growth iterative scheme did not converge\n');
          end
          break;
        end
        
        % mixing of solutions
        alphaRedEffNew =(alphaRedEffNew*mixingParam + (1-mixingParam)*alphaRedEffOld);
        
        % update cycle counter
        iteration = iteration+1;

      end
      
      % copy of terms used in power balance and swarm parameters calculation
      boltzmann.alphaRedEff = alphaRedEffOld;
      boltzmann.g_fieldSpatialGrowth = g_fieldSpatialGrowthAux;
      boltzmann.eedf = eedf;
      
      % copy of spatial growth model terms used on the electron-electron collisions routine
      if boltzmann.includeEECollisions
        boltzmann.ionSpatialGrowthD = D+D1;
        boltzmann.ionSpatialGrowthU = U;
        boltzmann.fieldMatrixSpatGrowth = MatrixFieldSpatialGrowth;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = solveTemporalGrowthMatrix(boltzmann)
      
      gamma = Constant.gamma;                     % gamma parameter (sqrt(2e/me))
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      WoN = boltzmann.workCond.reducedExcFreqSI;  % reduced angular exitation frequency (SI units)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      ionizationMatrixAux = boltzmann.ionizationMatrix;
      attachmentMatrixAux = boltzmann.attachmentMatrix;
      totalCrossSectionAux = boltzmann.totalCrossSection;
      totalCSI = zeros(size(totalCrossSectionAux));
      g_fieldTemporalGrowthAux = zeros(size(boltzmann.g_E));
      MatrixFieldTemporalGrowth = zeros(cellNumber);
      growthMatrix = zeros(cellNumber);
      
      % writing of the Boltzmann matrix without the growth model and the electron-electron collisional operator
      baseMatrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix + ...
        ionizationMatrixAux + attachmentMatrixAux;
      
      % electron-electron collisions terms (not updated in the density growth operator)
      Mee = zeros(cellNumber);
      if boltzmann.includeEECollisions
        alphaEE = boltzmann.alphaee;
        auxA=boltzmann.Aee;
        auxB=boltzmann.Bee;
        A = (alphaEE/energyStep)*(auxA*eedf');
        B = (alphaEE/energyStep)*(auxB*eedf');
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
      end
      
      mixingParam = boltzmann.mixingParameter;
      
      % evaluation of the effective ionization rate integrand
      integrandCI = gamma*energyStep*sum(ionizationMatrixAux+attachmentMatrixAux)';
      % evaluation of the initial value for the effective ionization rate
      CIEffNew = dot(eedf,integrandCI);
      
      % initialize cycle counter
      iteration = 0;

      while true
        % writing the total cross section plus the ionization rate divided by the electron velocity
        totalCSI(1)= totalCrossSectionAux(1);
        totalCSI(2:end) =  totalCrossSectionAux(2:end) + (CIEffNew/gamma)./sqrt(energyNode(2:end));
        
        % writing of the MatrixFieldTemporalGrowth which refers to the electric field operator of the temporal growth 
        % model and growthMatrix that refers to the time variation term (dn/dt) of the temporal growth model
        g_fieldTemporalGrowthAux = energyNode./(3*totalCSI.*(1+(WoN/gamma)^2./(energyNode.*totalCSI.*totalCSI)));
        g_fieldTemporalGrowthAux(1) = 0;
        g_fieldTemporalGrowthAux(end) = 0;
        for k = 1:cellNumber
          MatrixFieldTemporalGrowth(k,k) = -EoN^2*(g_fieldTemporalGrowthAux(k)+g_fieldTemporalGrowthAux(k+1))/energyStep^2;
          if k>1
            MatrixFieldTemporalGrowth(k,k-1) = EoN^2*g_fieldTemporalGrowthAux(k)/energyStep^2;
          end
          if k<cellNumber
            MatrixFieldTemporalGrowth(k,k+1) = EoN^2*g_fieldTemporalGrowthAux(k+1)/energyStep^2;
          end
          growthMatrix(k,k) = -(CIEffNew/gamma)*sqrt(energyCell(k));
        end
        
        % writting of the Boltzmann matrix (expansion of the following expression because of performance reasons)
%         matrixAux =  1.e20*(growthMatrix + MatrixI + baseMatrix + Mee);
        matrixAux = 1.e20*baseMatrix;
        matrixAux(1:cellNumber+1:cellNumber*cellNumber) = matrixAux(1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(1:cellNumber+1:cellNumber*cellNumber) + ...
          growthMatrix(1:cellNumber+1:cellNumber*cellNumber) + Mee(1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = ...
          matrixAux(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...
          Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) );
        matrixAux(2:cellNumber+1:cellNumber*cellNumber) = matrixAux(2:cellNumber+1:cellNumber*cellNumber) + 1.e20 * ...
          ( MatrixFieldTemporalGrowth(2:cellNumber+1:cellNumber*cellNumber) + ...
          Mee(2:cellNumber+1:cellNumber*cellNumber) );
        
        % save previous solution
        eedfOld = eedf;
        CIEffOld = CIEffNew;
        
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(matrixAux, boltzmann.energyGrid);
        
        % evaluate new effective ionization rate
        CIEffNew = dot(eedf,integrandCI);
        
        % evaluate convergence criteria
        if max(abs(eedf-eedfOld)./eedfOld)<boltzmann.maxEedfRelError && ...
            (CIEffNew==0 || abs((CIEffNew-CIEffOld)/CIEffOld)<1e-9)
          break;
        elseif iteration==400
          if ~boltzmann.includeEECollisions
            warning('Temporal growth iterative scheme did not converge');
          end
          break;
        end

        % mixing of solutions
        CIEffNew = mixingParam*CIEffNew + (1-mixingParam)*CIEffOld;
        
        % update cycle counter
        iteration = iteration+1;
      end
      
      % copy of terms used in power balance and swarm parameters calculation
      boltzmann.CIEff = CIEffOld;
      boltzmann.g_fieldTemporalGrowth = g_fieldTemporalGrowthAux;
      boltzmann.eedf = eedf;
      
      % copy of temporal growth model terms used on the electron-electron collisions routine
      if boltzmann.includeEECollisions
        boltzmann.ionTemporalGrowth = growthMatrix;
        boltzmann.fieldMatrixTempGrowth = MatrixFieldTemporalGrowth;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function eedf = solveEEColl(boltzmann)
      
      e = Constant.electronCharge;                % electron charge
      e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
      ne = boltzmann.workCond.electronDensity;    % electron density (SI units)
      n0 = boltzmann.workCond.gasDensity;         % gas density (SI units)
      EoN = boltzmann.workCond.reducedFieldSI;    % reduced electric field (SI units)
      EoNTd = boltzmann.workCond.reducedField;    % reduced electric field (Td)
      
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      cellNumber = boltzmann.energyGrid.cellNumber;
      energyNode = boltzmann.energyGrid.node;
      eedf = boltzmann.eedf;
      
      Mee = zeros(cellNumber);
      
      % writting of Boltzmann matrix without the electron-electron collisional operator
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'spatial'
            baseMatrix = boltzmann.ionizationMatrix + boltzmann.attachmentMatrix + boltzmann.elasticMatrix + ...
              boltzmann.CARMatrix + boltzmann.inelasticMatrix + EoN^2*boltzmann.fieldMatrix + ...
              boltzmann.ionSpatialGrowthD + boltzmann.ionSpatialGrowthU +  boltzmann.fieldMatrixSpatGrowth;
          case 'temporal'
            baseMatrix = boltzmann.ionizationMatrix + boltzmann.attachmentMatrix + boltzmann.elasticMatrix + ...
              boltzmann.CARMatrix + boltzmann.inelasticMatrix + boltzmann.fieldMatrixTempGrowth + ...
              boltzmann.ionTemporalGrowth;
        end
      else
        baseMatrix = boltzmann.ionizationConservativeMatrix + boltzmann.attachmentConservativeMatrix + ...
          boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix + EoN^2*boltzmann.fieldMatrix ;
      end
      
      % writing auxiliary matrix A, without constant alpha
      auxA = zeros(cellNumber);
      auxMatrix = zeros(cellNumber);
      auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
      % from power conservation, terms on the last row (cellNumber,:) and on the first column (:,1) are zero
      for k=1:cellNumber-1
        auxMatrix(k,2:k) = auxEnergyArray(2:k);
        auxMatrix(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
      end
      % detailed balance condition
      for k=1:cellNumber-1
        auxA(k,2:cellNumber) = sqrt(auxMatrix(k,2:cellNumber).*auxMatrix(1:cellNumber-1,k+1)'); 
      end
      
      % writing auxiliary matrix B, without constant alpha
      auxB = transpose(auxA);

      % initial value for eedfOld (calculated without electron-electron collisional operator)
      eedfOld = eedf;
      % initialize cycle counter
      iteration = 0;
      
      while true
        % Calculation of constant alpha
        meanEnergy = energyStep*dot(energyCell.^(3/2),eedf);  % mean energy
        Te = (2/3)*meanEnergy;                                % electron temperature in eV Te = (2/3)*meanEnergy
        logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));           % Coulomb logarithm
        alpha = (ne/n0)*(e^2/(8*pi*e0^2))*logC;               % alpha

        % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
        A = (alpha/energyStep)*(auxA*eedf');
        B = (alpha/energyStep)*(auxB*eedf');
        
        % writing of the electron-electron operator
        Mee(1:cellNumber+1:cellNumber*cellNumber) = -(A(1:cellNumber)+B(1:cellNumber));
        Mee(cellNumber+1:cellNumber+1:cellNumber*cellNumber) = B(2:cellNumber);
        Mee(2:cellNumber+1:cellNumber*cellNumber) = A(1:cellNumber-1);
        
        % sum all contributions to the Boltzmann matrix and reescale it to avoid very small numbers
        matrixAux = 1e20*(baseMatrix + Mee);
        
        % invert the matrix to obtain a new solution
        eedf = matrixInversion(matrixAux, boltzmann.energyGrid);
        
        % calculation of ratio = Pee/PRef
        boltzmann.eedf = eedf;
        boltzmann.Afinal = A;
        boltzmann.Bfinal = B;
        boltzmann.evaluatePower(false);
        Preference = boltzmann.power.reference;
        Pee = boltzmann.power.electronElectronNet;
        ratio = abs(Pee/Preference);
        
        % evaluate convergence criteria
        if  max(abs(eedf-eedfOld)./eedfOld)<boltzmann.maxEedfRelError && abs(ratio)<boltzmann.maxPowerBalanceRelError
          break;
        elseif max(abs(eedf-eedfOld)./eedfOld)<boltzmann.maxEedfRelError && iteration>200
          str = sprintf('\\t- e-e iterative scheme: EEDF has converged but abs(Pee/Pref)= %.16g > %.16g\\n', ...
          ratio,boltzmann.maxPowerBalanceRelError);
          notify(boltzmann, 'genericStatusMessage', StatusEventData(str, 'status'));
          str = sprintf('\\t- Ionization degree:%f; Reduced electric field:%f (Td) \\n',ne/n0,EoNTd);
          notify(boltzmann, 'genericStatusMessage', StatusEventData(str, 'status'));
          break;
        elseif iteration == 400
          warning('Electron-electron iterative scheme: EEDF did not converge\n');
          break;
        end
        
        % update "old" values
        if iteration == 0
          eedfOld = eedf;
        else  % acceleration scheme based on the Newton-Raphson
          % eedf estimate (avoid negative values)
          eedfAux = abs(eedf - (eedf-eedfOld)*ratio/(ratio-ratioOld));
          eedfOld = eedf;
          eedf = eedfAux;
        end
        ratioOld = ratio;
        % update cycle counter
        iteration = iteration+1;
        
      end
      
      % saving auxiliary matrix to be used on the growth model routine
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        boltzmann.Aee = auxA;
        boltzmann.Bee = auxB;
        boltzmann.alphaee = alpha;
      end
      
      % saving eedf
      boltzmann.eedf = eedf;
      
    end
    
    function power = evaluatePower(boltzmann, checkPowerBalance)
      
      % initialize power structure
      power = struct('field', 0, 'elasticNet', 0, 'elasticGain', 0, 'elasticLoss', 0, 'carNet', 0, 'carGain', 0, ...
        'carLoss', 0, 'excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
        'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
        'ionizationIne', 0, 'attachmentIne', 0, 'inelastic', 0, 'superelastic', 0, 'eDensGrowth', 0, ...
        'electronElectronNet', 0, 'electronElectronGain', 0, 'electronElectronLoss', 0, 'gases', '');
      
      % save a local copy of the EEDF because of performance reasons
      eedfLocal = boltzmann.eedf;
      
      % save a local copy of the energy grid information
      energyGridLocal = boltzmann.energyGrid;   % energyGrid object
      N = energyGridLocal.cellNumber;           % number of cells in the energy grid
      energyStep = energyGridLocal.step;        % energy step
      energyCell = energyGridLocal.cell;        % value of the energy at the cells
      energyNode = energyGridLocal.node;        % value of the energy at the nodes
      
      % multiplicative constant to obtain the right units
      gamma = Constant.gamma;
      
      % auxiliary quantities needed to evaluate the elastic and CAR powers
      kTg = Constant.boltzmannInEV*boltzmann.workCond.gasTemperature;
      aux1 = kTg+energyStep*0.5;
      aux2 = kTg-energyStep*0.5;
      
      % evaluate power absorved per electron at unit gas density due to elastic collisions
      g_cLocal = boltzmann.g_c; % elements of the elastic collision operator (local copy)
      power.elasticNet = gamma*sum(eedfLocal.*(g_cLocal(2:end)*aux2-g_cLocal(1:end-1)*aux1));
      power.elasticGain = gamma*kTg*sum(eedfLocal.*(g_cLocal(2:end)-g_cLocal(1:end-1)));
      power.elasticLoss = power.elasticNet-power.elasticGain;
      
      % evaluate power absorved per electron at unit gas density due to rotations CAR
      if ~isempty(boltzmann.CARgases)
        g_carLocal = boltzmann.g_car; % elements of the CAR operator (local copy)
        power.carNet = gamma*sum(eedfLocal.*(g_carLocal(2:end)*aux2-g_carLocal(1:end-1)*aux1));
        power.carGain = gamma*kTg*sum(eedfLocal.*(g_carLocal(2:end)-g_carLocal(1:end-1)));
        power.carLoss = power.carNet-power.carGain;
      end
      
      % evaluate power gained from electric field and lost due to electron density growth terms
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'temporal'
            g_ELocal = boltzmann.workCond.reducedFieldSI^2*boltzmann.g_fieldTemporalGrowth; % elements of the electric field operator (local copy)
            power.field = gamma*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
            % evaluate power gained from the field - alternative calculation for temporal case (DC case)     
%             EoN = boltzmann.workCond.reducedFieldSI;
%             boltzmann.evaluateSwarmParameters;
%             power.field = boltzmann.swarmParam.redMobility*EoN^2;      
            % evaluate power gained from the field - alternative calculation for temporal case (HF case)     
%             EoN = boltzmann.workCond.reducedFieldSI;
%             boltzmann.evaluateSwarmParameters;
%             power.field = real(boltzmann.swarmParam.redMobilityHF)*EoN^2;      

            power.eDensGrowth = - boltzmann.CIEff*energyStep*sum(eedfLocal.*energyCell.*sqrt(energyCell));
            % evaluate power lost due to density growth - alternative calculation for temporal case      
%             boltzmann.evaluateSwarmParameters;
%             power.eDensGrowth = - boltzmann.CIEff*boltzmann.swarmParam.meanEnergy;

            case 'spatial'
            totalCrossSectionLocal = boltzmann.totalCrossSection;
            cellTotalCrossSection =  0.5*(totalCrossSectionLocal(1:N) + totalCrossSectionLocal(2:N+1));
            alphaRedEffLocal = boltzmann.alphaRedEff;
            reducedFieldSILocal = boltzmann.workCond.reducedFieldSI;
            % evaluate general elements of the electric field operator (local copy)
            g_ELocal = reducedFieldSILocal^2*boltzmann.g_E; 
            % evaluate specific elements of the electric field operator due to the spatial growth model
            g_ELocalSpatial = reducedFieldSILocal*boltzmann.g_fieldSpatialGrowth;
            % evaluate power gained from the field
            power.field = gamma*(sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1) + ...
                (-g_ELocalSpatial(2:end)-g_ELocalSpatial(1:end-1))*energyStep)));            
            % evaluate power gained from the field - alternative calculation for spatial case      
%             EoN = boltzmann.workCond.reducedFieldSI;
%             boltzmann.evaluateSwarmParameters;
%             D0node = energyNode./(3*totalCrossSectionLocal);
%             power.field = boltzmann.swarmParam.redMobility*EoN^2 - ...
%                 0.5*gamma*energyStep*alphaRedEffLocal*EoN*sum(D0node(2:end-1).*(eedfLocal(2:end)+eedfLocal(1:end-1)));
                        
            % evaluate power lost due to spatial growth: diffusion contribution
            powerDiffusion = alphaRedEffLocal^2*gamma*energyStep/3*sum(energyCell(1:N).^2.*eedfLocal(1:N)./...
              cellTotalCrossSection(1:N));            
            % evaluate power lost due to spatial growth: mobility contribution
            powerMobility = gamma*alphaRedEffLocal*(reducedFieldSILocal/3)*(...
              sum(energyNode(2:N).^2.*(eedfLocal(2:N)-eedfLocal(1:N-1))./totalCrossSectionLocal(2:N)));            
            % evaluate total power lost due to spatial growth
            power.eDensGrowth =  powerDiffusion + powerMobility;
            % evaluate power lost due to density growth - alternative calculation for spatial case      
%             EoN = boltzmann.workCond.reducedFieldSI;
%             boltzmann.evaluateSwarmParameters;
%             power.eDensGrowth = - alphaRedEffLocal*(boltzmann.swarmParam.redMobilityEnergy*EoN - ...
%                 alphaRedEffLocal*boltzmann.swarmParam.redDiffCoeffEnergy);
            
        end
      else
        % evaluate general elements of the electric field operator (local copy)
        g_ELocal = boltzmann.workCond.reducedFieldSI^2*boltzmann.g_E; 
        % evaluate power gained from the field
        power.field = gamma*sum(eedfLocal.*(g_ELocal(2:end)-g_ELocal(1:end-1)));
      end
      
      % power absorved per electron at unit gas density due to electron-electron collisions (to be revised)
      if boltzmann.includeEECollisions
        A = boltzmann.Afinal;
        B = boltzmann.Bfinal;
        power.electronElectronNet = gamma*sum((A(1:N)-B(1:N)).*eedfLocal(1:N)')*energyStep^2;
        power.electronElectronGain = gamma*0.5*( sum((A(2:N-1)+B(3:N)-A(1:N-2)-B(2:N-1)).*eedfLocal(2:N-1)')+...
          (A(1)+B(2))*eedfLocal(1)-(A(N-1)+B(N))*eedfLocal(N) )*energyStep^2;
        power.electronElectronLoss = power.electronElectronNet - power.electronElectronGain;
      end
      
      % evaluate power absorved per electron at unit gas density due to the inelastic/super-elastic collisions
      % loop over each gas in the mixture
      for gas = boltzmann.gasArray
        gasName = gas.name;
        % initialize power balance information of this gas
        gasPower = struct('excitationIne', 0, 'excitationSup', 0, 'excitationNet', 0, 'vibrationalIne', 0, ...
          'vibrationalSup', 0, 'vibrationalNet', 0, 'rotationalIne', 0, 'rotationalSup', 0, 'rotationalNet', 0, ...
          'ionizationIne', 0, 'attachmentIne', 0);
        % loop over each collision with the gas
        for collision = gas.collisionArray
          % collision type
          collType = collision.type;
          % avoid Effective or Elastic collisions and collisions which threshold is larger than the maximum energy
          if strcmp(collType, 'Effective') || strcmp(collType, 'Elastic') || collision.threshold > energyNode(end)
            continue;
          elseif strcmp(collision.type, 'Ionization') && ~strcmp(boltzmann.ionCollOpType,'conservative')
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
            
            threshold = collision.threshold;
            lmin=floor(threshold/energyStep);
            
            switch boltzmann.ionCollOpType
              case 'equalSharing'
                ionizationIneAux = -gamma*collision.target.density*energyStep*(...
                  sum(energyCell(lmin:N).^2.*cellCrossSection(lmin:N).*eedfLocal(lmin:N)) + ...
                  2*energyCell(lmin+1)*sum(energyCell(2+lmin:2:N).*cellCrossSection(2+lmin:2:N).*...
                  eedfLocal(2+lmin:2:N))-2*sum(energyCell(2+lmin:2:N).^2.*cellCrossSection(2+lmin:2:N).*...
                  eedfLocal(2+lmin:2:N)));
                
              case 'oneTakesAll'
                ionizationIneAux = -gamma*collision.target.density*energyStep*energyCell(lmin)*...
                  sum(eedfLocal(lmin:N).*energyCell(lmin:N).*cellCrossSection(lmin:N));
                
              case 'usingSDCS'
                %calculation of the total (integrated) ionization cross section using the SDCS
                TICS = zeros(size(eedfLocal));
                W = gas.OPBParameter;
                if isempty(W)
                  W = threshold;
                end
                auxArray = 1./(1+(energyCell(1:N)/W).^2);
                for k=2:N
                  auxArray(k) = auxArray(k)+auxArray(k-1);
                  kmax = floor((k-lmin)/2);
                  if kmax>0
                    TICS(k) = TICS(k) + cellCrossSection(k)/(W*atan((energyCell(k)-threshold)/(2*W)))*auxArray(kmax);
                  end
                end
                TICS = energyStep*TICS;
                ionizationIneAux = -gamma*collision.target.density*energyCell(lmin+1)*energyStep*...
                  sum(eedfLocal(1:N).*energyCell(1:N).*TICS(1:N));
            end
            gasPower.ionizationIne = gasPower.ionizationIne+ionizationIneAux;
            continue;
          elseif strcmp(collision.type, 'Attachment') && boltzmann.includeNonConservativeAttachment
            % evaluate cross section at cell positions
            cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
            threshold = collision.threshold;
            lmin=floor(threshold/energyStep);
            gasPower.attachmentIne = gasPower.attachmentIne - gamma*collision.target.density*energyStep*...
              sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).^2.*cellCrossSection(1+lmin:N));
            continue;
            
          end
          % switch to lower case because of aesthetical reasons
          collType = lower(collType);
          % evaluate cross section at cell positions
          cellCrossSection = 0.5*(collision.crossSection(1:end-1)+collision.crossSection(2:end));
          % evaluate departure cell
          lmin=floor(collision.threshold/energyStep);
          % add contribution to the power due to the inelastic collisions
          gasPower.([lower(collType) 'Ine']) = gasPower.([lower(collType) 'Ine']) - gamma*collision.target.density*...
            energyStep*energyNode(lmin+1)*sum(eedfLocal(1+lmin:N).*energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
          % add contribution to the power due to the superelastic collisions
          if collision.isReverse
            statWeightRatio = collision.target.statisticalWeight/collision.productArray.statisticalWeight;
            gasPower.([lower(collType) 'Sup']) = gasPower.([lower(collType) 'Sup']) + gamma*statWeightRatio*...
              collision.productArray.density*energyStep*energyNode(lmin+1)*sum(eedfLocal(1:N-lmin).*...
              energyCell(1+lmin:N).*cellCrossSection(1+lmin:N));
          end
        end
        % evaluate net values (for each gas)
        gasPower.excitationNet = gasPower.excitationIne + gasPower.excitationSup;
        gasPower.vibrationalNet = gasPower.vibrationalIne + gasPower.vibrationalSup;
        gasPower.rotationalNet = gasPower.rotationalIne + gasPower.rotationalSup;
        gasPower.inelastic = gasPower.excitationIne + gasPower.vibrationalIne + gasPower.rotationalIne + ...
          gasPower.ionizationIne + gasPower.attachmentIne;
        gasPower.superelastic = gasPower.excitationSup + gasPower.vibrationalSup + gasPower.rotationalSup;
        
        % evaluate net values (for the gas mixture)
        power.excitationIne = power.excitationIne + gasPower.excitationIne;
        power.excitationSup = power.excitationSup + gasPower.excitationSup;
        power.vibrationalIne = power.vibrationalIne + gasPower.vibrationalIne;
        power.vibrationalSup = power.vibrationalSup + gasPower.vibrationalSup;
        power.rotationalIne = power.rotationalIne + gasPower.rotationalIne;
        power.rotationalSup = power.rotationalSup + gasPower.rotationalSup;
        power.ionizationIne = power.ionizationIne + gasPower.ionizationIne;
        power.attachmentIne = power.attachmentIne + gasPower.attachmentIne;
        
        % store power balance information of this gas in the main structure
        power.gases.(gasName) = gasPower;
      end
      power.excitationNet = power.excitationIne + power.excitationSup;
      power.vibrationalNet = power.vibrationalIne + power.vibrationalSup;
      power.rotationalNet = power.rotationalIne + power.rotationalSup;
      power.inelastic = power.excitationIne + power.vibrationalIne + power.rotationalIne + power.ionizationIne + ...
        power.attachmentIne;
      power.superelastic = power.excitationSup + power.vibrationalSup + power.rotationalSup;
      
      % evaluate power balance
      powerValues = [power.field power.elasticGain power.elasticLoss power.carGain power.carLoss ...
        power.excitationSup power.excitationIne power.vibrationalSup power.vibrationalIne ...
        power.rotationalSup power.rotationalIne power.eDensGrowth power.electronElectronGain ...
        power.electronElectronLoss];
      totalGain = 0;
      totalLoss = 0;
      for powerValue = powerValues
        if powerValue > 0 
          totalGain = totalGain+powerValue;
        else
          totalLoss = totalLoss+powerValue;
        end
      end
      power.balance = power.field + power.elasticNet + power.carNet + power.inelastic + power.superelastic + ...
        power.eDensGrowth + power.electronElectronNet;
      power.relativeBalance = abs(power.balance)/totalGain;
      power.reference = totalGain;
      
      % store power balance information in the boltzmann properties
      boltzmann.power = power;
      
      % check for errors in the power balance for final solution
      if checkPowerBalance && power.relativeBalance > boltzmann.maxPowerBalanceRelError
        warning(sprintf(['Relative power balance greater than %e.\n' ...
          'Results may be wrong, please check input/output of the simulation'], boltzmann.maxPowerBalanceRelError));
      end
      
    end
    
    function swarmParam = evaluateSwarmParameters(boltzmann)
      
      % save local copies of different constants and variables
      gamma = Constant.gamma;
      energyNode = boltzmann.energyGrid.node;
      energyCell = boltzmann.energyGrid.cell;
      energyStep = boltzmann.energyGrid.step;
      eedfLocal = boltzmann.eedf;
      WoN = boltzmann.workCond.reducedExcFreqSI;
      
      % evaluate auxiliary total momentum transfer cross section function (see documentation)
      totalCrossSectionAux = boltzmann.totalCrossSection;
      cellTotalCrossSectionAux = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
      if strcmp(boltzmann.eDensGrowthModel,'temporal') && ...
          ( boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment )
        totalCrossSectionAux(2:end) = totalCrossSectionAux(2:end) + (boltzmann.CIEff/gamma)./sqrt(energyNode(2:end));
      end
      
      % initialize transport parameters structure
      swarmParam = struct('redDiffCoeff', [], 'redMobility', [], 'redMobilityHF', [], 'redDiffCoeffEnergy', [], ...
        'redMobilityEnergy', [], 'redTownsendCoeff', [], 'redAttCoeff', [], 'meanEnergy', [], 'characEnergy', [], ...
        'Te', [], 'driftVelocity', []);
      
      % evaluate reduced diffusion coefficient
      swarmParam.redDiffCoeff = (2*gamma/3)*energyStep*sum(energyCell.*eedfLocal./...
        (totalCrossSectionAux(1:end-1)+totalCrossSectionAux(2:end)));
      
      % evaluate reduced mobility (DC expression)
      swarmParam.redMobility = -gamma/3*sum(energyNode(2:end-1).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
        (totalCrossSectionAux(2:end-1)));

      % evaluate complex HF reduced mobility (in case the excitation frequency is not zero)
      if WoN ~= 0
        swarmParam.redMobilityHF = -gamma/3*sum(energyNode(2:end-1).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
          (totalCrossSectionAux(2:end-1)+(WoN/gamma)^2./(energyNode(2:end-1).*totalCrossSectionAux(2:end-1))));
        swarmParam.redMobilityHF = swarmParam.redMobilityHF + 1i*(1/3)*sum(sqrt(energyNode(2:end-1)).*...
          (WoN./totalCrossSectionAux(2:end-1)).*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
          (totalCrossSectionAux(2:end-1)+(WoN/gamma)^2./(energyNode(2:end-1).*totalCrossSectionAux(2:end-1))));
      end
      
      % evaluate reduced energy diffusion coefficient
      % this coefficient corresponds to the heat diffusion (G) defined in Allis 1956
      % here, and contrary to other references, the reduced energy diffusion coefficient (in eV m-1 s-1) is not normalized to the electron mean energy
      swarmParam.redDiffCoeffEnergy = 2*gamma/3*energyStep*sum(energyCell.^2.*eedfLocal./...
        (totalCrossSectionAux(1:end-1)+totalCrossSectionAux(2:end)));
    
      % evaluate reduced energy mobility
      % this coefficient corresponds to the thermoelectricity (beta) defined in Allis 1956
      % here, and contrary to other references, the reduced energy mobility (in m-1 s-1) is not normalized to the electron mean energy
      swarmParam.redMobilityEnergy = -gamma/3*sum(energyNode(2:end-1).^2.*(eedfLocal(2:end)-eedfLocal(1:end-1))./...
        totalCrossSectionAux(2:end-1));
    
      % evaluate drift velocity
      if strcmp(boltzmann.eDensGrowthModel,'spatial') && ...
          ( boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment )
        swarmParam.driftVelocity = - swarmParam.redDiffCoeff*boltzmann.alphaRedEff + ...
          swarmParam.redMobility*boltzmann.workCond.reducedFieldSI;
      else
        swarmParam.driftVelocity = swarmParam.redMobility*boltzmann.workCond.reducedFieldSI;
      end
      
      % evaluate reduced Townsend coefficient
      totalIonRateCoeff = 0;
      for gas = boltzmann.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Ionization')
            totalIonRateCoeff = totalIonRateCoeff + collision.target.density*collision.ineRateCoeff;
          end
        end
      end
      swarmParam.redTownsendCoeff = totalIonRateCoeff / swarmParam.driftVelocity;
      % evaluate reduced attachment coefficient
      totalAttRateCoeff = 0;
      for gas = boltzmann.gasArray
        for collision = gas.collisionArray
          if strcmp(collision.type, 'Attachment')
            totalAttRateCoeff = totalAttRateCoeff + collision.target.density*collision.ineRateCoeff;
          end
        end
      end
      swarmParam.redAttCoeff = totalAttRateCoeff / swarmParam.driftVelocity;
      
      % evaluate mean energy
      swarmParam.meanEnergy = sum(energyCell(1:end).^(1.5).*eedfLocal(1:end))*energyStep;
      
      % evaluate characteristic energy
      swarmParam.characEnergy = swarmParam.redDiffCoeff/swarmParam.redMobility;
      
      % evaluate electron temperature
      swarmParam.Te = (2./3.)*swarmParam.meanEnergy;
      boltzmann.workCond.update('electronTemperature', swarmParam.Te);
      
      % store swarm parameters information in the boltzmann properties
      boltzmann.swarmParam = swarmParam;
      
    end
    
    function [rateCoeffAll, rateCoeffExtra] = evaluateRateCoeff(boltzmann)
      
      % initialize rateCoeffAll and rateCoeffExtra structures
      rateCoeffAll = struct.empty;
      rateCoeffExtra = struct.empty;
      
      % evaluate rate coefficient for all collision 
      for gas = boltzmann.gasArray
        % collisions taken into account for solving the eedf
        for collision = gas.collisionArray
          rateCoeffAll(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(boltzmann.eedf);
          rateCoeffAll(end).value = [ineRate, supRate];
          rateCoeffAll(end).energy = collision.threshold;
          rateCoeffAll(end).collDescription = collision.description;
        end
        % collisions not taken into account for solving the eedf
        for collision = gas.collisionArrayExtra
          rateCoeffExtra(end+1).collID = collision.ID;
          [ineRate, supRate] = collision.evaluateRateCoeff(boltzmann.eedf);
          rateCoeffExtra(end).value = [ineRate, supRate];
          rateCoeffExtra(end).energy = collision.threshold;
          rateCoeffExtra(end).collDescription = collision.description;
        end
      end
      
      % store rate coefficients information in the boltzmann properties
      boltzmann.rateCoeffAll = rateCoeffAll;
      boltzmann.rateCoeffExtra = rateCoeffExtra;
      
    end
    
    function evaluateFirstAnisotropy(boltzmann)
    % evaluateFirstAnisotropy evaluates the value of the first anisotropy of the electron distribution function in
    % energy space in the framework of the two term approximation.
      
      % local copy of variables
      localEedf = boltzmann.eedf;                             % electron energy distribution function (isotropic)
      energyCell = boltzmann.energyGrid.cell;                 % values of energy at cell position (same as eedf)
      energyStep = boltzmann.energyGrid.step;                 % energy step of the energy grid
      localTotalCrossSection = boltzmann.totalCrossSection;   % total momentum transfer cross section
      gamma = Constant.gamma;                                 % gamma parameter (sqrt(2e/me))
      EoN = boltzmann.workCond.reducedFieldSI;                % reduced electric field (SI units)
      WoN = boltzmann.workCond.reducedExcFreqSI;              % reduced angular exitation frequency (SI units)
      
      % evaluate derivative of the eedf
      eedfDerivative = zeros(size(localEedf));
      eedfDerivative(1) = (localEedf(2)-localEedf(1))/energyStep;                     % 1st order forward approximation
      eedfDerivative(end) = (localEedf(end)-localEedf(end-1))/energyStep;             % 1st order backward approximation
      eedfDerivative(2:end-1) = (localEedf(3:end)-localEedf(1:end-2))/(energyStep*2); % 2nd order centered approximation
      
      % evaluate total momentum transfer cross section at cell positions
      totalCrossSectionCell = (localTotalCrossSection(1:end-1)+localTotalCrossSection(2:end))/2.0;
      
      % evaluate the first anisotropy
      if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
        switch boltzmann.eDensGrowthModel
          case 'temporal'
            totalCrossSectionCell = totalCrossSectionCell + (boltzmann.CIEff/gamma)./(sqrt(energyCell));
            if WoN == 0
              boltzmann.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
            else
              boltzmann.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+(WoN/gamma)^2./...
                (energyCell.*totalCrossSectionCell));
            end
          case 'spatial'
            boltzmann.firstAnisotropy = -(boltzmann.alphaRedEff.*localEedf+EoN.*eedfDerivative)./totalCrossSectionCell;
        end
      elseif WoN == 0
        boltzmann.firstAnisotropy = -EoN.*eedfDerivative./totalCrossSectionCell;
      else
        boltzmann.firstAnisotropy = -EoN*sqrt(2).*eedfDerivative./(totalCrossSectionCell+(WoN/gamma)^2./...
          (energyCell.*totalCrossSectionCell));
      end
      
    end
    
  end
  
end

function eedf = matrixInversion(matrix, energyGrid)
% matrixInversion introduce the normalization condition intro de matrix in the arguments and obtain the solution by
% direct inverion (mldivide matlab routine) and returns the corresponding eedf

  % local copies of energy grid variables
  energyCell = energyGrid.cell;
  energyStep = energyGrid.step;
  cellNumber = energyGrid.cellNumber;

  % include normalization condition for the EEDF in the Boltzmann matrix
  matrix(1,:) = matrix(1,:) + energyStep*sqrt(energyCell);

  % invert the Boltzmann matrix
  eedf = (matrix\([1 zeros(1,cellNumber-1)]'))';

  % renormalize of the EEDF
  eedf = eedf/dot(eedf,sqrt(energyCell)*energyStep);

  end

function derivatives = eedfTimeDerivative(time,variables,boltzmann,clearPersistentVars, electronDensityIsTimeDependent)
% eedfTimeDerivative evaluates the time derivative of the EEDF at a given time t, assuming that the EEDF at this precise
% time is given by eedf. 
  
  % separate variables into components
  eedf = variables(1:end-1);    % electron energy distribution function
  ne = variables(end);          % electron density (SI units)
  
  % local copy of fundamental constants
  persistent e;
  persistent e0;
  persistent gamma;
  if isempty(e)
    e = Constant.electronCharge;                % electron charge (SI units)
    gamma = Constant.gamma;                     % gamma parameter (sqrt(2e/me))
    e0 = Constant.vacuumPermittivity;           % vacuum permitivity (SI units)
  end
  
  % local copy of other simulation constants and values
  persistent cellNumber;
  persistent energyStep;
  persistent energyCell;
  persistent energyNode;
  persistent N;
  persistent WoN;
  if isempty(cellNumber)
    cellNumber = boltzmann.energyGrid.cellNumber;   % number of cells used in the energy grid
    energyStep = boltzmann.energyGrid.step;         % energy step of the energy grid
    energyCell = boltzmann.energyGrid.cell;         % energy at cells of the energy grid
    energyNode = boltzmann.energyGrid.node;         % energy at nodes of the energy grid
    N = boltzmann.workCond.gasDensity;              % gas density (SI units)
    WoN = boltzmann.workCond.reducedExcFreqSI;      % reduced angular exitation frequency (SI units)
  end
  
  % evaluation of the reduced electric field (pulsed or constant value)
  if boltzmann.isTimeDependent
    EoN = boltzmann.pulseFunction(time,boltzmann.pulseFunctionParameters)*1e-21;
  else
    EoN = boltzmann.workCond.reducedFieldSI;
  end
    
  % renormalize EEDF
  eedf = eedf/sum(eedf.*sqrt(energyCell')*energyStep);
  
  % evaluate time idependent elements of the Boltzmann equation (persistent for performance reasons)
  persistent matrix;
  persistent fieldMatrix;
  persistent nonConservativeMatrix;
  persistent CIEffIntegrand;
  persistent totalCrossSectionAux;
  persistent cellTotalCrossSection;
  persistent growthMatrixDiagElements;
  persistent growthMatrixSupElements;
  persistent growthMatrixInfElements;
  persistent g_extraFieldSpatialGrowth;
  persistent eeMatrixAuxA;
  persistent eeMatrixAuxB;
  if isempty(matrix)
    % evaluate basic boltzmann matrix (without field, ionization, attachment or e-e collisions operators)
    matrix = boltzmann.elasticMatrix + boltzmann.CARMatrix + boltzmann.inelasticMatrix;
    % save local copy of the regular field operator (without growth models)
    fieldMatrix = boltzmann.fieldMatrix;
    
    % include ionization operator (either conservative or not)
    nonConservativeMatrix = zeros(cellNumber);
    if boltzmann.includeNonConservativeIonization
      matrix = matrix + boltzmann.ionizationMatrix;
      nonConservativeMatrix = nonConservativeMatrix + boltzmann.ionizationMatrix;
    else
      matrix = matrix + boltzmann.ionizationConservativeMatrix;
    end
    
    % include attachment operator (either conservative or not)
    if boltzmann.includeNonConservativeAttachment
      matrix = matrix + boltzmann.attachmentMatrix;
      nonConservativeMatrix = nonConservativeMatrix + boltzmann.attachmentMatrix;
    else
      matrix = matrix + boltzmann.attachmentConservativeMatrix;
    end
    
    % evaluate time independent elements of the growth model operators (in case they are activated)
    if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
      % evaluate integrand to calculate the effective ionization rate
      CIEffIntegrand = sum(nonConservativeMatrix)';
      % local copies of total momentum transfer cross sections (at nodes and cells)
      totalCrossSectionAux = boltzmann.totalCrossSection;
      cellTotalCrossSection = 0.5*(totalCrossSectionAux(1:end-1) + totalCrossSectionAux(2:end));
      % choose growth model for the electron density (either spatial or temporal)
      switch boltzmann.eDensGrowthModel
        case 'temporal'
          % evaluate matrix elements of the temporal growth operator (e-density time variation term, dn/dt)
          growthMatrixDiagElements = -sqrt(energyCell)/gamma;
        case 'spatial'
          % evaluate the diffusion and mobility components of the spatial growth operator
          growthMatrixDiagElements = ...    % diffusion component
            energyCell./(3*cellTotalCrossSection);
          growthMatrixSupElements = ...     % mobility component (diag sup)
            1/(6*energyStep)*[0 energyCell(1:cellNumber-1)./(cellTotalCrossSection(1:cellNumber-1))];
          growthMatrixInfElements = ...     % mobility component (diag inf)
            -1/(6*energyStep)*[energyCell(2:cellNumber)./(cellTotalCrossSection(2:cellNumber)) 0];
          % evaluate components of the extra electric field operator of the spatial growth model
          g_extraFieldSpatialGrowth = energyNode/energyStep./(6*totalCrossSectionAux);
          g_extraFieldSpatialGrowth(1) = 0;
          g_extraFieldSpatialGrowth(end) = 0;
      end
    end
    
    % evaluate time idependent elements of the e-e collision operator (in case that e-e collisions are activated)
    if boltzmann.includeEECollisions
      % evaluating auxiliary matrix used in the evaluation of the ee upflux vector without multiplicative constant
      eeMatrixAuxA = zeros(cellNumber);
      auxEnergyArray = -(energyStep/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
      % because of power conservation, terms on the last row (cellNumber,:) and on the first column (:,1) are zero
      for k=1:cellNumber-1
        eeMatrixAuxA(k,2:k) = auxEnergyArray(2:k);
        eeMatrixAuxA(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
      end
      % detailed balance condition
      for k=1:cellNumber-1
        eeMatrixAuxA(k,2:cellNumber) = sqrt(eeMatrixAuxA(k,2:cellNumber).*eeMatrixAuxA(1:cellNumber-1,k+1)');
      end
      % evaluating auxiliary matrix used in the evaluation of the ee downflux vector without multiplicative constant
      eeMatrixAuxB = transpose(eeMatrixAuxA);
    end
    
  end
  
  % evaluate time derivative (time dependent) of each discrete component of the eedf (except e-e collisions operator)
  if boltzmann.includeNonConservativeIonization || boltzmann.includeNonConservativeAttachment
    % calculation of the current effective ionization rate (needed for both growth models, spatial and temporal)
    CIEff = gamma*energyStep*dot(eedf,CIEffIntegrand);
    % copy CIEff into boltzmann properties (needed for power balance and swarm parameters calculation)
    boltzmann.CIEff = CIEff;
    % choose growth model for the electron density (either spatial or temporal)
    switch boltzmann.eDensGrowthModel
      case 'temporal'
        % writing the modified total momentum transfer cross section (total momentum transfer cross section plus the
        % ionization rate divided by the electron velocity) needed for the reevaluation of the field operator
        totalCrossSectionMod = totalCrossSectionAux + (CIEff/gamma)./sqrt(energyNode);
        % writing the electric field operator matrix of the temporal growth model (with totalCrossSectionMod)
        g_fieldTemporalGrowthAux = energyNode./( 3*energyStep^2*totalCrossSectionMod.* ...
          (1+(WoN/gamma)^2./(energyNode.*totalCrossSectionMod.*totalCrossSectionMod)) );
        g_fieldTemporalGrowthAux(1) = 0;
        g_fieldTemporalGrowthAux(end) = 0;
        % copy field terms into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.g_fieldTemporalGrowth = g_fieldTemporalGrowthAux*energyStep^2;
        % evaluate time derivative of each discrete component of the eedf (case with temporal growth models)
        dfdt = N*gamma*( ( ... % multiplicative constant to obtain proper units of time
          matrix*eedf + ...           % full basic boltzmann matrix
          (EoN^2*( -g_fieldTemporalGrowthAux(1:cellNumber)-g_fieldTemporalGrowthAux(2:cellNumber+1)) + ... % time dependent diagonal component
          CIEff*growthMatrixDiagElements)'.*eedf + ...
          EoN^2*[ 0 g_fieldTemporalGrowthAux(2:cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...  % time dependent inf. diagonal component
          EoN^2*[ g_fieldTemporalGrowthAux(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...      % time dependent sup. diagonal component
          )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
      case 'spatial'
        % evaluation of gas density times diffusion coefficient and mobility times the electric field
        ND = gamma*energyStep*dot(growthMatrixDiagElements,eedf);
        muE = -gamma*energyStep*EoN*dot(growthMatrixSupElements+growthMatrixInfElements,eedf);
        % calculation of the effective reduced first Townsend coefficient
        if muE^2-4*CIEff*ND < 0
          alphaRedEff = CIEff/muE;
        else
          alphaRedEff = (muE - sqrt(muE^2-4*CIEff*ND))/(2*ND);
        end
        % copy alphaRedEff into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.alphaRedEff = alphaRedEff;
        % copy field terms into boltzmann properties (needed for power balance and swarm parameters calculation)
        boltzmann.g_fieldSpatialGrowth = alphaRedEff*energyStep*g_extraFieldSpatialGrowth;
        % evaluate time derivative of each discrete component of the eedf (case with spatial growth models)
        dfdt = N*gamma*( ( ...   % multiplicative constant to obtain proper units of time
          matrix*eedf+...               % full basic boltzmann matrix (without field nor growth operators)
          (EoN^2*fieldMatrix(1:cellNumber+1:cellNumber*cellNumber) - ...              % time dependent diagonal component
          alphaRedEff*EoN*(g_extraFieldSpatialGrowth(1:cellNumber)-g_extraFieldSpatialGrowth(2:cellNumber+1)) + ...
          alphaRedEff^2*growthMatrixDiagElements)'.*eedf + ...
          [ 0 EoN^2*(fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) - ...          % time dependent inf. diagonal component
          alphaRedEff*EoN*g_extraFieldSpatialGrowth(2:cellNumber) + ...
          alphaRedEff*EoN*growthMatrixInfElements(1:cellNumber-1)) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
          [ EoN^2*fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) + ...  % time dependent sup. diagonal component
          alphaRedEff*EoN*g_extraFieldSpatialGrowth(2:cellNumber) + ...
          alphaRedEff*EoN*growthMatrixSupElements(2:cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
          )./sqrt(energyCell'));        % divide by the square root of the energy to obtain derivative of the EEDF
    end
  else
    % evaluate time derivative of each discrete component of the eedf (case without growth models)
    dfdt = N*gamma*( ( ... % multiplicative constant to obtain proper units of time
      matrix*eedf + ...           % full basic boltzmann matrix (without field operator contribution)
      EoN^2*(fieldMatrix(1:cellNumber+1:cellNumber*cellNumber))'.*eedf + ...
      EoN^2*[ 0 fieldMatrix(2:cellNumber+1:cellNumber*cellNumber) ]'.*[ 0; eedf(1:cellNumber-1) ] + ...
      EoN^2*[ fieldMatrix(cellNumber+1:cellNumber+1:cellNumber*cellNumber) 0 ]'.*[ eedf(2:cellNumber); 0 ] ...
      )./sqrt(energyCell'));      % divide by the square root of the energy to obtain derivative of the EEDF
  end
  
  % evaluate e-e contribution to the time derivative (time dependent) of each discrete component of the eedf
  if boltzmann.includeEECollisions
    % electron temperature in eV Te = (2/3)*meanEnergy
    Te = (2/3)*energyStep*dot(energyCell.^(3/2),eedf);
    % Coulomb logarithm
    logC = log(12*pi*(e0*Te/e)^(3/2)/sqrt(ne));
    % multiplicative constant
    eeConstant = (ne/N)*(e^2/(8*pi*e0^2))*logC;
    % calculation of electron-electron collisions vectors of upflux (A) and downflux (B)
    A = (eeConstant/energyStep)*(eeMatrixAuxA*eedf);
    B = (eeConstant/energyStep)*(eeMatrixAuxB*eedf);
    % copy A and B into boltzmann properties (needed for power balance calculations)
    boltzmann.Afinal = A;
    boltzmann.Bfinal = B;
    % add contribution to time derivative of each discrete component of the eedf due to e-e collisions
    dfdt = dfdt + N*gamma*( ( ... % multiplicative constant to obtain proper units of time
      (-A(1:cellNumber)-B(1:cellNumber)).*eedf + ...                       % time dependent diagonal component
      [ 0; A(1:cellNumber-1) ].*[ 0; eedf(1:cellNumber-1) ] + ...          % time dependent inf. diagonal component
      [ B(2:cellNumber); 0 ].*[ eedf(2:cellNumber); 0 ] ...                % time dependent sup. diagonal component
      )./sqrt(energyCell'));    % divide by the square root of the energy to obtain derivative of the EEDF
  end
  
  % flush persistent memory for a new integration of the Boltzmann equation
  if clearPersistentVars
    vars = whos;
    vars = vars([vars.persistent]);
    varName = {vars.name};
    clear(varName{:});
    derivatives = [];
    return
  end
  
  % collect derivatives of the different variables
  if electronDensityIsTimeDependent
    derivatives = [dfdt; ne*boltzmann.workCond.gasDensity*CIEff];
  else
    derivatives = [dfdt; 0];
  end
  
end

function [value,isTerminal,direction] = steadyStateEventFcn(~,eedf,boltzmann,~,~)
% evaluate time derivative of each discrete component of the eedf

persistent eedfOld;
persistent maxEedfDifference;
if isempty(eedfOld)
  eedfOld = eedf;
  maxEedfDifference = boltzmann.maxEedfRelError;
  value = 1;
elseif max(abs(eedf-eedfOld)./eedfOld) < maxEedfDifference
  value = 0;
else
  eedfOld = eedf;
  value = 1;
end

isTerminal = 1;   % stop the integration
direction = 0;    % stop in any direction

end

function status = odeProgressBar(t,~,flag,varargin)

  persistent progressFigure;
  persistent progressGraph;
  persistent integrationTimeStr;
  persistent progressBar;
  persistent initialClock;
  
  switch(flag)
    case 'firstInit'
      firstTimeStep = varargin{1};
      finalTime = varargin{2};
      screenSize = get(groot,'ScreenSize');
      progressFigure = figure('Name', 'Pulse temporal integration progress bar', 'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'Position', [floor(screenSize(3)/4) floor(screenSize(4)*5/12) floor(screenSize(3)/2) floor(screenSize(4)/6)]);
      integrationTimeStr = uicontrol('Parent', progressFigure, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.11 0.65 0.89 0.15], 'HorizontalAlignment', 'left', 'String', 'Computational time: 0 s');
      progressGraph = axes('Parent', progressFigure, 'Units', 'normalized', 'OuterPosition', [0 0 1 0.6], ...
        'Box', 'on', 'XScale', 'log', 'Ytick', [], 'Xlim', [firstTimeStep finalTime]); 
      xlabel('Integration time (s)');
      hold on;
      progressBar = area(progressGraph, [firstTimeStep firstTimeStep], [1 1]);
      initialClock = clock;
    case 'lastDone'
      close(progressFigure)
      clear progressBar integrationTimeStr progressGraph progressFigure initialClock
    case 'init'
      
    case 'done'
      
    otherwise
      progressBar.XData = [1e-18 t(end)];
      integrationTimeStr.String = sprintf('Computational time: %.1f s', etime(clock, initialClock));
  end
  status = 0;
  drawnow;
end
