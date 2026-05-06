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

classdef Setup < handle
  
  properties
    
    % ---- configuration properties ----
    
    fileName = '';                  % setup file name
    unparsedInfo = cell.empty;      % cell array with the unparsed setup information
    info = struct.empty;            % parsed setup information needed to set up the components of the simulation
    
    enableGui = false;              % determines if the gui must be activated (default value false)
    enableOutput = false;           % determines if the output must be activated (default value false)
    enableElectronKinetics = false; % determines if the electron kinetics module must be activated (default value false)
    pulsedSimulation = false;       % determines if the simulation must be run in steady-state or pulsed configuration
    pulseInfo;                      % information about the pulse to be simulated
    
    batches = struct.empty;         % information about the different jobs to be run
    numberOfBatches = 1;            % total number of batches of jobs to be run (one by default)
    jobMatrixSize = [];             % dimensions of the matrix of jobs to be run
    numberOfJobs = 1;               % total number of jobs to be run (one by default)
    currentJobID = 1;               % ID of the job currently running (one by default, initial value)
    
    % ---- objects of the simulation ----
    
    workCond;                       %
    cli;                            % -> general objectes of the simulation
    gui;                            %
    output;                         %
    
    electronKinetics;                     %
    electronKineticsGasArray;             %
    electronKineticsStateArray;           % -> objects related with the electron kinetics module
    electronKineticsCollisionArray;       %
    energyGrid;                           %
    
    % ---- time counters ----
    simulationStartTime;
    
  end

  events
    genericStatusMessage;
  end
  
  methods (Access = public)
    
    function setup = Setup(fileName)
      
      % start (and save) stopwatch for the simulation
      setup.simulationStartTime = tic;

      % create CLI object
      setup.cli = CLI(setup);
      
      % parse setup file
      setup.fileName = fileName;
      notify(setup, 'genericStatusMessage', StatusEventData('Parsing setup file ...\n', 'status'));
      start = tic;

      if endsWith(fileName, '.json')

        [setupStruct, setupCell] = Parse.setupFileJson(fileName);
        str = sprintf('  Finished (%f seconds).\\n', toc(start));
        notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
        % save unstructured setup info for output and GUI
        setup.unparsedInfo = setupCell;
        % save structured info to set up the different components of the simulation
        setup.info = setupStruct;

      else

        [setupStruct, setupCell] = Parse.setupFile(fileName);
        str = sprintf('  Finished (%f seconds).\\n', toc(start));
        notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
        % save unstructured setup info for output and GUI
        setup.unparsedInfo = setupCell;
        % save structured info to set up the different components of the simulation
        setup.info = setupStruct;
      
      end

      % Perform a diagnostic of the correctness of the configuration provided by the user
      notify(setup, 'genericStatusMessage', StatusEventData('Performing selfdiagnostic of the setup file ...\n', 'status'));
      setup.selfDiagnostic();
      str = sprintf('  Finished (%f seconds).\\n', toc(start));
      notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      
    end
    
    function electronKinetics = initializeSimulation(setup)
    % initializeSimulation creates the objects necessary to run the simulation specified in a particular setup.
      
      % ----- INITIAL MESSAGE -----
      initializationStartTime = tic;
      notify(setup, 'genericStatusMessage', StatusEventData('Initializing simulation:\n', 'status'));
      
      % ----- SET SIMULATION PATH -----
      path(path, [pwd filesep 'PropertyFunctions']);
      dirContent = dir('PropertyFunctions');
      for idx = find([dirContent.isdir]==1)
        if ~any(strcmp(dirContent(idx).name, {'.' '..'}))
          path(path, [dirContent(idx).folder filesep dirContent(idx).name]);
        end
      end
      path(path, [pwd filesep 'OtherAuxFunctions']);
      dirContent = dir('OtherAuxFunctions');
      for idx = find([dirContent.isdir]==1)
        if ~any(strcmp(dirContent(idx).name, {'.' '..'}))
          path(path, [dirContent(idx).folder filesep dirContent(idx).name]);
        end
      end
      
      % ----- SETTING UP THE WORKING CONDITIONS OF THE SIMULATION -----
      % setting up the working conditions
      setup.workCond = WorkingConditions(setup);

      % ----- SETTING UP THE ELECTRON KINETICS (LoKI-B) -----
      if isfield(setup.info, 'electronKinetics') && setup.info.electronKinetics.isOn
        setup.enableElectronKinetics = true;
        start = tic;
        notify(setup, 'genericStatusMessage', StatusEventData('\t- Setting up electron kinetics ...\n', 'status'));
        % setting up the gas mixture that is going to be used to solve the electron kinetics
        setup.createElectronKineticsGasMixture();
        % create energy grid object and save handle for later use
        setup.energyGrid = Grid(setup.info.electronKinetics.numerics.energyGrid);
        % adjusting (and linking) the collision's cross sections to the energy grid
        Collision.adjustToEnergyGrid(setup.energyGrid, setup.electronKineticsCollisionArray);
        % setting up (and linking) the electron kinetics solver
        electronKinetics = setup.createElectronKinetics();
        str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
        notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      else
        electronKinetics = [];
      end
      
      % ----- SETTING UP JOBS INFORMATION -----
      start = tic;
      notify(setup, 'genericStatusMessage', StatusEventData('\t- Setting up other components ...\n', 'status'));
      % (partially implemented if electron kinetics is enabled and chemistry is disabled)
      if setup.enableElectronKinetics
        % select working conditions eligeble for jobs
        if strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann')
          auxWorkingConditions.reducedField = setup.info.workingConditions.reducedField;
        elseif strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
          auxWorkingConditions.electronTemperature = setup.info.workingConditions.electronTemperature;
        else
          auxWorkingConditions = struct();
        end
        % evaluate jobs
        for field = fieldnames(auxWorkingConditions)'
          jobs = length(auxWorkingConditions.(field{1}));
          if jobs > 1
            setup.batches(end+1).jobs = jobs;
            setup.batches(end).property = field{1};
            setup.batches(end).value = auxWorkingConditions.(field{1});
          end
        end
        % initialize jobs related variables
        setup.numberOfBatches = length(setup.batches);
        for batch = setup.batches
          setup.jobMatrixSize(end+1) = batch.jobs;
          setup.numberOfJobs = setup.numberOfJobs*setup.jobMatrixSize(end);
        end
      end
      str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
      notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      
      % ----- SETTING UP THE CLI -----
      start = tic;
      notify(setup, 'genericStatusMessage', StatusEventData('\t- Configuring CLI ...\n', 'status'));
      setup.cli.configure();
      str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
      notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));

      % ----- SETTING UP THE GUI AND OUTPUT OF THE SIMULATION -----
      % setting up the GUI object
      if isfield(setup.info, 'gui') && setup.info.gui.isOn
        setup.enableGui = true;
        start = tic;
        notify(setup, 'genericStatusMessage', StatusEventData('\t- Setting up GUI ...\n', 'status'));
        % create the object and save handle for later use
        setup.gui = GUI(setup);
        str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
        notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      end
      % setting up the output object
      if isfield(setup.info, 'output') && setup.info.output.isOn
        setup.enableOutput = true;
        start = tic;
        notify(setup, 'genericStatusMessage', StatusEventData('\t- Setting up Output ...\n', 'status'));
        % create the object and save handle for later use
        setup.output = Output(setup);
        str = sprintf('\\t    Finished (%f seconds).\\n', toc(start));
        notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      end
      
      % ----- NOTIFY END OF THE INITIALIZATION AND START OF THE CALCULATIONS -----
      str = sprintf('  Finished (%f seconds).\\n', toc(initializationStartTime));
      notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      notify(setup, 'genericStatusMessage', StatusEventData('Starting calculations: ...\n', 'status'));
      
    end
    
    function nextJob(setup)
      
      % increase the current job ID
      setup.currentJobID = setup.currentJobID + 1;
      
      % setup next job (if there are)
      if setup.currentJobID <= setup.numberOfJobs
        % obtain old and new job indices
        if isscalar(setup.jobMatrixSize) % single batch case (runs over different values of a single working condition)
          oldJobIndeces{1} = setup.currentJobID-1;
          newJobIndeces{1} = setup.currentJobID;
        else % multiple batch case (not allowed for now, runs over different values of multiple working conditions)
          [oldJobIndeces{1:setup.numberOfBatches}] = ind2sub(setup.jobMatrixSize, setup.currentJobID-1);
          [newJobIndeces{1:setup.numberOfBatches}] = ind2sub(setup.jobMatrixSize, setup.currentJobID);
        end
        
        % obtain properties that needs to be updated and the updated values
        propertiesToUpdate = cell.empty;
        newValues = [];
        for batchID = 1:setup.numberOfBatches
          if oldJobIndeces{batchID}~=newJobIndeces{batchID}
            propertiesToUpdate{end+1} = setup.batches(batchID).property;
            newValues(end+1) = setup.batches(batchID).value(newJobIndeces{batchID});
          end
        end
        % set properties for the next job
        setup.workCond.update(propertiesToUpdate, newValues);
        
        % obtain subFolder of the new run
        outputSubFolder = '';
        for i = setup.numberOfBatches:-1:1
          outputSubFolder = sprintf('%s%s%s_%g', outputSubFolder, filesep, setup.batches(i).property, ...
            setup.batches(i).value(newJobIndeces{i}));
        end
        % set subFolder for the output of the next job
        setup.output.subFolder = outputSubFolder;
        
      end
      
    end
    
    function finishSimulation(setup)
      
      % ----- RESTORE SEARCH PATH -----
      rmpath([pwd filesep 'PropertyFunctions']);
      dirContent = dir('PropertyFunctions');
      for idx = find([dirContent.isdir]==1)
        if ~any(strcmp(dirContent(idx).name, {'.' '..'}))
          rmpath([dirContent(idx).folder filesep dirContent(idx).name]);
        end
      end
      rmpath([pwd filesep 'OtherAuxFunctions']);
      dirContent = dir('OtherAuxFunctions');
      for idx = find([dirContent.isdir]==1)
        if ~any(strcmp(dirContent(idx).name, {'.' '..'}))
          rmpath([dirContent(idx).folder filesep dirContent(idx).name]);
        end
      end
      
      % ----- FINISH MESSAGE -----
      str = sprintf('  Finished (%f seconds).\n', toc(setup.simulationStartTime));
      notify(setup, 'genericStatusMessage', StatusEventData(str, 'status'));
      
    end

  end
  
  methods (Access = private)
    
    function LXCatEntryArray = parseLXCatFiles(~, fileList) 
      % parseLXCatFiles parses LXCat files (both .txt and .json) from a file list
      % Input: fileList - cell array of file paths (or single file)
      % Output: LXCatEntryArray - array of parsed LXCat entries
      
      % if it's not empty and is not a cell array, turn into cell array
      if ~isempty(fileList) & ~iscell(fileList)
        fileList = {fileList};
      end

      txtFiles = {};
      jsonFiles = {};
      
      % Distinguish between .txt and .json files
      for i = 1:length(fileList)
        if endsWith(fileList{i}, '.json')
          jsonFiles{end+1} = fileList{i};  
        else
          txtFiles{end+1} = fileList{i};
        end
      end

      % Parse .txt files
      if ~isempty(txtFiles)
        LXCatEntryArray = Parse.LXCatFiles(txtFiles);
      else
        LXCatEntryArray = [];
      end
      % Parse .json files
      if ~isempty(jsonFiles)
        LXCatEntryArray = [LXCatEntryArray, Parse.LXCatFilesJson(jsonFiles)];
      end
    end
    
    function [gasArray, stateArray, collisionArray] = LXCatData(setup)
      % LXCatData parses the LXCat files (regular or extra) included in the 
      % setup and create the corresponding arrays of gases, states and 
      % collisions. Those objects are created only with the information 
      % available in the LXCat files, so their properties (like, mases, 
      % populations, etc.) must be filled later on with the information 
      % available in the input file.
      
      % initialise arrays of gases, states and collisions
      gasArray = EedfGas.empty;
      stateArray = EedfState.empty;
      collisionArray = Collision.empty;
      

      % Parse the LXCat files using the helper function
      LXCatEntryArray = [];
      if ~isempty(setup.info.electronKinetics.LXCatFiles)
        % if it's not empty and is not a cell array, turn into cell array
        if ~iscell(setup.info.electronKinetics.LXCatFiles)
          setup.info.electronKinetics.LXCatFiles = {setup.info.electronKinetics.LXCatFiles};
        end
        LXCatEntryArray = setup.parseLXCatFiles(setup.info.electronKinetics.LXCatFiles);
      end
      
      % create gases, states and collisions from the LXCat parsed info
      for LXCatEntry = LXCatEntryArray
        [gasArray, gasID] = Gas.add(LXCatEntry.target.gasName, gasArray);
        [stateArray, targetID] = State.add(gasArray(gasID), LXCatEntry.target.ionCharg, LXCatEntry.target.eleLevel, ...
          LXCatEntry.target.vibLevel, LXCatEntry.target.rotLevel, stateArray);
        target = stateArray(targetID);
        productArray = EedfState.empty;
        productStoiCoeff = [];
        for product = LXCatEntry.productArray
          [gasArray, gasID] = Gas.add(product.gasName, gasArray);
          [stateArray, productID] = State.add(gasArray(gasID), product.ionCharg, product.eleLevel, product.vibLevel, ...
            product.rotLevel, stateArray);
          productArray(end+1) = stateArray(productID);
          productStoiCoeff(end+1) = product.quantity;
        end
        [collisionArray, ~] = Collision.add(LXCatEntry.type, LXCatEntry.isMomentumTransfer, target, productArray, ...
          productStoiCoeff, LXCatEntry.isReverse, LXCatEntry.threshold, LXCatEntry.rawCrossSection, collisionArray, ...
          false);
      end
      
      % create "extra" collisions (and correspponding gases/states) in case they are specified in the setup file
      if isfield(setup.info.electronKinetics, 'LXCatFilesExtra')
        % if it's not empty and is not a cell array, turn into cell array
        if ~isempty(setup.info.electronKinetics.LXCatFilesExtra) & ...
            ~iscell(setup.info.electronKinetics.LXCatFilesExtra)
          setup.info.electronKinetics.LXCatFilesExtra = {setup.info.electronKinetics.LXCatFilesExtra};
        end
        % Parse LXCat extra files using the helper function
        LXCatEntryArrayExtra = setup.parseLXCatFiles(setup.info.electronKinetics.LXCatFilesExtra);

                
        % create gases, states and collisions from the LXCat parsed info
        for LXCatEntry = LXCatEntryArrayExtra
          [gasArray, gasID] = Gas.add(LXCatEntry.target.gasName, gasArray);
          [stateArray, targetID] = State.add(gasArray(gasID), LXCatEntry.target.ionCharg, ...
            LXCatEntry.target.eleLevel, LXCatEntry.target.vibLevel, LXCatEntry.target.rotLevel, stateArray);
          target = stateArray(targetID);
          productArray = EedfState.empty;
          productStoiCoeff = [];
          for product = LXCatEntry.productArray
            [gasArray, gasID] = Gas.add(product.gasName, gasArray);
            [stateArray, productID] = State.add(gasArray(gasID), product.ionCharg, product.eleLevel, ...
              product.vibLevel, product.rotLevel, stateArray);
            productArray(end+1) = stateArray(productID);
            productStoiCoeff(end+1) = product.quantity;
          end
          [collisionArray, ~] = Collision.add(LXCatEntry.type, LXCatEntry.isMomentumTransfer, target, productArray, ...
            productStoiCoeff, LXCatEntry.isReverse, LXCatEntry.threshold, LXCatEntry.rawCrossSection, collisionArray, ...
            true);
        end
      end
      
      % create states needed to fix orphan states
      stateArray = State.fixOrphanStates(stateArray);
      
    end
    
    function gasProperties(setup, gasArray)
    % gasProperties fill the properties of the gases that are going to be used for solving either the electron kinetics 
    % or the heavy species kinetics
      
      % seletect type of gas (electronKinetics)
      switch class(gasArray)
        case 'EedfGas'
          gasType = 'electronKinetics';
      end
      
      % obtain properties from the setup file
      gasPropertiesAll = setup.info.(gasType).gasProperties;
      for property = fieldnames(gasPropertiesAll)'
        if ischar(gasPropertiesAll.(property{1}))
          entriesAll = {gasPropertiesAll.(property{1})};
        else
          entriesAll = gasPropertiesAll.(property{1});
        end
        for entry = entriesAll
          parsedEntry = Parse.gasPropertyEntry(entry{1});
          if isfield(parsedEntry, 'fileName')
            gasAndValueArray = Parse.gasPropertyFile(entry{1});
            for gasAndValue = gasAndValueArray
              gasID = Gas.find(gasAndValue.gasName, gasArray);
              if gasID == -1
                continue;
              end
              gasArray(gasID).(property{1}) = gasAndValue.value;
            end
          else
            gasID = Gas.find(parsedEntry.gasName, gasArray);
            if gasID == -1
              error(['Trying to assign property ''%s'' to non existing gas '...
                '''%s''.\nCheck input file'], property{1}, parsedEntry.gasName);
            end
            if isempty(parsedEntry.constant)
              gasArray(gasID).([property{1} 'Func']) = str2func(parsedEntry.function);
              gasArray(gasID).([property{1} 'Params']) = parsedEntry.argument;
              gasArray(gasID).(['evaluate' upper(property{1}(1)) property{1}(2:end)])(setup.workCond);
            else
              gasArray(gasID).(property{1}) = parsedEntry.constant;
            end
          end
        end
      end
      
      % check for the proper normalisation of the gas fractions
      Gas.checkFractionNorm(gasArray)
      
    end
    
    function stateProperties(setup, stateArray)
      % stateProperties fill the properties of the states that are necesary
      % for solving the Boltzmann equation. In particular their populations.
      
      % seletect type of state (electronKinetics)
      switch class(stateArray)
        case 'EedfState'
          stateType = 'electronKinetics';
      end
      
      % obtain properties from the setup file
      statePropertiesAll = setup.info.(stateType).stateProperties;
      for property = fieldnames(statePropertiesAll)'
        if ischar(statePropertiesAll.(property{1}))
          entriesAll = {statePropertiesAll.(property{1})};
        else
          entriesAll = statePropertiesAll.(property{1});
        end
        for entry = entriesAll
          parsedEntry = Parse.statePropertyEntry(entry{1});
          if isfield(parsedEntry, 'fileName')
            stateAndValueArray = Parse.statePropertyFile(parsedEntry.fileName);
            for stateAndValue = stateAndValueArray
              stateID = State.find(stateAndValue.gasName, stateAndValue.ionCharg, ...
                stateAndValue.eleLevel, stateAndValue.vibLevel, ...
                stateAndValue.rotLevel, stateArray);
              if stateID == -1
                continue;
              end
              stateArray(stateID).(property{1}) = stateAndValue.value;
            end
          else
            IDArray = State.find(parsedEntry.gasName, parsedEntry.ionCharg, parsedEntry.eleLevel, ...
              parsedEntry.vibLevel, parsedEntry.rotLevel, stateArray);
            if IDArray == -1
              error('Trying to asing property %s to a state which doesn''t exist.\n( %s ).\nCheck input file', ...
                property{1}, entry{1});
            end
            if isempty(parsedEntry.constant)
              for stateID = IDArray
                stateArray(stateID).([property{1} 'Func']) = str2func(parsedEntry.function);
                stateArray(stateID).([property{1} 'Params']) = parsedEntry.argument;
                stateArray(stateID).(['evaluate' upper(property{1}(1)) property{1}(2:end)])(setup.workCond);
              end
            else
              for stateID = IDArray
                stateArray(stateID).(property{1}) = parsedEntry.constant;
              end
            end
          end
        end
      end

      % evaluate state densities
      for state = stateArray
        state.evaluateDensity;
      end
      
    end
    
    function createElectronKineticsGasMixture(setup)
      % createElectronKineticsGasMixture creates the gas mixture that is going to be considered to solve the electron
      % kinetics. This includes: a gasArray with all the information related to the gas mixture, a stateArray with all
      % the states of the different gases and a collisionArray with all the electron collisions of all gases.
      %
      % This is done for a particular setup specified in the input file.
      
      % create arrays of gases, states and collisions with the information of the LXCat files
      [gasArray, stateArray, collisionArray] = setup.LXCatData();
      
      % fill properties of gases with the information of the input file
      setup.gasProperties(gasArray);
      
      % fill properties of states with the information of the input file
      setup.stateProperties(stateArray);
      
      % check for gases for which CAR is activated to meet the appropiate conditions
      if isfield(setup.info.electronKinetics, 'CARgases')
        if ischar(setup.info.electronKinetics.CARgases)
          setup.info.electronKinetics.CARgases = {setup.info.electronKinetics.CARgases};
        end
        for gasName = setup.info.electronKinetics.CARgases
          gasID = Gas.find(gasName, gasArray);
          gasArray(gasID).checkCARconditions;
        end
      end
      
      % check for non standard populations in effective cross sections
      if isfield(setup.info.electronKinetics, 'effectiveCrossSectionPopulations')
        entriesAll = setup.info.electronKinetics.effectiveCrossSectionPopulations;
        if ischar(entriesAll)
          entriesAll = {entriesAll};
        end
        for entry = entriesAll
          parsedEntry = Parse.statePropertyEntry(entry{1});
          if ~isfield(parsedEntry, 'fileName')
            error('Effective collision populations should be specified through an input file', entry{1});
          end
          stateAndValueArray = Parse.statePropertyFile(parsedEntry.fileName);
          for stateAndValue = stateAndValueArray
            stateID = State.find(stateAndValue.gasName, ...
              stateAndValue.ionCharg, stateAndValue.eleLevel, ...
              stateAndValue.vibLevel, stateAndValue.rotLevel, ...
              stateArray);
            if stateID == -1
              continue;
            end
            stateArray(stateID).gas.effectivePopulations(stateID) = stateAndValue.value;
          end
        end
      end
      
      % different checks for each gas in gasArray
      for gas = gasArray
        % avoid dummy gases (gases created because of extra cross sections or for the sake of a pretty output)
        if isempty(gas.collisionArray)
          continue;
        end
        % check for each (non dummy) gas to have its mass defined
        if isempty(gas.mass)
          error(['Mass of gas %s not found.\nNeeded for the evaluation of the elastic collision operator (Boltzmann).' ...
            '\nCheck input file'], gas.name);
        end
        % check for the distribution of states to be properly normalised
        gas.checkPopulationNorms();
        % check for Elastic collisions to be defined, for each target state
        collisionArray = gas.checkElasticCollisions(collisionArray);
      end

      % check for every collision to have non empty rawCrossSection property
      for collision = collisionArray
        if isempty(collision.rawCrossSection)
          error(['Raw cross section not found for collision %s.\nCheck the corresponding LXCat file(s).' ...
            '\nPlease, fix the problem and run the code again.'], collision.description);
        end
      end
      
      % save handles for later use
      setup.electronKineticsGasArray = gasArray;
      setup.electronKineticsStateArray = stateArray;
      setup.electronKineticsCollisionArray = collisionArray;
      
    end
    
    function electronKinetics = createElectronKinetics(setup)
      % createElectronKinetics is in charge of the creation of the eedf solver specified by the user in the setup file.
      % For the moment, the different options are:
      %   - boltzmann => Boltzmann equation solver under the classical two term approximation
      %   - prescribedEedf  => generalized eedf Maxwellian/Druyvesteyn
      
      % create the electron kinetics object
      switch lower(setup.info.electronKinetics.eedfType)
        case 'boltzmann'
          electronKinetics = Boltzmann(setup);
        case 'prescribedeedf'
          if ~isfield(setup.info.electronKinetics, 'shapeParameter')
            error(['When selectinc a ''prescribedEedf'' the user must include the '...
              '''electronKinetics->shapeParameter'' property in the setup file.']);
          end
          electronKinetics = PrescribedEedf(setup);
        otherwise
          error('''%s'' type of EEDF is not currently supported by LoKI-B\n', ...
            setup.info.electronKinetics.eedfType);
      end
      
      % save handle for later use
      setup.electronKinetics = electronKinetics;
      
    end
    
    function selfDiagnostic(setup)
      % selfDiagnostic is a function that performs a diagnostic of the values provided by the user in the setup file,
      % checking for the correctness of the configuration of the simulation.
      
      % local copy of the setup configuration
      setupInfo = setup.info;

      % define strings used in all error messages
      str1 = '\nError found in the configuration of the setup file.\n';
      str1w = '\nWarning about the configuration of the setup file.\n';
      str2 = '\nPlease, fix the problem and run the code again.';
      
      % check working conditions
      % check whether the mandatory fields of the working conditions module are present 
      workCondStruct = setupInfo.workingConditions;
      if ~isfield(workCondStruct, 'reducedField') && ~isfield(workCondStruct, 'electronTemperature')
        error([str1 '''reducedField'' and ''electronTemperature'' mandatory fields not found\n' ...
            'in the ''workingConditions'' section of the setup file.' str2],1);
      elseif ~isfield(workCondStruct, 'excitationFrequency')
        error([str1 '''excitationFrequency'' mandatory field not found' ...
            ' in the ''workingConditions'' section of the setup file.' str2],1);
      elseif ~isfield(workCondStruct, 'gasTemperature')
        error([str1 '''gasTemperature'' mandatory field not found' ...
            ' in the ''workingConditions'' section of the setup file.' str2],1);
      end  


      % check input for 'reducedField' field in working conditions to know if simulation is in pulsed mode
      if isfield(workCondStruct, 'reducedField') && ischar(workCondStruct.reducedField)
        % configuration for pulsed simulations: pulse@functionName, firstStep, finalTime, samplingType, samplingPoints, additionalParameters
        pulseInfoAux = regexp(workCondStruct.reducedField, ['pulse@\s*(?<function>\w+)\s*,\s*' ...
          '(?<firstStep>[\(\)^\d.eE+-]+)\s*,\s*(?<finalTime>[\(\)^\d.eE+-]+)\s*,\s*(?<samplingType>linspace|logspace)\s*,\s*' ...
          '(?<samplingPoints>\d+)\s*,?\s*(?<functionParameters>.*)?'], 'names', 'once');
        if isempty(pulseInfoAux)
          error([str1 'Wrong configuration for the field ' ...
            '''workingConditions>reducedField''.' str2],'foo');
        end
        pulseInfoAux.function = str2func(pulseInfoAux.function);
        pulseInfoAux.firstStep = str2num(pulseInfoAux.firstStep);
        pulseInfoAux.finalTime = str2num(pulseInfoAux.finalTime);
        pulseInfoAux.samplingPoints = str2double(pulseInfoAux.samplingPoints);
        if ~isnumeric(pulseInfoAux.firstStep) || isnan(pulseInfoAux.firstStep)
          error([str1 'Wrong configuration for the field ' ...
            '''workingConditions>reducedField''.\n''firstStep'' parameter should be numeric.' str2],'foo');
        elseif ~isnumeric(pulseInfoAux.finalTime) || isnan(pulseInfoAux.finalTime)
          error([str1 'Wrong configuration for the field ' ...
            '''workingConditions>reducedField''.\n''finalTime'' parameter should be numeric.' str2], 'foo');
        elseif ~isnumeric(pulseInfoAux.samplingPoints) || isnan(pulseInfoAux.samplingPoints)
          error([str1 'Wrong configuration for the field ' ...
            '''workingConditions>reducedField''.\n''samplingPoints'' parameter should be numeric.' str2],'foo');
        end
        pulseInfoAux.functionParameters = regexp(pulseInfoAux.functionParameters, ',', 'split');
        for idx = 1:length(pulseInfoAux.functionParameters)
          value = str2double(pulseInfoAux.functionParameters{idx});
          if ~isnan(value)
            pulseInfoAux.functionParameters{idx} = value;
          end
        end
        % save pulse information in the setup object
        setup.pulsedSimulation = true;
        setup.pulseInfo = pulseInfoAux;
        % set reduced field in the working conditions to zero (later updated to initial value of the pulse)
        workCondStruct.reducedField = 0;
        setupInfo.workingConditions.reducedField = 0;
        setup.info.workingConditions.reducedField = 0;
      end
      % check working conditions to be positive quantities
      for field = fieldnames(workCondStruct)'
        if any(workCondStruct.(field{1})<0)
          error([str1 'Negative value found for the field' ...
            '''workinConditions>%s''.' str2], field{1});
        end
      end
      % check for gas pressure (needed to evaluate gas density) to be defined when performing HF simulations
      if workCondStruct.excitationFrequency ~= 0
        if ~isfield(workCondStruct, 'gasPressure')
          error([str1 '''gasPressure'' field not found in the ' ...
            '''workingConditions'' section of the setup file.' str2], field{1});
        elseif workCondStruct.gasPressure == 0
          error([str1 '''gasPressure'' must be different than ' ...
            'zero for HF simulations.' str2], field{1});
        end
      end

      
      % check configuration of the electron kinetic module (in case it is present in the setup file)
      if isfield(setupInfo, 'electronKinetics')
        % check whether the isOn field is present and logical. Then in case it is true the checking continues
        if ~isfield(setupInfo.electronKinetics, 'isOn')
          error([str1 '''isOn'' field not found in the ' ...
            '''electronKinetics'' section of the setup file.' str2],1);
        elseif ~islogical(setupInfo.electronKinetics.isOn)
          error([str1 'Wrong value for the field ' ...
            '''electronKinetics>isOn''.\nValue should be logical (''true'' or ''false'').' str2],1);
        elseif setupInfo.electronKinetics.isOn
          % check whether the mandatory fields of the electron kinetic module (apart from isOn) are present when the
          % electron kinetic module is activated
          % --- 'eedfType' field
          if ~isfield(setupInfo.electronKinetics, 'eedfType')
            error([str1 '''eedfType'' field not found in the ' ...
              '''electronKinetics'' section of the setup file.' str2],1);
          elseif ~any(strcmp({'boltzmann' 'prescribedEedf'}, setupInfo.electronKinetics.eedfType))
            error(['Error found in the configuration of the setup file.\nWrong value for the field ' ...
              '''electronKinetics>eedfType''.\nValue should be either: ''boltzmann'' or ''prescribedEedf''.\n' ...
              'Please, fix the problem and run the code again.'],1);
          elseif strcmp(setupInfo.electronKinetics.eedfType, 'prescribedEedf')
            % check whether the shapeParameter field is present and the value is correct (double between 1 and 2)
            if ~isfield(setupInfo.electronKinetics, 'shapeParameter')
              error([str1 '''shapeParameter'' field not found in ' ...
                'the ''electronKinetics'' section of the setup file.' str2],1);
            elseif ~isnumeric(setupInfo.electronKinetics.shapeParameter) || ...
                length(setupInfo.electronKinetics.shapeParameter)~=1 || ...
                setupInfo.electronKinetics.shapeParameter < 1 || ...
                setupInfo.electronKinetics.shapeParameter > 2
              error([str1 'Wrong value for the field ' ...
                '''electronKinetics>shapeParameter''.\nValue should a number between 1 and 2 (1 for Maxwellian and ' ...
                '2 for Druyvesteyn).' str2],1);
            end
          end
          % --- 'ionizationOperatorType' field
          if ~isfield(setupInfo.electronKinetics, 'ionizationOperatorType')
            error([str1 '''ionizationOperatorType'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~any(strcmp({'conservative' 'oneTakesAll' 'equalSharing' 'usingSDCS'}, ...
              setupInfo.electronKinetics.ionizationOperatorType))
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>ionizationOperatorType''.\nValue should be either: ''conservative'', ' ...
              '''oneTakesAll'', ''equalSharing'' or ''usingSDCS''.' str2],1);
          end
          % --- 'growthModelType' field
          if ~isfield(setupInfo.electronKinetics, 'growthModelType')
            error([str1 '''growthModelType'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~any(strcmp({'temporal' 'spatial'}, setupInfo.electronKinetics.growthModelType))
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>growthModelType''.\nValue should be either: ''temporal'' or ''spatial''.' str2],1);
          end
          % --- 'includeEECollisions' field
          if ~isfield(setupInfo.electronKinetics, 'includeEECollisions')
            error([str1 '''includeEECollisions'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~islogical(setupInfo.electronKinetics.includeEECollisions)
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>includeEECollisions''.\nValue should be logical (''true'' or ''false'').' str2],1);
          elseif setupInfo.electronKinetics.includeEECollisions
            if ~isfield(workCondStruct, 'electronDensity')
              error([str1 '''electronDensity'' (mandatory field for ''electronKinetics.includeEECollisions=true'')\n' ...
                'not found in the ''workingConditions'' section of the setup file.' str2],1);
            elseif ~isfield(workCondStruct, 'gasPressure')
              error([str1 '''gasPressure'' (mandatory field for ''electronKinetics.includeEECollisions=true'')\n' ...
                'not found in the ''workingConditions'' section of the setup file.' str2],1);              
            end
          end
          % --- 'LXCatFiles' field
          if ~isfield(setupInfo.electronKinetics, 'LXCatFiles')
            error([str1 '''LXCatFiles'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          end
          % --- 'gasProperties' field
          if ~isfield(setupInfo.electronKinetics, 'gasProperties')
            error([str1 '''gasProperties'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.gasProperties, 'mass')
            error([str1 '''mass'' field not found in the ' ...
              '''electronKinetics>gasProperties'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.gasProperties, 'fraction')
            error([str1 '''fraction'' field not found in the ' ...
              '''electronKinetics>gasProperties'' section of the setup file.' str2],1);
          end
          % --- 'stateProperties' field
          if ~isfield(setupInfo.electronKinetics, 'stateProperties')
            error([str1 '''stateProperties'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.stateProperties, 'population')
            error([str1 '''population'' field not found in the ' ...
              '''electronKinetics>stateProperties'' section of the setup file.' str2],1);
          end
          % --- 'numerics' field
          if ~isfield(setupInfo.electronKinetics, 'numerics')
            error([str1 '''numerics'' field not ' ...
              'found in the ''electronKinetics'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics, 'energyGrid')
            error([str1 '''energyGrid'' field not found in the ' ...
              '''electronKinetics>numerics'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics.energyGrid, 'maxEnergy')
            error([str1 '''maxEnergy'' field not found in the ' ...
              '''electronKinetics>numerics>energyGrid'' section of the setup file.' str2],1);
          elseif ~isnumeric(setupInfo.electronKinetics.numerics.energyGrid.maxEnergy) || ...
              length(setupInfo.electronKinetics.numerics.energyGrid.maxEnergy)~=1 || ...
              setupInfo.electronKinetics.numerics.energyGrid.maxEnergy <= 0
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>numerics>energyGrid>maxEnergy''.\nValue should be a single positive number.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics.energyGrid, 'cellNumber')
            error([str1 '''cellNumber'' field not found in the ' ...
              '''electronKinetics>numerics>energyGrid'' section of the setup file.' str2],1);
          elseif ~isnumeric(setupInfo.electronKinetics.numerics.energyGrid.cellNumber) || ...
              length(setupInfo.electronKinetics.numerics.energyGrid.cellNumber)~=1 || ...
              setupInfo.electronKinetics.numerics.energyGrid.cellNumber <= 0 || ...
              mod(setupInfo.electronKinetics.numerics.energyGrid.cellNumber,1) ~= 0
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>numerics>energyGrid>cellNumber''.\nValue should be a single positive integer.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics, 'maxPowerBalanceRelError')
            error([str1 '''maxPowerBalanceRelError'' field not ' ...
              'found in the ''electronKinetics>numerics'' section of the setup file.' str2],1);
          elseif ~isnumeric(setupInfo.electronKinetics.numerics.maxPowerBalanceRelError) || ...
              length(setupInfo.electronKinetics.numerics.maxPowerBalanceRelError) ~= 1 || ...
              setupInfo.electronKinetics.numerics.maxPowerBalanceRelError < 1e-15
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>numerics>energyGrid>maxPowerBalanceRelError''.\nValue should be a number ' ...
              'larger than 1e-15.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics, 'nonLinearRoutines')
            error([str1 '''nonLinearRoutines'' field not ' ...
              'found in the ''electronKinetics>numerics'' section of the setup file.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics.nonLinearRoutines, 'algorithm')
            error([str1 '''algorithm'' field not ' ...
              'found in the ''electronKinetics>numerics>nonLinearRoutines'' section of the setup file.' str2],1);
          elseif ~any(strcmp({'mixingDirectSolutions' 'temporalIntegration'}, ...
              setupInfo.electronKinetics.numerics.nonLinearRoutines.algorithm))
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>numerics>nonLinearRoutines>algorithm''.\nValue should be either: ' ...
              '''mixingDirectSolutions'' or ''temporalIntegration''.' str2],1);
          elseif ~isfield(setupInfo.electronKinetics.numerics.nonLinearRoutines, 'maxEedfRelError')
            error([str1 '''maxEedfRelError'' field not ' ...
              'found in the ''electronKinetics>numerics>nonLinearRoutines'' section of the setup file.' str2],1);
          elseif ~isnumeric(setupInfo.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError) || ...
              length(setupInfo.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError) ~= 1 || ...
              setupInfo.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError < 1e-15
            error([str1 'Wrong value for the field ' ...
              '''electronKinetics>numerics>nonLinearRoutines>maxEedfRelError''.\nValue should be a number ' ...
              'larger than 1e-15.' str2],1);
          end
          % check configuration of the CAR in case it is activated
          if isfield(setupInfo.electronKinetics, 'CARgases')
            if ~isfield(setupInfo.electronKinetics.gasProperties, 'rotationalConstant')
              error([str1 '''rotationalConstant'' field not found ' ...
                'in the ''electronKinetics>gasProperties'' section of the setup file.' str2],1);
            elseif ~isfield(setupInfo.electronKinetics.gasProperties, 'electricQuadrupoleMoment')
              error([str1 '''electricQuadrupoleMoment'' field not ' ...
                'found in the ''electronKinetics>gasProperties'' section of the setup file.' str2],1);
            end
          end
          % check the configuration of pulsed simulations (when simulating field pulses)
          if setup.pulsedSimulation
            % pulsed simulations only implemented for the electron kinetics
            if isfield(setupInfo, 'chemistry') && setupInfo.chemistry.isOn
              error([str1 'Pulsed simulations are not implemented for the chemistry module.' str2],1);
            end
            % pulsed simulations not allowed for oscillating fields
            if any(workCondStruct.excitationFrequency~=0)
              error([str1 'Pulsed simulations are not allowed for oscillating fields.' str2],1);
            end
            % smart grid not inplemented for pulsed simulations
            if isfield(setupInfo.electronKinetics.numerics.energyGrid, 'smartGrid')
              error([str1 '''smartGrid'' feature not implemented for pulsed simulations.' str2],1);
            end
            % pulsed simulations must use the "temporalIntegration" algorithm
            if ~strcmp(setupInfo.electronKinetics.numerics.nonLinearRoutines.algorithm, 'temporalIntegration')
              error([str1 'When doing pulsed simulations the ' ...
                'parameter ''electronKinetics>numerics>nonLinearRoutines>algorithm'' must be set to ' ...
                '''temporalIntegration''.' str2],1);
            end
            % pulsed simulations must use either: conservative ionization or temporal growth for the electron density
            if ~strcmp(setupInfo.electronKinetics.ionizationOperatorType, 'conservative') && ...
                strcmp(setupInfo.electronKinetics.growthModelType, 'spatial')
              error([str1 'When doing pulsed simulations the ' ...
                'parameter ''electronKinetics>growthModelType'' must be set to ''temporal''.' str2],1);
            end
          end
          % check configuration of the smart grid in case it is activated
          if isfield(setupInfo.electronKinetics.numerics.energyGrid, 'smartGrid')
            if ~isfield(setupInfo.electronKinetics.numerics.energyGrid.smartGrid, 'minEedfDecay')
              error([str1 '''minEedfDecay'' field not found in the ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.' str2],1);
            elseif ~isnumeric(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay) || ...
                length(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay)~=1 || ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay <= 0 || ...
                mod(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay,1) ~= 0
              error([str1 'Wrong value for the field ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid>minEedfDecay''.\nValue should be a single ' ...
                'positive integer.' str2],1);
            elseif ~isfield(setupInfo.electronKinetics.numerics.energyGrid.smartGrid, 'maxEedfDecay')
              error([str1 '''maxEedfDecay'' field not found in the ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.' str2],1);
            elseif ~isnumeric(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay) || ...
                length(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay)~=1 || ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay <= 0 || ...
                mod(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay,1) ~= 0 || ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay >= ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay
              error([str1 'Wrong value for the field ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid>maxEedfDecay''.\nValue should be a single ' ...
                'positive integer (larger than minEedfDecay).' str2],1);
            elseif ~isfield(setupInfo.electronKinetics.numerics.energyGrid.smartGrid, 'updateFactor')
              error([str1 '''updateFactor'' field not found in the ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid'' section of the setup file.' str2],1);
            elseif ~isnumeric(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.updateFactor) || ...
                length(setupInfo.electronKinetics.numerics.energyGrid.smartGrid.updateFactor)~=1 || ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.updateFactor <= 0 || ...
                setupInfo.electronKinetics.numerics.energyGrid.smartGrid.updateFactor >= 1
              error([str1 'Wrong value for the field ' ...
                '''electronKinetics>numerics>energyGrid>smartGrid>updateFactor''.\nValue should be a number ' ...
                'larger than 0 and smaller than 1.' str2],1);
            end
          end
          % check configuration of the mixingDirectSolutions algorithm in case it is selected
          if strcmp(setupInfo.electronKinetics.numerics.nonLinearRoutines.algorithm, 'mixingDirectSolutions')
            if ~isfield(setupInfo.electronKinetics.numerics.nonLinearRoutines, 'mixingParameter')
              error([str1 '''mixingParameter'' field not ' ...
                'found in the ''electronKinetics>numerics>nonLinearRoutines'' section of the setup file.' str2],1);
            elseif ~isnumeric(setupInfo.electronKinetics.numerics.nonLinearRoutines.mixingParameter) || ...
                length(setupInfo.electronKinetics.numerics.nonLinearRoutines.mixingParameter) ~= 1 || ...
                setupInfo.electronKinetics.numerics.nonLinearRoutines.mixingParameter < 0 || ...
                setupInfo.electronKinetics.numerics.nonLinearRoutines.mixingParameter > 1
              error([str1 'Wrong value for the field ' ...
                '''electronKinetics>numerics>energyGrid>numerics>nonLinearRoutines>mixingParameter''.\nValue ' ...
                'should be a number between 0 and 1.' str2],1);
            end
          end
          % check the model for electron density growth when simulating time oscillating fields
          if strcmp(setupInfo.electronKinetics.growthModelType, 'spatial') && any(workCondStruct.excitationFrequency~=0)
            error([str1 'The spatial electron density growth model can not be selected with an HF field.' str2],1);
          end
        end
      end
      
      % check for 'empty' simulations (simulations with no module activated)
      if (~isfield(setupInfo, 'electronKinetics') || ~setupInfo.electronKinetics.isOn) 
        error([str1 'Module ''electronKinetics'' is not activated.' str2],1);
      end
      
      % check configuration of the graphical user interface (in case it is present in the setup file)
      if isfield(setupInfo, 'gui')
        % check whether the isOn field is present and logical. Then in case it is true the checking continues
        if ~isfield(setupInfo.gui, 'isOn')
          error([str1 '''isOn'' field not found in the ''gui'' section of the setup file.' str2],1);
        elseif ~islogical(setupInfo.gui.isOn)
          error([str1 'Wrong value for the field ''gui>isOn''.\nValue should be logical (''true'' or ' ...
            '''false'').' str2],1);
        elseif setupInfo.gui.isOn
          if ~isfield(setupInfo.gui, 'refreshFrequency')
            error([str1 '''refreshFrequency'' field not found in the ''gui'' section of the setup file.' str2],1);
          elseif ~isnumeric(setupInfo.gui.refreshFrequency) || length(setupInfo.gui.refreshFrequency)~=1 || ...
              setupInfo.gui.refreshFrequency <= 0 || mod(setupInfo.gui.refreshFrequency,1) ~= 0
            error([str1 'Wrong value for the field ''gui>refreshFrequency''.\nValue should be a single positive ' ...
              'integer.' str2],1);
          end
        end
      end
      
      % check configuration of the output (in case it is present in the setup file)
      if isfield(setupInfo, 'output')
        % check whether the isOn field is present and logical. Then in case it is true the checking continues
        if ~isfield(setupInfo.output, 'isOn')
          error([str1 '''isOn'' field not found in the ''output'' section of the setup file.' str2],1);
        elseif ~islogical(setupInfo.output.isOn)
          error([str1 'Wrong value for the field ''output>isOn''.\nValue should be logical (''true'' or ' ...
            '''false'').' str2],1);
        elseif setupInfo.output.isOn
          if ~isfield(setupInfo.output, 'dataFormat')
            error([str1 'dataFormat'' field not found in the ''output'' section of the setup file.' str2],1);
          elseif ~contains(setupInfo.output.dataFormat, 'txt') && ~contains(setupInfo.output.dataFormat, 'hdf5')
            error([str1 'Wrong value for the field ''dataFormat''.\nValue should be ''txt'', ''hdf5'' or both.' str2],1);
          end
          if ~isfield(setupInfo.output, 'dataSets')
            error([str1 '''dataSets'' field not found in the ''output'' section of the setup file.' str2],1);
          else
            dataSets = setupInfo.output.dataSets;
            if ischar(dataSets)
              dataSets = {dataSets};
            end
            possibleDataSets = {'none' 'inputs' 'log'};
            if isfield(setupInfo, 'electronKinetics') && setupInfo.electronKinetics.isOn
              possibleDataSets = [possibleDataSets 'eedf' 'swarmParameters' 'rateCoefficients' 'powerBalance' ...
                'lookUpTables'];
            end
            possibleDataSetsStr = possibleDataSets{1};
            for idx = 2:length(possibleDataSets)-1
              possibleDataSetsStr = [possibleDataSetsStr ', ' possibleDataSets{idx}];
            end
            possibleDataSetsStr = [possibleDataSetsStr ' or ' possibleDataSets{end}];
            for dataSet = dataSets
              if ~any(strcmp(possibleDataSets, dataSet))
                error([str1 'Wrong value for the field ''output>dataSets''.\nPossible data sets are: %s.' str2], ...
                  possibleDataSetsStr);
              end
            end
          end
        end
      end
      
    end
    
  end
  
end

function [gasNames, gasID] = addGas(gasName, gasNames)

  if isempty(gasNames)
    gasNames{1} = gasName;
    gasID = 1;
  else
    for i = 1:length(gasNames)
      if strcmp(gasNames{i}, gasName)
        gasID = i;
        break;
      elseif i == length(gasNames)
        gasNames{end+1} = gasName;
        gasID = length(gasNames);
      end
    end
  end

end

function [stateGasID, stateIonLevel, stateEleLevel, stateVibLevel, stateRotLevel, stateID] = ...
  addState(gasID, ionCharg, eleLevel, vibLevel, rotLevel, ...
  stateGasID, stateIonLevel, stateEleLevel, stateVibLevel, stateRotLevel)

if isempty(stateGasID)
  stateGasID = gasID;
  stateIonLevel = {ionCharg};
  stateEleLevel = {eleLevel};
  stateVibLevel = {vibLevel};
  stateRotLevel = {rotLevel};
  stateID = 1;
else
  for i = 1:length(stateGasID)
    if (stateGasID(i) == gasID && ...
        strcmp(stateIonLevel{i}, ionCharg) && ...
        strcmp(stateEleLevel{i}, eleLevel) && ...
        strcmp(stateVibLevel{i}, vibLevel) && ...
        strcmp(stateRotLevel{i}, rotLevel))
      stateID = i;
      break;
    elseif i == length(stateGasID)
      stateGasID(end+1) = gasID;
      stateIonLevel{end+1} = ionCharg;
      stateEleLevel{end+1} = eleLevel;
      stateVibLevel{end+1} = vibLevel;
      stateRotLevel{end+1} = rotLevel;
      stateID = i+1;
      break;
    end
  end
end

end
