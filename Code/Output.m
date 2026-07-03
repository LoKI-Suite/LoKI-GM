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

classdef Output < handle
  
  properties
    folder = '';                        % main output folder
    subFolder = '';                     % sub folder (deeper level) for output of different jobs
    subFolderBatches = '';              % sub folder (higher level) for output of different jobs    
    h5file = '';                        % output file name if dataFormat is hdf5
    dataFormat = '';                    % data format to save results. Options are: 'txt' and 'hdf5'

    isSimulationHF = [];                % boolean vector (dim = number of electricField jobs) to know if the electron kinetics (Boltzmann only) is HF
    isExcitationFrequencyBatch = false; % boolean value. True if ExcitationFrequency has several values
    isBoltzmann = true;                 % boolean to know if the electron kinetics is Boltzmann (true) of prescribedEedf (false)
    isPulse = false;

    logIsToBeSaved = false;             % boolean to know if the log (as written by the CLI) must be saved
    inputsAreToBeSaved = false;         % boolean to know if the input files must be saved
    currentJobID = 1;                   % cummulative index of the job value used for data writing
    currentJobIndeces = [1];            % indeces of the job value used for data writing
    numberOfEoverNJobs = 1;             % number of distinct E/N values (or Te values in case of prescribedEedf)
    extraDims = [];                     % Number of tasks for each property, in addition to E/N, with multiple values.
    eedfIsToBeSaved = false;            % boolean to know if the eedf must be saved
    swarmParamsIsToBeSaved = false;     % boolean to know if the swarm parameters info must be saved
    rateCoeffsIsToBeSaved = false;      % boolean to know if the rate coefficients info must be saved
    powerBalanceIsToBeSaved = false;    % boolean to know if the power balance info must be saved
    lookUpTablesAreToBeSaved = false;   % boolean to know if look-up tables with results must be saved
    finalDensitiesIsToBeSaved = false;  % boolean to know if the final densities of the chemistry species must be saved
    finalTemperaturesIsToBeSaved = false;     % boolean to know if the gas-related temperatures must be saved
    finalParticleBalanceIsToBeSaved = false;  % boolean to know if the final particle balances must be saved
    finalThermalBalanceIsToBeSaved = false;   % boolean to know if the final thermal balance must be saved
    chemParamsIsToBeSaved = false;            % boolean to know if the chemistry parameters info must be saved
    chemSolutionTimeIsToBeSaved = false;      % boolean to know if the solution vs time of the chemistry must be saved
  end
  
  methods (Access = public)
    
    function output = Output(setup)
      
      setupInfo = setup.info;
      output.dataFormat = setupInfo.output.dataFormat;
      output.isPulse = setup.pulsedSimulation;

      % check if electronKinetics calculations are boltzmann of prescribedEedf
      if setup.enableElectronKinetics && strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
          output.isBoltzmann = false;
      end
      
      if contains(output.dataFormat, 'txt')
        % set output folder (if not specified in the setup, a generic folder with a timestamp is created)
        if isfield(setup.info.output, 'folder')
          output.folder = ['Output' filesep setup.info.output.folder];
        else
          output.folder = ['Output' filesep 'Simulation ' datestr(datetime, 'dd mmm yyyy HHMMSS')];
        end
        % create output folder in case it doesn't exist
        if ~isfolder(output.folder)
          mkdir(output.folder);
        end

        % set initial output subfolder (in case multiple jobs are to be run)
        outputSubFolder = '';
        iBatches = setup.numberOfBatchTypes;
        if setup.numberOfJobs > 1
          for i = setup.numberOfBatchTypes:-1:1
            outputSubFolder = sprintf('%s%s%s_%g', outputSubFolder, filesep, setup.batches(i).property, ...
            setup.batches(i).value(1));
          end
          % locate the output subfolder at next level, 
          % in case multiple jobs refer to a parameter different from
          % 'reduced field' or 'electron temperature'
          outputSubFolderBatches = '';
          if ~strcmp(setup.batches(iBatches).property, 'reducedField') && ...
              ~strcmp(setup.batches(iBatches).property, 'electronTemperature')
           % set higher-order folder for a single job or a single 'reduced field' / 'electron temperature' value
            if setup.numberOfBatchTypes == 1 || ...
              (output.isBoltzmann && isscalar(setup.info.workingConditions.reducedField)) || ...
              (~output.isBoltzmann && isscalar(setup.info.workingConditions.electronTemperature))
                firstFolder = 1;
            % set higher-order folder in other cases                
            else
                firstFolder = 2;
            end
            for i = setup.numberOfBatchTypes:-1:firstFolder
                outputSubFolderBatches = sprintf('%s%s%s_%g', outputSubFolderBatches, filesep, ...
                    setup.batches(i).property, setup.batches(i).value(1));
            end   
          end           
        end
        % save output subfolder info (folder inside the output.folder folder)
        output.subFolder = outputSubFolder;
        if iBatches ~= 0
            if ~strcmp(setup.batches(iBatches).property, 'reducedField') && ...
                ~strcmp(setup.batches(iBatches).property, 'electronTemperature')
                output.subFolderBatches = outputSubFolderBatches;
            end
        end
      end
      
      if contains(output.dataFormat, 'hdf5')
        %{
        In this case we only need the root 'Output' folder where the hdfFile
        is written. This file contains the data organized in groups which are
        the equivalent to subfoldes. The hdfFile filename is saved in output.folder.
        The hdfFile identifier is saved in output h5fid.
        %}

        % create root output folder in case it doesn't exist
        % If we also have 'txt' output format we already have an output.folder
        if ~contains(output.dataFormat, 'txt')
          if isfield(setup.info.output, 'folder')
            output.folder = ['Output' filesep setup.info.output.folder];
          else
            output.folder = ['Output' filesep 'Simulation ' datestr(datetime, 'dd mmm yyyy HHMMSS')];
          end
          if ~isfolder(output.folder)
            mkdir(output.folder);
          end
        end

        % The hdfFile name is the setup.info.output.folder with the h5 extension
        if isfield(setup.info.output, 'folder')
          hdfFile = [output.folder filesep setup.info.output.folder '.h5'];
        else
          hdfFile = [output.folder filesep 'Simulation ' datestr(datetime, 'dd mmm yyyy HHMMSS') '.h5'];
        end
        output.h5file = hdfFile;

        % Choose the number of E/N values
        for i = setup.numberOfBatchTypes:-1:1
          if (strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann') && ...
            strcmp(setup.batches(i).property, 'reducedField')) || ...
            (strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf') && ...
            strcmp(setup.batches(i).property, 'electronTemperature'))
%           if strcmp(setup.batches(i).property, 'reducedField')
            if output.isPulse
              output.numberOfEoverNJobs = setup.pulseInfo.samplingPoints+1;
            else
              output.numberOfEoverNJobs = setup.batches(i).jobs;
            end
            break;
          end
        end

        % Common constants for hdf5 calls
        dcpl = "H5P_DEFAULT";
        doubleType = H5T.copy("H5T_NATIVE_DOUBLE");
        intType = H5T.copy("H5T_NATIVE_INT");
        % creates the file with default library properties (overwrite, ...)
        fID = H5F.create(hdfFile, "H5F_ACC_TRUNC", dcpl, dcpl);

        % Saves working conditions are single valued, saved them as attributes
        workingConditions = setup.info.workingConditions;
        spaceID = H5S.create("H5S_SCALAR");
        acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
        % - Temperatures
        found = false;
        for i = 1:setup.numberOfBatchTypes
          if strcmpi(setup.batches(i).property, 'gasTemperature')
            found = true;
            if setup.batches(i).jobs == 1
              attrID = H5A.create(fID,"Gas temperature (K)",doubleType,spaceID,acpl);
              H5A.write(attrID,"H5ML_DEFAULT",workingConditions.gasTemperature);
            end
            break
          end
        end
        if ~found   % gasTemperature not in setup.batches
          attrID = H5A.create(fID,"Gas temperature (K)",doubleType,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",workingConditions.gasTemperature);
        end
        attrID = H5A.create(fID,"Wall temperature (K)",doubleType,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workingConditions.wallTemperature);
        attrID = H5A.create(fID,"External temperature (K)",doubleType,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workingConditions.extTemperature);
        % - Gas pressure
        % If chemistry is on, the gas pressure is mandatory.
        if (isfield(setupInfo, 'chemistry') && setupInfo.chemistry.isOn) || ...
          isfield(workingConditions, 'gasPressure')
          found = false;
          for i = 1:setup.numberOfBatchTypes
            if strcmpi(setup.batches(i).property, 'gasPressure')
              found = true;
              if setup.batches(i).jobs == 1
                attrID = H5A.create(fID,"Gas pressure (Pa)",doubleType,spaceID,acpl);
                H5A.write(attrID,"H5ML_DEFAULT",workingConditions.gasPressure);
              end
              break
            end
          end
          if ~found   % gasPressure not in setup.batches
            attrID = H5A.create(fID,"Gas pressure (Pa)",doubleType,spaceID,acpl);
            H5A.write(attrID,"H5ML_DEFAULT",workingConditions.gasPressure);
          end
        end
        % - Electron density
        % If chemistry or setupInfo.electronKinetics.includeEECollisions are on, ...
        % the electronDensity is mandatory; in other cases writes electronDensity if present
        if (isfield(setupInfo, 'chemistry') && setupInfo.chemistry.isOn) || ...
          setupInfo.electronKinetics.includeEECollisions || ...
          isfield(workingConditions, 'electronDensity')
          found = false;
          for i = 1:setup.numberOfBatchTypes
            if strcmpi(setup.batches(i).property, 'electronDensity')
              found = true;
              if setup.batches(i).jobs == 1
                attrID = H5A.create(fID,"Electron density (m-3)",doubleType,spaceID,acpl);
                H5A.write(attrID,"H5ML_DEFAULT",workingConditions.electronDensity);
              end
              break
            end
          end
          if ~found   % electronDensity not in setup.batches
            attrID = H5A.create(fID,"Electron density (m-3)",doubleType,spaceID,acpl);
            H5A.write(attrID,"H5ML_DEFAULT",workingConditions.electronDensity);
          end
        end

        % - Excitation frequency
          found = false;
          for i = 1:setup.numberOfBatchTypes
            if strcmpi(setup.batches(i).property, 'excitationFrequency')
              found = true;
              if setup.batches(i).jobs == 1
                attrID = H5A.create(fID,"Excitation frequency (Hz)",doubleType,spaceID,acpl);
                H5A.write(attrID,"H5ML_DEFAULT",workingConditions.excitationFrequency);
              end
              break
            end
          end
        if ~found   % excitationFrequency not in setup.batches
          attrID = H5A.create(fID,"Excitation frequency (Hz)",doubleType,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",workingConditions.excitationFrequency);
        end

        % - Surface site density
        attrID = H5A.create(fID,"Surface site density (m-2)",doubleType,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workingConditions.surfaceSiteDensity);
        % - Chamber dimensions (cylindric)
        attrID = H5A.create(fID,"Chamber length (m)",doubleType,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workingConditions.chamberLength);
        attrID = H5A.create(fID,"Chamber radius (m)",doubleType,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workingConditions.chamberRadius);
        % - Discharge current
        if isfield(workingConditions,'dischargeCurrent') 
          attrID = H5A.create(fID,"Discharge current (A)",doubleType,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",workingConditions.dischargeCurrent);
        end  
        % - Discharge power density
        if isfield(workingConditions,'dischargePowerDensity') 
          attrID = H5A.create(fID,"Discharge power density (W m-3)",doubleType,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",workingConditions.dischargePowerDensity);
        end
        % - Input gas flow
        % Only if chemistry is on, the totalSccmInFlow is written (is mandatory in this case).
        if isfield(setupInfo, 'chemistry') && setupInfo.chemistry.isOn
          if isfield(workingConditions,'totalSccmInFlow')
            attrID = H5A.create(fID,"Total input flow (sccm)",doubleType,spaceID,acpl);
            H5A.write(attrID,"H5ML_DEFAULT",workingConditions.totalSccmInFlow);
            if (~isempty(setup.workCond.totalSccmInFlow) || setup.workCond.totalSccmInFlow ~= 0) && ...
                    setup.enableChemistry
              % Gas composition
              % Write inFlowFraction for each gas
              inFlowFraction = setup.info.chemistry.gasProperties.inFlowFraction;
              for gasFlow = inFlowFraction
                name_flow = strsplit(gasFlow{1}, ' = ');
                attrID = H5A.create(fID,name_flow{1}+" fraction in input flow",doubleType,spaceID,acpl);
                H5A.write(attrID,"H5ML_DEFAULT",str2double(name_flow{2}));
              end
            end
          end
          if isfield(workingConditions,'totalSccmOutFlow')
            if isnumeric(workingConditions.totalSccmOutFlow)
              attrID = H5A.create(fID,"Total output flow (sccm)",doubleType,spaceID,acpl);
              H5A.write(attrID,"H5ML_DEFAULT",workingConditions.totalSccmOutFlow);
            elseif strcmp(workingConditions.totalSccmOutFlow, 'totalSccmInFlow')
              attrID = H5A.create(fID,"Total output flow (sccm)",doubleType,spaceID,acpl);
              H5A.write(attrID,"H5ML_DEFAULT",workingConditions.totalSccmInFlow);
              % Check what to write (if anything) in this case...
%             elseif strcmp(workCondStruct.totalSccmOutFlow, 'ensureIsobaric')
            end
          end
        end
        H5A.close(attrID);
        H5S.close(spaceID);

        % create groups for electron kinetics and chemistry
        if setup.enableElectronKinetics
          geid = H5G.create(fID,"electronKinetics",dcpl,dcpl,dcpl);
        end
        if setup.enableChemistry
          gcid = H5G.create(fID,"chemistry",dcpl,dcpl,dcpl);
        end
      end

      % save what information must be saved
      dataSets = setup.info.output.dataSets;
      if ischar(dataSets)
        dataSets = {dataSets};
      end

      % save the information if the electron kinetics is HF
      if output.isPulse
        numberOfJobs = setup.numberOfJobs*(setup.pulseInfo.samplingPoints+1);
      else
        numberOfJobs = setup.numberOfJobs;
      end
      workCond = setup.info.workingConditions;
      numberOfJobsFreq = length(workCond.excitationFrequency);
      numberOfJobsModFreq = numberOfJobs/numberOfJobsFreq;
      % DEBUG
      for idx = 1:numberOfJobsFreq
        if workCond.excitationFrequency(idx) > 0
            for idxJobs = 1:numberOfJobsModFreq
                output.isSimulationHF(end+1) = true;
            end 
        else
            for idxJobs = 1:numberOfJobsModFreq  
                output.isSimulationHF(end+1) = false;
            end
        end
      end

      % If we study the electron kinetics, we write the reducedField values,
      % electronTemperature values and all the properties with multiple values,
      % corresponding to multiple jobs, on geid.
      % If we have a time-dependent reducedField pulse, the E/N values are
      % written in the end.
      % If we study the chemistry, the maintenance field is written in the end.
      if contains(output.dataFormat, 'hdf5')
        if setup.enableElectronKinetics && strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann')
          % First saves the reducedField values if we don't do the chemistry.
          reducedField = setup.info.workingConditions.reducedField;
          % set dims dimensions
          if setup.pulsedSimulation   % reducedField(t) -> columns for t and E/N
            dims = [2 output.numberOfEoverNJobs];
          else
            dims = [1 output.numberOfEoverNJobs];
          end
          spaceID = H5S.create_simple(2,fliplr(dims),[]);
          dsID = H5D.create(geid,'reducedField',doubleType,spaceID,dcpl);
          if setup.pulsedSimulation
            H5DS.set_label(dsID,0,'time')
            H5DS.set_label(dsID,1,'E/N')
          else
            H5DS.set_label(dsID,0,'E/N')
          end
          if ~setup.enableChemistry && ~setup.pulsedSimulation
            H5D.write(dsID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', ...
              'H5P_DEFAULT',reducedField);
          end
          H5S.close(spaceID);
          % units attribute
          units = "Td";
          stypeID = H5T.copy("H5T_C_S1");
          H5T.set_size(stypeID, "H5T_VARIABLE");
          spaceID = H5S.create_simple(1,1,[]);
          acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
          attrID = H5A.create(dsID,'units',stypeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",units);
          % set as scale
          H5DS.set_scale(dsID,"E/N");
          H5T.close(stypeID);
        elseif setup.enableElectronKinetics && strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
          output.isBoltzmann = false;
          % First saves the electronTemperature values
          % get dims
          electronTemperature = setup.info.workingConditions.electronTemperature;
          if isscalar(electronTemperature)
            dims = [1 1];
          else
            dims = size(electronTemperature);
          end
          spaceID = H5S.create_simple(2,fliplr(dims),[]);
          dsID = H5D.create(geid,'electronTemperature',doubleType,spaceID,dcpl);
          H5D.write(dsID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', ...
            'H5P_DEFAULT',electronTemperature);
          H5S.close(spaceID);
          % units attribute
          units = "eV";
          stypeID = H5T.copy("H5T_C_S1");
          H5T.set_size (stypeID, "H5T_VARIABLE");
          spaceID = H5S.create_simple(1,1,[]);
          acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
          attrID = H5A.create(dsID,'units',stypeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",units);
          % set as scale
          H5DS.set_scale(dsID,"Te");
        end
        for i = 1:setup.numberOfBatchTypes
          if (strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann') && ...
            strcmp(setup.batches(i).property, 'reducedField')) || ...
            (strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf') && ...
            strcmp(setup.batches(i).property, 'electronTemperature'))
            continue    % We already wrote this dataset
          end
          dims = [1 setup.batches(i).jobs];
          spaceID = H5S.create_simple(2,fliplr(dims),[]);
          dsbID = H5D.create(geid,setup.batches(i).property,doubleType,spaceID,dcpl);
          H5DS.set_label(dsbID,0,setup.batches(i).property);
          H5D.write(dsbID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', ...
              'H5P_DEFAULT',setup.batches(i).value);
          H5S.close(spaceID);
          % units attribute
          units = setup.batches(i).units;
          stypeID = H5T.copy("H5T_C_S1");
          H5T.set_size(stypeID, "H5T_VARIABLE");
          spaceID = H5S.create_simple(1,1,[]);
          acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
          attrID = H5A.create(dsID,units,stypeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",units);
%           % set as scale
%           % DEBUG: If we want to set the property as scale, we need a name...
%           H5DS.set_scale(dsID,setup.batches(i).property);
%           H5T.close(stypeID);
        end
      end

      % Sets the common extraDims (except for reducedField) for the datasets
      for i = 1:length(setup.batches)
        if setup.enableElectronKinetics && ...
          ((strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann') && ...
          strcmp(setup.batches(i).property, 'reducedField')) || ...
          (strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf') && ...
          strcmp(setup.batches(i).property, 'electronTemperature')))
          continue
        else
          output.extraDims(end+1) = setup.batches(i).jobs;
        end
      end
      output.currentJobIndeces = ones(1,length(setup.batches));

      % Writes the datasets
      for dataSet = dataSets
        switch dataSet{1}
          case 'log'
            % Log file is always written as a txt file
            output.logIsToBeSaved = true;
            output.initializeLogFile(setup.cli.logStr);
          case 'inputs'
            % Input file is always written as a txt file
            output.inputsAreToBeSaved = true;
            output.saveInputFiles(setup);
          case 'eedf'
            output.eedfIsToBeSaved = true;
            if contains(output.dataFormat, 'hdf5')
              if strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann')
                sz(1:3) = H5T.get_size(doubleType);
                offset(1) = 0;
                offset(2:3) = cumsum(sz(1:2));
                name = ["Energy" "EEDF" "Anisotropy"];
              elseif strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
                sz(1:2) = H5T.get_size(doubleType);
                offset(1) = 0;
                offset(2) = sz(1);
                name = ["Energy" "EEDF"];
              end
              ctypeID = H5T.create ('H5T_COMPOUND', sum(sz));
              for i = 1:length(sz)
                H5T.insert(ctypeID,name(i),offset(i),doubleType);
              end
%               dims = [length(setup.energyGrid.cell) 1 output.numberOfEoverNJobs output.extraDims-1];
              dims = [length(setup.energyGrid.cell) 1 output.numberOfEoverNJobs output.extraDims];
              h5_dims = fliplr(dims);
              spaceID = H5S.create_simple(length(dims),h5_dims,h5_dims);
              dsfID = H5D.create(geid,'eedf',ctypeID,spaceID,dcpl);
              H5DS.attach_scale(dsfID,dsID,0);
              H5S.close(spaceID);

              % now create the attributes: variable and units in each column
              if strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann')
                units = ['eV       '; 'eV^-(3/2)'; 'eV^-(3/2)'];
                atdims = 3;
              elseif strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
                units = ['eV       '; 'eV^-(3/2)'];
                atdims = 2;
              end
              filetype = H5T.copy('H5T_FORTRAN_S1');
              H5T.set_size(filetype, 9);
              memtype = H5T.copy('H5T_C_S1');
              H5T.set_size(memtype, 9);
              space = H5S.create_simple(1,fliplr(atdims), []);
              attr = H5A.create(dsfID, 'Units', filetype, space, 'H5P_DEFAULT');
              H5A.write(attr, memtype, units');

              % finally close the workspaces
              H5A.close(attr);
              H5S.close(space);
              H5T.close(filetype);
              H5T.close(memtype);
              H5D.close(dsfID);
              clear sz;
            end
          case 'swarmParameters'
            output.swarmParamsIsToBeSaved = true;
            if contains(output.dataFormat, 'hdf5')
              sz(1:9) = H5T.get_size(doubleType);
              offset(1)=0;
              % get dims
             % if output.isSimulationHF(output.currentJobID)
              output.isExcitationFrequencyBatch = contains([setup.batches(:).property], 'excitationFrequency');
              if output.isSimulationHF(output.currentJobIndeces(1)) || ...
                output.isExcitationFrequencyBatch
                offset(2:9)=cumsum(sz(1:8));
                if output.isBoltzmann
                  name = ["meanEnergy" "characEnergy" "Te" "redMobility" ...
                    "redMobilityHFr" "redMobilityHFi" "redDiffCoeff" ...
                    "redMobilityEnergy" "redDiffCoeffEnergy"];
                  units = ['eV      '; 'eV      '; 'eV      '; '1/(msV) '; '1/(msV) '; ...
                    '1/(msV) '; '1/(ms)  '; 'eV/(msV)'; 'eV/(ms) '];
                else
                  name = ["meanEnergy" "characEnergy" "reducedField" "redMobility" ...
                  "redMobilityHFr" "redMobilityHFi" "redDiffCoeff" ...
                  "redMobilityEnergy" "redDiffCoeffEnergy"];
                  units = ['eV      '; 'eV      '; 'Td      '; '1/(msV) '; '1/(msV) '; ...
                    '1/(msV) '; '1/(ms)  '; 'eV/(msV)'; 'eV/(ms) '];
                end 
              else
                sz(10) = H5T.get_size(doubleType);
                offset(2:10)=cumsum(sz(1:9));
                if output.isBoltzmann
                  name = ["meanEnergy" "characEnergy" "Te" "driftVelocity" ...
                    "redMobility" "redDiffCoeff" "redMobilityEnergy" ...
                    "redDiffCoeffEnergy" "redTownsendCoeff" "redAttCoeff"];
                  units = ['eV      '; 'eV      '; 'eV      '; 'm/s     '; '1/(msV) '; ...
                    '1/(ms)  '; 'eV/(msV)'; 'eV/(ms) '; 'm2      '; 'm2      '];
                else
                  name = ["meanEnergy" "characEnergy" "reducedField" "driftVelocity" ...
                    "redMobility" "redDiffCoeff" "redMobilityEnergy" ...
                    "redDiffCoeffEnergy" "redTownsendCoeff" "redAttCoeff"];
                  units = ['eV      '; 'eV      '; 'Td      '; 'm/s     '; '1/(msV) '; ...
                    '1/(ms)  '; 'eV/(msV)'; 'eV/(ms) '; 'm2      '; 'm2      '];
                end                   
              end
              atdims = length(name);
              ctypeID = H5T.create ('H5T_COMPOUND', sum(sz));
              for i = 1:length(sz)
                H5T.insert(ctypeID,name(i),offset(i),doubleType);
              end
              dims = [output.numberOfEoverNJobs 1 output.extraDims];
              spaceID = H5S.create_simple(length(dims),fliplr(dims),[]);
              dssID = H5D.create(geid,'swarmParameters',ctypeID,spaceID,dcpl);
              H5DS.attach_scale(dssID,dsID,1);
              % attributes: variable and units in each column
              filetype = H5T.copy('H5T_FORTRAN_S1');
              H5T.set_size(filetype, 8);
              memtype = H5T.copy('H5T_C_S1');
              H5T.set_size(memtype, 8);
              space = H5S.create_simple(1,fliplr(atdims), []);
              attr = H5A.create(dssID, 'Units', filetype, space, 'H5P_DEFAULT');
              H5A.write(attr, memtype, units');
              %
              H5A.close(attr);
              H5S.close(space);
              H5T.close(filetype);
              H5T.close(memtype);
              H5T.close(ctypeID);
              H5S.close(spaceID);
              H5D.close(dssID);
              clear sz;
            end
          case 'powerBalance'
            output.powerBalanceIsToBeSaved = true;
            if contains(output.dataFormat, 'hdf5')
              % we set two datasets: powerBalanceSummary and powerBalanceGases
              % powerBalanceGases
              sz(1:5) = H5T.get_size(doubleType);
              offset(1) = 0;
              offset(2:5) = cumsum(sz(1:4));
              name = ["rotCol" "vibCol" "eleCol" "ionCol" "attCol"];
              ctypeID = H5T.create('H5T_COMPOUND', sum(sz));
              for i = 1:length(sz)
                H5T.insert(ctypeID,name(i),offset(i),doubleType);
              end
              ngas = length(setup.electronKineticsGasArray);
              dims = [output.numberOfEoverNJobs 1 3 output.extraDims ngas];
              spaceID = H5S.create_simple(length(dims),fliplr(dims),[]);
              dspID = H5D.create(geid,'powerBalanceGases',ctypeID,spaceID,dcpl);
              H5DS.attach_scale(dspID,dsID,3);
              % attributes: variable and units in each column
              units = ['eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'];
              atdims = 5;
              filetype = H5T.copy('H5T_FORTRAN_S1');
              H5T.set_size(filetype, 7);
              memtype = H5T.copy('H5T_C_S1');
              H5T.set_size(memtype, 7);
              space = H5S.create_simple(1,fliplr(atdims), []);
              attr = H5A.create(dspID, 'Units', filetype, space, 'H5P_DEFAULT');
              H5A.write(attr, memtype, units');
              %
              H5A.close(attr);
              H5S.close(space);
              H5T.close(filetype);
              H5T.close(memtype);
              H5T.close(ctypeID);
              H5S.close(spaceID);
              H5D.close(dspID);
              clear sz;
              % powerBalanceSummary
              sz(1:10) = H5T.get_size(doubleType);
              offset(1) = 0;
              offset(2:10) = cumsum(sz(1:9));
              name = ["Field" "Elastic" "CAR" "Rotational" "Vibrational" ...
                "Electronic" "Ionization" "Attachment" "eDensGrowth" "Balance"];
              ctypeID = H5T.create('H5T_COMPOUND', sum(sz));
              for i = 1:length(sz)
                H5T.insert(ctypeID,name(i),offset(i),doubleType);
              end
              dims = [output.numberOfEoverNJobs 1 3 output.extraDims];
              spaceID = H5S.create_simple(length(dims),fliplr(dims),[]);
              dspID = H5D.create(geid,'powerBalanceSummary',ctypeID,spaceID,dcpl);
              H5DS.attach_scale(dspID,dsID,2);
              % attributes: variable and units in each column
              units = ['eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; ...
                'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'; 'eVm^3/s'];
              atdims = 10;
              filetype = H5T.copy('H5T_FORTRAN_S1');
              H5T.set_size(filetype, 7);
              memtype = H5T.copy('H5T_C_S1');
              H5T.set_size(memtype, 7);
              space = H5S.create_simple(1,fliplr(atdims), []);
              attr = H5A.create(dspID, 'Units', filetype, space, 'H5P_DEFAULT');
              H5A.write(attr, memtype, units');
              %
              H5A.close(attr);
              H5S.close(space);
              H5T.close(filetype);
              H5T.close(memtype);
              H5D.close(dspID);
              clear sz;
            end
          case 'rateCoefficients'
            output.rateCoeffsIsToBeSaved = true;
            if contains(output.dataFormat, 'hdf5')
              % e-collision reactions
              % builds two tables for each gas (collisions and extraCollisions)
              temp = [setup.electronKineticsCollisionArray(1:end).isExtra];
              nReactions = length(temp(~temp));
              % Create the base types
              sz(1) = H5T.get_size(intType);
              sz(2:4) = H5T.get_size(doubleType);
              strType = H5T.copy ('H5T_C_S1');
              H5T.set_size (strType, 'H5T_VARIABLE');
              sz(5) = H5T.get_size(strType);
              % Compute the offsets to each field. The first offset is always zero.
              offset(1)=0;
              offset(2:5)=cumsum(sz(1:4));
              % Create the compound datatype for the file.
              ctypeID = H5T.create ('H5T_COMPOUND', sum(sz));
              H5T.insert(ctypeID,'rate_id',offset(1),intType);
              H5T.insert(ctypeID,'ine_coeff',offset(2),doubleType);
              H5T.insert(ctypeID,'sup_coeff',offset(3),doubleType);
              H5T.insert(ctypeID,'threshold',offset(4),doubleType);
              H5T.insert(ctypeID,'description',offset(5),strType);
              % get dims
              dims = [nReactions 1 output.numberOfEoverNJobs output.extraDims];
              spaceID = H5S.create_simple(length(dims),fliplr(dims),[]);
              dsrID = H5D.create(geid,'rateCoefficients',ctypeID,spaceID,dcpl);
              H5DS.attach_scale(dsrID,dsID,1);
              % clean-up
              H5S.close(spaceID);
              % attributes: variable and units in each column
              units = ['  -  '; 'm^3/s'; 'm^3/s'; 'eV   '; '  -  '];
              atdims = 5;
              filetype = H5T.copy('H5T_FORTRAN_S1');
              H5T.set_size(filetype, 5);
              memtype = H5T.copy('H5T_C_S1');
              H5T.set_size(memtype, 5);
              space = H5S.create_simple(1,fliplr(atdims), []);
              attr = H5A.create(dsrID, 'Units', filetype, space, 'H5P_DEFAULT');
              H5A.write(attr, memtype, units');
              %
              H5A.close(attr);
              H5S.close(space);
              H5T.close(filetype);
              H5T.close(memtype);
              H5D.close(dsrID);
              clear sz;

              % extra rate coefficients dataset
              nxReactions = length(setup.electronKineticsCollisionArray) - nReactions;
              if nxReactions > 0
                % get dims
                dims = [nxReactions 1 output.numberOfEoverNJobs output.extraDims];
                spaceID = H5S.create_simple(length(dims),fliplr(dims),[]);
                dsxID = H5D.create(geid,'extraRateCoefficients',ctypeID,spaceID,dcpl);
                H5DS.attach_scale(dsxID,dsID,1);
                % clean-up
                H5S.close(spaceID);
                % attributes: variable and units in each column
                units = ['  -  '; 'm^3/s'; 'm^3/s'; 'eV   '; '  -  '];
                atdims = 5;
                filetype = H5T.copy('H5T_FORTRAN_S1');
                H5T.set_size(filetype, 5);
                memtype = H5T.copy('H5T_C_S1');
                H5T.set_size(memtype, 5);
                space = H5S.create_simple(1,fliplr(atdims), []);
                attr = H5A.create(dsxID, 'Units', filetype, space, 'H5P_DEFAULT');
                H5A.write(attr, memtype, units');
                %
                H5A.close(attr);
                H5S.close(space);
                H5T.close(filetype);
                H5T.close(memtype);
                H5D.close(dsxID);
              end
              H5T.close(ctypeID);

              if setup.enableChemistry
                % rate coefficients for the chemistry            
                nReactions = setup.chemistry.numberOfReactions;
                % Create the base types
                sz(1) = H5T.get_size(intType);
                sz(2:5) = H5T.get_size(doubleType);
                sz(6) = H5T.get_size(strType);
                % Compute the offsets to each field. The first offset is always zero.
                offset(1)=0;
                offset(2:6)=cumsum(sz(1:5));
                % Create the compound datatype for the file.
                ctypeID = H5T.create ('H5T_COMPOUND', sum(sz));
                H5T.insert(ctypeID,'rate_id',offset(1),intType);
                H5T.insert(ctypeID,'dir_coeff',offset(2),doubleType);
                H5T.insert(ctypeID,'inv_coeff',offset(3),doubleType);
                H5T.insert(ctypeID,'enthalpy',offset(4),doubleType);
                H5T.insert(ctypeID,'net_rate',offset(5),doubleType);
                H5T.insert(ctypeID,'description',offset(6),strType);
                % get dims, spaceID and writes dataset
                dims = [nReactions 1];
                spaceID = H5S.create_simple(2,fliplr(dims),[]);
                dsrID = H5D.create(gcid,'rateCoefficients',ctypeID,spaceID,dcpl);
                % clean-up
                H5S.close(spaceID);
                H5T.close(ctypeID);
                % attributes: variable and units in each column
                units = ['  -   '; 'S.I.  '; 'S.I.  '; 'eV    '; 'm^-3/s'; '  -   '];
                atdims = 6;
                filetype = H5T.copy('H5T_FORTRAN_S1');
                H5T.set_size(filetype, 6);
                memtype = H5T.copy('H5T_C_S1');
                H5T.set_size(memtype, 6);
                space = H5S.create_simple(1,fliplr(atdims), []);
                attr = H5A.create(dsrID, 'Units', filetype, space, 'H5P_DEFAULT');
                H5A.write(attr, memtype, units');
                %
                H5A.close(attr);
                H5S.close(space);
                H5T.close(filetype);
                H5T.close(memtype);
                H5D.close(dsrID);
                clear sz;
              end
              % final clean-up
              H5T.close(strType);
              H5T.close(intType);
            end
          case 'lookUpTables'
            % lookUpTables is always written as a txt file
            output.lookUpTablesAreToBeSaved = true;
          case 'finalDensities'
            % as they are only known at the end, the datasets finalDensities,
            % electronKineticsPopulations, particleBalance and chemSolutionTime
            % are only created at the end.
            output.finalDensitiesIsToBeSaved = true;
          case 'finalTemperatures'
            output.finalTemperaturesIsToBeSaved = true;
          case 'finalParticleBalance'
            output.finalParticleBalanceIsToBeSaved = true;
          case 'finalThermalBalance'
            output.finalThermalBalanceIsToBeSaved = true;
          case 'chemSolutionTime'
            output.chemSolutionTimeIsToBeSaved = true;
          case 'chemParameters'
            output.chemParamsIsToBeSaved = true;
        end     % switch dataSet{1}
      end     % for dataSet = dataSets

      % closes the hdf5 objects
      if contains(output.dataFormat, 'hdf5')
        H5D.close(dsID);    % We still had the scale dataset open...
        if setup.enableElectronKinetics
          H5G.close(geid);
        end
        if setup.enableChemistry
          H5G.close(gcid);
        end
%         H5S.close(spaceID);            % this is not realy needed as all ...
        H5T.close(doubleType);         % identifiers are closed when they ...
        H5F.close(fID);                % go out of scope if inside a function.
      end

      % save the setup information for reference (always saved)
      output.saveSetupInfo(setup.unparsedInfo);

      % add listener to status messages of the setup object
      addlistener(setup, 'genericStatusMessage', @output.genericStatusMessage);
      % add listener of the working conditions object
      addlistener(setup.workCond, 'genericStatusMessage', @output.genericStatusMessage);
      
      if setup.enableChemistry
        % add listener to status messages of the chemistry object
        addlistener(setup.chemistry, 'genericStatusMessage', @output.genericStatusMessage);
        % add listener to output log info when a new iteration of the pressure cycle is found
        addlistener(setup.chemistry, 'newPressureCycleIteration', @output.newPressureCycleIteration);
        % add listener to output log info when a new iteration of the neutrality cycle is found
        addlistener(setup.chemistry, 'newNeutralityCycleIteration', @output.newNeutralityCycleIteration);
        % add listener to output log info when a new iteration of the global cycle is found
        addlistener(setup.chemistry, 'newGlobalCycleIteration', @output.newGlobalCycleIteration);
        
        % add listener to output results when a new solution for the Chemistry is found
        addlistener(setup.chemistry, 'obtainedNewChemistrySolution', @output.chemistrySolution);
        if setup.enableElectronKinetics
          % add listener to status messages of the electron kinetics object
          addlistener(setup.electronKinetics, 'genericStatusMessage', @output.genericStatusMessage);
        end
      elseif setup.enableElectronKinetics
        % add listener to status messages of the electron kinetics object
        addlistener(setup.electronKinetics, 'genericStatusMessage', @output.genericStatusMessage);
        % add listener to output results when a new solution for the EEDF is found
        addlistener(setup.electronKinetics, 'obtainedNewEedf', @output.electronKineticsSolution);
      end

    end
    
%  end
  
%  methods (Access = private)
    
    function saveSetupInfo(output, setupCellArray)
    % saveSetupInfo saves the setup of the current simulation
    
      fileName = [output.folder filesep 'setup.txt'];
      fileID = fopen(fileName, 'wt');
      for cell = setupCellArray
          fprintf(fileID, '%s\n', cell{1});
      end
      fclose(fileID);

    end % saveSetupInfo

    function saveInputFiles(output, setup)
    % saveInputFiles saves all the input files found in the setup of the simulation
    % inside an Input folder in the Output folder
      
      % find setup file
      files = {['Input' filesep setup.fileName]};

      % find electron kinetics input files
      if setup.enableElectronKinetics
        % find cross-section files (regular)
        for file = setup.info.electronKinetics.LXCatFiles
          files{end+1} = ['Input' filesep file{1}];
        end
        % find cross-section files (extra)
        if isfield(setup.info.electronKinetics, 'LXCatFilesExtra')
          for file = setup.info.electronKinetics.LXCatFilesExtra
            files{end+1} = ['Input' filesep file{1}];
          end
        end
        % find gas property files
        for field = fieldnames(setup.info.electronKinetics.gasProperties)'
          entries = setup.info.electronKinetics.gasProperties.(field{1});
          if ischar(entries)
            entries = {entries};
          end
          for entry = entries
            file = ['Input' filesep entry{1}];
            if isfile(file)
              files{end+1} = file;
            end
          end
        end
        % find state property files
        for field = fieldnames(setup.info.electronKinetics.stateProperties)'
          entries = setup.info.electronKinetics.stateProperties.(field{1});
          if ischar(entries)
            entries = {entries};
          end
          for entry = entries
            file = ['Input' filesep entry{1}];
            if isfile(file)
              files{end+1} = file;
            end
          end
        end
      end

      % find heavy-species kinetics input files
      if setup.enableChemistry
        % find chemistry files (files containing the reaction mechanism)
        for file = setup.info.chemistry.chemFiles
          files{end+1} = ['Input' filesep file{1}];
        end
        % find gas property files
        if isfield('gasProperties', setup.info.chemistry)
          for field = fieldnames(setup.info.chemistry.gasProperties)'
            entries = setup.info.chemistry.gasProperties.(field{1});
            if ischar(entries)
              entries = {entries};
            end
            for entry = entries
              file = ['Input' filesep entry{1}];
              if isfile(file)
                files{end+1} = file;
              end
            end
          end
        end
        % find state property files
        if isfield('stateProperties', setup.info.chemistry)
          for field = fieldnames(setup.info.chemistry.stateProperties)'
            entries = setup.info.chemistry.stateProperties.(field{1});
            if ischar(entries)
              entries = {entries};
            end
            for entry = entries
              file = ['Input' filesep entry{1}];
              if isfile(file)
                files{end+1} = file;
              end
            end
          end
        end
      end

      % create Input folder inside the current output folder
      inputOutputFolder = [output.folder filesep 'Input'];
      if ~isfolder(inputOutputFolder)
        mkdir(inputOutputFolder);
      end

      % copy input files to output folder
      for file = files
        finalFile = [output.folder filesep file{1}];
        [finalFolder, fileName, ~] = fileparts(finalFile);
        if ~isfolder(finalFolder)
          mkdir(finalFolder);
        end
        copyfile(file{1}, finalFile);
      end

    end
    
    function initializeLogFile(output, logCellArray)
    % initializeLogFile initialized the output file containing the log of the
    % simulation and writes previous messages of the log produced before the
    % creation of the output object
    
      fileName = [output.folder filesep 'log.txt'];
      fileID = fopen(fileName, 'wt');
      
      for cell = logCellArray
        fprintf(fileID, '%s\n', cell{1});
      end
      
      fclose(fileID);
      
    end
    
    function genericStatusMessage(output, ~, statusEventData)
      
      if output.logIsToBeSaved
        fileName = [output.folder filesep 'log.txt'];
        fileID = fopen(fileName, 'at');
        fprintf(fileID, statusEventData.message);
        fclose(fileID);
      end

    end

    function newPressureCycleIteration(output, chemistry, ~)

      if output.logIsToBeSaved
        fileName = [output.folder filesep 'log.txt'];
        fileID = fopen(fileName, 'at');
        fprintf(fileID, '\t- New pressure cycle iteration (%d): relative error = %e\n', ...
          chemistry.pressureIterationCurrent, chemistry.pressureRelErrorCurrent);
        fclose(fileID);
      end

    end

    function newNeutralityCycleIteration(output, chemistry, ~)

      if output.logIsToBeSaved
        fileName = [output.folder filesep 'log.txt'];
        fileID = fopen(fileName, 'at');
        fprintf(fileID, '\t- New neutrality cycle iteration (%d): relative error = %e\n', ...
          chemistry.neutralityIterationCurrent, chemistry.neutralityRelErrorCurrent);
        fclose(fileID);
      end

    end

    function newGlobalCycleIteration(output, chemistry, ~)

      if output.logIsToBeSaved
        fileName = [output.folder filesep 'log.txt'];
        fileID = fopen(fileName, 'at');
        fprintf(fileID, '\t- New global cycle iteration (%d): relative error = %e\n', ...
          chemistry.globalIterationCurrent, chemistry.globalRelErrorCurrent);
        fclose(fileID);
      end

    end

    function electronKineticsSolution(output, electronKinetics, ~)
    
      % create subfolder name in case of time-dependent boltzmann calculations
      if isa(electronKinetics, 'Boltzmann') && electronKinetics.isTimeDependent
        output.subFolder = sprintf('%s%s%stime_%e', filesep, output.subFolderBatches, filesep, electronKinetics.workCond.currentTime);
      end
      % create subfolder in case it is needed (when performing runs of simmulations or in time-dependent Boltzmann)
      if ~isempty(output.subFolder) && (output.eedfIsToBeSaved || output.powerBalanceIsToBeSaved || ...
          output.swarmParamsIsToBeSaved || output.rateCoeffsIsToBeSaved )
        if contains(output.dataFormat, 'txt')
          % By now output.folder SHOULD exist as per Output() function
          if ~isfolder([output.folder filesep output.subFolder])
            mkdir(output.folder,output.subFolder);
          end
        end
      end
      
      % if output format is hdf5 and isPulse, save the time and reducedField values
      if contains(output.dataFormat, 'hdf5')
        % DEBUG: Check if the E/N is written for 'chemistry'
        if output.isPulse
          % we only need to write the first output.numberOfEoverNJobs values...
          if output.currentJobID <= output.numberOfEoverNJobs
            data = [electronKinetics.workCond.currentTime, electronKinetics.workCond.reducedField];
            % write the values
            fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
            dseID = H5D.open(fID,'/electronKinetics/reducedField');
%             if output.isPulse
              start = [0 output.currentJobID-1];
%             else
%               start = [0 output.currentJobIndeces(1)-1];
%             end
            h5_block = [1 2];
            memSpaceID = H5S.create_simple(2,h5_block,[]);
            dspaceID = H5D.get_space(dseID);
            H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
            H5D.write(dseID,"H5ML_DEFAULT",memSpaceID,dspaceID,"H5P_DEFAULT",data);
            %
            H5S.close(dspaceID);
            H5S.close(memSpaceID);
            H5D.close(dseID);
            H5F.close(fID);
          end
        end
      end

      % save selected results of the electron kinetics
      if output.eedfIsToBeSaved
        if isa(electronKinetics, 'Boltzmann')
          % Check if 'Boltzmann' is a datatype. Otherwise change to strcmp function
          output.saveEedf(electronKinetics.eedf, electronKinetics.firstAnisotropy, electronKinetics.energyGrid.cell);
        else
          output.saveEedf(electronKinetics.eedf, [], electronKinetics.energyGrid.cell);
        end
      end
      if output.swarmParamsIsToBeSaved
        output.saveSwarm(electronKinetics.swarmParam, electronKinetics.workCond.reducedField, ...
          electronKinetics.workCond.electronDensity);
      end
      if output.rateCoeffsIsToBeSaved
        output.saveRateCoefficients(electronKinetics.rateCoeffAll, electronKinetics.rateCoeffExtra, []);
      end
      if output.powerBalanceIsToBeSaved
        output.savePower(electronKinetics.power);
      end
      if output.lookUpTablesAreToBeSaved
        if contains(output.dataFormat, 'txt')        % hdf5 format is already a lookUptable!
          output.saveLookUpTables(electronKinetics);
        end
      end
      
      output.currentJobID = output.currentJobID + 1;
    end
    
    function saveEedf(output, eedf, firstAnisotropy, energy)
    % saveEedf saves the eedf information of the current simulation
      
      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'eedf.txt'];

        % open file
        fileID = fopen(fileName, 'wt');

        % save information into the file
        if isempty(firstAnisotropy)
          fprintf(fileID, 'Energy(eV)           EEDF(eV^-(3/2))\n');
          values(2:2:2*length(eedf)) = eedf;
          values(1:2:2*length(eedf)) = energy;
          fprintf(fileID, '%#.14e %#.14e \n', values);
        else
          fprintf(fileID, 'Energy(eV)           EEDF(eV^-(3/2))      Anisotropy(eV^-(3/2))\n');
          values(3:3:3*length(eedf)) = firstAnisotropy;
          values(2:3:3*length(eedf)) = eedf;
          values(1:3:3*length(eedf)) = energy;
          fprintf(fileID, '%#.14e %#.14e %#.14e \n', values);
        end

        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        % write dataset on outputFile
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        dsfID = H5D.open(fID,'/electronKinetics/eedf');
        doubleType = H5T.copy('H5T_NATIVE_DOUBLE');
        offset(1)=0;
        data.Energy = energy;
        data.f0 = eedf;
        if ~isempty(firstAnisotropy)
          sz(1:3) = H5T.get_size(doubleType);
          offset(2:3)=cumsum(sz(1:2));
          data.f1 = firstAnisotropy;
          name = ["Energy" "EEDF" "Anisotropy"];
        else
          sz(1:2) = H5T.get_size(doubleType);
          offset(2:2)=cumsum(sz(1:1));
          name = ["Energy" "EEDF"];
        end
        memtype = H5T.create ('H5T_COMPOUND', sum(sz));
        for i = 1:length(sz)
          H5T.insert(memtype,name(i),offset(i),doubleType);
        end
        extraDims = length(output.currentJobIndeces);
        extraStart = output.currentJobIndeces - ones(1,extraDims);
        if output.isPulse
          reducedFieldStart = mod(output.currentJobID,output.numberOfEoverNJobs);
          if reducedFieldStart == 0
            reducedFieldStart = output.numberOfEoverNJobs;
          end
          start = [0 0 reducedFieldStart-1 extraStart(2:end)];
          block = [length(energy) 1 ones(1,extraDims)];
        else
          start = [0 0 extraStart];
          block = [length(energy) 1 ones(1,extraDims)];
        end
        h5_block = fliplr(block);
        memSpaceID = H5S.create_simple(length(block),h5_block,[]);
        dspaceID = H5D.get_space(dsfID);
        H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
        H5D.write(dsfID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
        % clean-up...
        H5T.close(doubleType);
        H5S.close(dspaceID);
        H5S.close(memSpaceID);
        H5T.close(memtype);
        H5D.close(dsfID);
        H5F.close(fID);
        clear sz;
      end
      
    end
    
    function saveSwarm(output, swarmParam, reducedField, electronDensity)
    % Saves the swarm parameters information of the current simulation
    
      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'swarmParameters.txt'];

        % open file
        fileID = fopen(fileName, 'wt');

        % save information into the file
        fprintf(fileID, '                     Electron density = %#.14e (m^-3)\n', electronDensity);
        fprintf(fileID, '               Reduced electric field = %#.14e (Td)\n', reducedField);
        fprintf(fileID, '                          Mean energy = %#.14e (eV)\n', swarmParam.meanEnergy);
        fprintf(fileID, '                Characteristic energy = %#.14e (eV)\n', swarmParam.characEnergy);
        fprintf(fileID, '                 Electron temperature = %#.14e (eV)\n', swarmParam.Te);
        if ~output.isSimulationHF(output.currentJobID)
          fprintf(fileID, '                       Drift velocity = %#.14e (ms^-1)\n', swarmParam.driftVelocity);
        end
        fprintf(fileID, '                     Reduced mobility = %#.14e ((msV)^-1)\n', swarmParam.redMobility);
        if output.isSimulationHF(output.currentJobID)
          fprintf(fileID, '                  Reduced mobility HF = %#.14e%+#.14ei ((msV)^-1)\n', ...
            real(swarmParam.redMobilityHF), imag(swarmParam.redMobilityHF));
        end
        fprintf(fileID, '        Reduced diffusion coefficient = %#.14e ((ms)^-1)\n', swarmParam.redDiffCoeff);
        fprintf(fileID, '              Reduced energy mobility = %#.14e (eV(msV)^-1)\n', swarmParam.redMobilityEnergy);
        fprintf(fileID, ' Reduced energy diffusion coefficient = %#.14e (eV(ms)^-1)\n', swarmParam.redDiffCoeffEnergy);
        if ~output.isSimulationHF(output.currentJobID)
          fprintf(fileID, '         Reduced Townsend coefficient = %#.14e (m^2)\n', swarmParam.redTownsendCoeff);
          fprintf(fileID, '       Reduced attachment coefficient = %#.14e (m^2)\n', swarmParam.redAttCoeff);
        end

        % close file
        fclose(fileID);
      end
      
      if contains(output.dataFormat, 'hdf5')
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        doubleType = H5T.copy('H5T_NATIVE_DOUBLE');
        dssID = H5D.open(fID,'/electronKinetics/swarmParameters');
        sz(1:9) = H5T.get_size(doubleType);
        offset(1)=0;
        % get offset and name
        if output.isSimulationHF(output.currentJobIndeces(1)) || ...
          output.isExcitationFrequencyBatch
          offset(2:9)=cumsum(sz(1:8));
          if output.isBoltzmann
            name = ["meanEnergy" "characEnergy" "Te" "redMobility" ...
                "redMobilityHFr" "redMobilityHFi" "redDiffCoeff" ...
                "redMobilityEnergy" "redDiffCoeffEnergy"];
            data.meanEnergy = swarmParam.meanEnergy;
            data.characEnergy = swarmParam.characEnergy;
            data.Te = swarmParam.Te;
          else
            name = ["meanEnergy" "characEnergy" "reducedField" "redMobility" ...
                "redMobilityHFr" "redMobilityHFi" "redDiffCoeff" ...
                "redMobilityEnergy" "redDiffCoeffEnergy"];
            data.meanEnergy = swarmParam.meanEnergy;
            data.characEnergy = swarmParam.characEnergy;
            data.reducedField = reducedField;
          end
          data.redMobility = swarmParam.redMobility;
          if ~isempty(swarmParam.redMobilityHF)
            data.redMobilityHFr = real(swarmParam.redMobilityHF);
            data.redMobilityHFi = imag(swarmParam.redMobilityHF);
          else
            data.redMobilityHFr = 0.0;
            data.redMobilityHFi = 0.0;
          end
          data.redDiffCoeff = swarmParam.redDiffCoeff;
          data.redMobilityEnergy = swarmParam.redMobilityEnergy;
          data.redDiffCoeffEnergy = swarmParam.redDiffCoeffEnergy;
        else
          sz(10) = H5T.get_size(doubleType);
          offset(2:10)=cumsum(sz(1:9));
          if output.isBoltzmann
            name = ["meanEnergy" "characEnergy" "Te" "driftVelocity" ...
                "redMobility" "redDiffCoeff" "redMobilityEnergy" ...
                "redDiffCoeffEnergy" "redTownsendCoeff" "redAttCoeff"];
            data.meanEnergy = swarmParam.meanEnergy;
            data.characEnergy = swarmParam.characEnergy;
            data.Te = swarmParam.Te;
          else
            name = ["meanEnergy" "characEnergy" "reducedField" "driftVelocity" ...
                "redMobility" "redDiffCoeff" "redMobilityEnergy" ...
                "redDiffCoeffEnergy" "redTownsendCoeff" "redAttCoeff"];
            data.meanEnergy = swarmParam.meanEnergy;
            data.characEnergy = swarmParam.characEnergy;
            data.reducedField = reducedField;
          end   
          data.driftVelocity = swarmParam.driftVelocity;
          data.redMobility = swarmParam.redMobility;
          data.redDiffCoeff = swarmParam.redDiffCoeff;
          data.redMobilityEnergy = swarmParam.redMobilityEnergy;
          data.redDiffCoeffEnergy = swarmParam.redDiffCoeffEnergy;
          data.redTownsendCoeff = swarmParam.redTownsendCoeff;
          data.redAttCoeff = swarmParam.redAttCoeff;
        end
        memtype = H5T.create ('H5T_COMPOUND', sum(sz));
        for i = 1:length(sz)
          H5T.insert(memtype,name(i),offset(i),doubleType);
        end
        if output.isPulse
          extraStart = output.currentJobIndeces(2:end) - ones(1,length(output.extraDims));
          reducedFieldStart = mod(output.currentJobID,output.numberOfEoverNJobs);
          if reducedFieldStart == 0
            reducedFieldStart = output.numberOfEoverNJobs;
          end
          start = [reducedFieldStart-1 0 extraStart];
          h5_block = [1 1 ones(1,length(output.extraDims))];
        else
          extraStart = output.currentJobIndeces(2:end) - ones(1,length(output.extraDims));
          start = [output.currentJobIndeces(1)-1 0 extraStart];  % Note: location is 0-based, not 1-based!
          h5_block = [1 1 ones(1,length(output.extraDims))];
        end
        memSpaceID = H5S.create_simple(length(h5_block),h5_block,[]);
        dspaceID = H5D.get_space(dssID);
        H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
        H5D.write(dssID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
        % clean-up...
        H5S.close(dspaceID);
        H5S.close(memSpaceID);
        H5T.close(memtype);
        H5T.close(doubleType);
        H5D.close(dssID);
        H5F.close(fID);
        clear sz;
      end
    end
    
    function saveRateCoefficients(output, eKineticsRateCoeffs, eKineticsRateCoeffsExtra, reactionsInfo)
    % saveRateCoefficients saves the rate coefficients obtained in the current simulation
      
      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'rateCoefficients.txt'];
      
        % open file
        fileID = fopen(fileName, 'wt');
      
        % save information into the file
        if ~isempty(eKineticsRateCoeffs)
          fprintf(fileID, '%s\n*    e-Kinetics Rate Coefficients    *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
          fprintf(fileID, 'ID   Ine.R.Coeff.(m^3s^-1) Sup.R.Coeff.(m^3s^-1) Threshold(eV)         Description\n');
          for rateCoeff = eKineticsRateCoeffs
            if length(rateCoeff.value) == 1
              fprintf(fileID, '%4d %20.14e  (N/A)                 %20.14e  %s\n', rateCoeff.collID, rateCoeff.value, ...
                rateCoeff.energy, rateCoeff.collDescription);
            else
              fprintf(fileID, '%4d %20.14e  %20.14e  %20.14e  %s\n', rateCoeff.collID, rateCoeff.value(1), ...
                rateCoeff.value(2), rateCoeff.energy, rateCoeff.collDescription);
            end
          end
        end
        if ~isempty(eKineticsRateCoeffsExtra)
          fprintf(fileID, '\n%s\n* e-Kinetics Extra Rate Coefficients *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
          fprintf(fileID, 'ID   Ine.R.Coeff.(m^3s^-1) Sup.R.Coeff.(m^3s^-1) Threshold(eV)         Description\n');
          for rateCoeff = eKineticsRateCoeffsExtra
            if length(rateCoeff.value) == 1
              fprintf(fileID, '%4d %20.14e  (N/A)                 %20.14e  %s\n', rateCoeff.collID, rateCoeff.value, ...
                rateCoeff.energy, rateCoeff.collDescription);
            else
              fprintf(fileID, '%4d %20.14e  %20.14e  %20.14e  %s\n', rateCoeff.collID, rateCoeff.value(1), ...
                rateCoeff.value(2), rateCoeff.energy, rateCoeff.collDescription);
            end
          end
        end
        if ~isempty(reactionsInfo)
          fprintf(fileID, '\n%s\n*     Chemistry Rate Coefficients    *\n%s\n\n', repmat('*', 1,38), repmat('*', 1,38));
          fprintf(fileID, ['ID   Dir.R.Coeff.(S.I.)    Inv.R.Coeff.(S.I.)    Enthalpy(eV)          ' ...
            'Net.Reac.Rate(m^-3s^-1) Description\n']);
          for reaction = reactionsInfo
            if length(reaction.rateCoeff) == 1
              fprintf(fileID, '%4d %20.14e  (N/A)                 %+20.14e %+20.14e   %s\n', reaction.reactID, ...
                reaction.rateCoeff, reaction.energy, reaction.netRate, reaction.description);
            else
              fprintf(fileID, '%4d %20.14e  %20.14e  %+20.14e %+20.14e   %s\n', reaction.reactID, ...
                reaction.rateCoeff(1), reaction.rateCoeff(2), reaction.energy, reaction.netRate, reaction.description);
            end
          end
        end
        fclose(fileID);
      end
      
      if contains(output.dataFormat, 'hdf5')
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        % Create the base types that will be used is the datasets
        intType     = H5T.copy('H5T_NATIVE_INT');
        doubleType  = H5T.copy('H5T_NATIVE_DOUBLE');
        strType     = H5T.copy ('H5T_C_S1');
        H5T.set_size(strType, 'H5T_VARIABLE');
        %
        if ~isempty(eKineticsRateCoeffs)
          % process the hdf5 file
          dsrID = H5D.open(fID,'/electronKinetics/rateCoefficients');
          % Create the base types
          sz(1) = H5T.get_size(intType);
          sz(2:4) = H5T.get_size(doubleType);
          sz(5) = H5T.get_size(strType);
          % Compute the offsets to each field. The first offset is always zero.
          offset(1)=0;
          offset(2:5)=cumsum(sz(1:4));
          % Create the compound datatype for memory.
          memtype = H5T.create ('H5T_COMPOUND', sum(sz));
          H5T.insert(memtype,'rate_id',offset(1),intType);
          H5T.insert(memtype,'ine_coeff',offset(2),doubleType);
          H5T.insert(memtype,'sup_coeff',offset(3),doubleType);
          H5T.insert(memtype,'threshold',offset(4),doubleType);
          H5T.insert(memtype,'description',offset(5),strType);
          % Get the data values
          ratePosition = -1;
          for rateCoeff = eKineticsRateCoeffs
            ratePosition = ratePosition + 1;
            data.rate_id    = int32(rateCoeff.collID);
            if length(rateCoeff.value) == 1
              data.ine_coeff = rateCoeff.value;
              data.sup_coeff = 0.0;
            else
              data.ine_coeff = rateCoeff.value(1);
              data.sup_coeff = rateCoeff.value(2);
            end
            data.threshold  = rateCoeff.energy;
            data.reaction   = rateCoeff.collDescription;
            %
            if output.isPulse
              extraDims = length(output.currentJobIndeces);
              extraStart = output.currentJobIndeces - ones(1,extraDims);
              reducedFieldStart = mod(output.currentJobID,output.numberOfEoverNJobs);
              if reducedFieldStart == 0
                reducedFieldStart = output.numberOfEoverNJobs;
              end
              start = [ratePosition 0 reducedFieldStart-1 extraStart(2:end)];
%               start = [ratePosition 0 extraStart];
              h5_block = [1 1 ones(1,extraDims)];
            else
              extraDims = length(output.currentJobIndeces);
              extraStart = output.currentJobIndeces - ones(1,extraDims);
              start = [ratePosition 0 extraStart];
              h5_block = [1 1 ones(1,extraDims)];
            end
            memSpaceID = H5S.create_simple(length(h5_block),h5_block,[]);
            dspaceID = H5D.get_space(dsrID);
            H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
            H5D.write(dsrID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
          end
          clear data;
          H5S.close(dspaceID);
          H5D.close(dsrID);
          % We don't close memtype and memSpaceID as we may need them for eKineticsRateCoeffsExtra
          if ~isempty(eKineticsRateCoeffsExtra)
            % process the extra e-collisions rates
            dsrxID = H5D.open(fID,'/electronKinetics/extraRateCoefficients');
            ratePosition = -1;
            for rateCoeff = eKineticsRateCoeffsExtra
              ratePosition = ratePosition + 1;
              if length(rateCoeff.value) == 1
                data.rate_id = int32(rateCoeff.collID);
                data.ine_coeff = rateCoeff.value;
                data.sup_coeff = 0.0;
                data.threshold = rateCoeff.energy;
                data.reaction = rateCoeff.collDescription;
              else
                data.rate_id = int32(rateCoeff.collID);
                data.ine_coeff = rateCoeff.value(1);
                data.sup_coeff = rateCoeff.value(2);
                data.threshold = rateCoeff.energy;
                data.reaction = rateCoeff.collDescription;
              end
              extraDims = length(output.currentJobIndeces);
              extraStart = output.currentJobIndeces - ones(1,extraDims);
              start = [ratePosition 0 extraStart];
              h5_block = [1 1 ones(1,extraDims)];
%               memSpaceID = H5S.create_simple(length(h5_block),h5_block,[]);  % We use the one defined for eKineticsRateCoeffs
              dspaceID = H5D.get_space(dsrxID);
              H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
              H5D.write(dsrxID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
            end
            H5S.close(dspaceID);
            H5D.close(dsrxID);
            clear data sz;
          end
          H5S.close(memSpaceID);  % We close memSpaceID and memtype open in ...
          H5T.close(memtype);     % eKineticsRateCoeffs
        end
        %
        if ~isempty(reactionsInfo)
          % rate coefficients for the chemistry
          dsrID = H5D.open(fID,'/chemistry/rateCoefficients');
          % Create the base types
          sz(1) = H5T.get_size(intType);
          sz(2:5) = H5T.get_size(doubleType);
          sz(6) = H5T.get_size(strType);
          % Compute the offsets to each field. The first offset is always zero.
          offset(1)=0;
          offset(2:6)=cumsum(sz(1:5));
          % Create the compound datatype for memory.
          memtype = H5T.create ('H5T_COMPOUND', sum(sz));
          H5T.insert(memtype,'rate_id',offset(1),intType);
          H5T.insert(memtype,'dir_coeff',offset(2),doubleType);
          H5T.insert(memtype,'inv_coeff',offset(3),doubleType);
          H5T.insert(memtype,'enthalpy',offset(4),doubleType);
          H5T.insert(memtype,'net_rate',offset(5),doubleType);
          H5T.insert(memtype,'description',offset(6),strType);
          % Get the data values
          ratePosition = -1;
          for reaction = reactionsInfo
            ratePosition    = ratePosition + 1;
            data.rate_id    = int32(reaction.reactID);
            if length(reaction.rateCoeff) == 1
              data.dir_coeff  = reaction.rateCoeff;
              data.inv_coeff  = 0.0;
            else
              data.dir_coeff  = reaction.rateCoeff(1);
              data.inv_coeff  = reaction.rateCoeff(2);
            end
            data.enthalpy   = reaction.energy;
            data.net_rate  = reaction.netRate;
            data.description  = reaction.description;
            % get dims, spaceID and writes dataset
            start = [ratePosition 0];
            h5_block = [1 1];
            memSpaceID = H5S.create_simple(2,h5_block,[]);
            dspaceID = H5D.get_space(dsrID);
            H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
            H5D.write(dsrID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
          end
          H5D.close(dsrID);
          H5T.close(memtype);
          H5S.close(dspaceID);
          H5S.close(memSpaceID);
          H5F.close(fID);
        end
      end

    end
    
    function savePower(output, power)
    % savePower saves the power balance information of the current simulation
      
      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'powerBalance.txt'];
        
        % open file
        fileID = fopen(fileName, 'wt');
        
        % save information into the file
        fprintf(fileID, '                               Field = %#+.14e (eVm^3s^-1)\n', power.field);
        fprintf(fileID, '           Elastic collisions (gain) = %#+.14e (eVm^3s^-1)\n', power.elasticGain);
        fprintf(fileID, '           Elastic collisions (loss) = %#+.14e (eVm^3s^-1)\n', power.elasticLoss);
        fprintf(fileID, '                          CAR (gain) = %#+.14e (eVm^3s^-1)\n', power.carGain);
        fprintf(fileID, '                          CAR (loss) = %#+.14e (eVm^3s^-1)\n', power.carLoss);
        fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationIne);
        fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationSup);
        fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalIne);
        fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalSup);
        fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalIne);
        fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalSup);
        fprintf(fileID, '               Ionization collisions = %#+.14e (eVm^3s^-1)\n', power.ionizationIne);
        fprintf(fileID, '               Attachment collisions = %#+.14e (eVm^3s^-1)\n', power.attachmentIne);
        fprintf(fileID, '             Electron density growth = %#+.14e (eVm^3s^-1) +\n', power.eDensGrowth);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '                       Power Balance = %#+.14e (eVm^3s^-1)\n', power.balance);
        fprintf(fileID, '              Relative Power Balance = % #.14e%%\n\n', power.relativeBalance*100);
        fprintf(fileID, '           Elastic collisions (gain) = %#+.14e (eVm^3s^-1)\n', power.elasticGain);
        fprintf(fileID, '           Elastic collisions (loss) = %#+.14e (eVm^3s^-1) +\n', power.elasticLoss);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '            Elastic collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.elasticNet);
        fprintf(fileID, '                          CAR (gain) = %#+.14e (eVm^3s^-1)\n', power.carGain);
        fprintf(fileID, '                          CAR (loss) = %#+.14e (eVm^3s^-1) +\n', power.carLoss);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '                           CAR (net) = %#+.14e (eVm^3s^-1)\n\n', power.carNet);
        fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.excitationIne);
        fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.excitationSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '         Excitation collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.excitationNet);
        fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.vibrationalIne);
        fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.vibrationalSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '        Vibrational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', power.vibrationalNet);
        fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', power.rotationalIne);
        fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', power.rotationalSup);
        fprintf(fileID, ' %s\n', repmat('-', 1, 73));
        fprintf(fileID, '         Rotational collisions (net) = %#+.14e (eVm^3s^-1)\n', power.rotationalNet);

        % power balance by gases
        gases = fields(power.gases);
        powerByGas = power.gases;
        for i = 1:length(gases)
          gas = gases{i};
          fprintf(fileID, '\n%s\n\n', [repmat('*', 1, 37) ' ' gas ' ' repmat('*', 1, 39-length(gas))]);
          fprintf(fileID, '     Excitation inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).excitationIne);
          fprintf(fileID, '  Excitation superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).excitationSup);
          fprintf(fileID, ' %s\n', repmat('-', 1, 73));
          fprintf(fileID, '         Excitation collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).excitationNet);
          fprintf(fileID, '    Vibrational inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).vibrationalIne);
          fprintf(fileID, ' Vibrational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).vibrationalSup);
          fprintf(fileID, ' %s\n', repmat('-', 1, 73));
          fprintf(fileID, '        Vibrational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).vibrationalNet);
          fprintf(fileID, '     Rotational inelastic collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).rotationalIne);
          fprintf(fileID, '  Rotational superelastic collisions = %#+.14e (eVm^3s^-1) +\n', powerByGas.(gas).rotationalSup);
          fprintf(fileID, ' %s\n', repmat('-', 1, 73));
          fprintf(fileID, '         Rotational collisions (net) = %#+.14e (eVm^3s^-1)\n\n', powerByGas.(gas).rotationalNet);
          fprintf(fileID, '               Ionization collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).ionizationIne);
          fprintf(fileID, '               Attachment collisions = %#+.14e (eVm^3s^-1)\n', powerByGas.(gas).attachmentIne);
        end
        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        % Convert data in struct to mat
        temp = struct2cell(power);
        dataSummary = [cell2mat(temp(1:24)); cell2mat(temp(26:end))];
        % write powerBalanceSummary
        doubleType = H5T.copy("H5T_NATIVE_DOUBLE");
        sz(1:10) = H5T.get_size(doubleType);
        offset(1) = 0;
        offset(2:10) = cumsum(sz(1:9));
        name = ["Field" "Elastic" "CAR" "Rotational" "Vibrational" "Electronic" ...
          "Ionization" "Attachment" "eDensGrowth" "Balance"];
        memtype = H5T.create('H5T_COMPOUND', sum(sz));
        for i = 1:length(sz)
          H5T.insert(memtype,name(i),offset(i),doubleType);
        end
        powerSummary.field = [dataSummary(1) dataSummary(1) 0];                        % field
        powerSummary.elastic = dataSummary(2:4);                                       % elastic
        powerSummary.CAR = dataSummary(5:7);                                           % CAR
        powerSummary.rotational = [dataSummary(16) dataSummary(15) dataSummary(14)];   % rotational
        powerSummary.vibrational = [dataSummary(13) dataSummary(12) dataSummary(11)];  % vibrational
        powerSummary.electronic = [dataSummary(10) dataSummary(9) dataSummary(8)];     % electronic
        powerSummary.ionization = [dataSummary(17) 0 dataSummary(17)];                 % ionization
        powerSummary.attachment = [dataSummary(18) 0 dataSummary(18)];                 % attachment
        powerSummary.eDensGrowth = [dataSummary(21) dataSummary(21) 0];                % eDensGrowth
        powerSummary.balance = [dataSummary(25) dataSummary(26) 0];                    % balance and relativeBalance
        % Process the hdf5 file
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        % powerBalanceSummary
        dspID = H5D.open(fID,'/electronKinetics/powerBalanceSummary');
        extraStart = output.currentJobIndeces(2:end) - ones(1,length(output.extraDims));
        if output.isPulse
          reducedFieldStart = mod(output.currentJobID,output.numberOfEoverNJobs);
          if reducedFieldStart == 0
            reducedFieldStart = output.numberOfEoverNJobs;
          end
          start = [reducedFieldStart-1 0 0 extraStart];
%           start = [output.currentJobIndeces(1)-1 0 0 extraStart];
        else
          start = [output.currentJobIndeces(1)-1 0 0 extraStart];
        end
        block = [1 1 3 ones(1,length(output.extraDims))];
        h5_block = fliplr(block);
        memSpaceID = H5S.create_simple(length(block),h5_block,[]);
        dspaceID = H5D.get_space(dspID);
        H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
        H5D.write(dspID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",powerSummary);
        H5S.close(memSpaceID);
        H5D.close(dspID);
        H5T.close(memtype);
        % powerBalanceGases
        dspID = H5D.open(fID,'/electronKinetics/powerBalanceGases');
        gases = fields(power.gases);
        clear sz;
        sz(1:5) = H5T.get_size(doubleType);
        offset(1) = 0;
        offset(2:5) = cumsum(sz(1:4));
        name = ["rotCol" "vibCol" "eleCol" "ionCol" "attCol"];
        memtype = H5T.create('H5T_COMPOUND', sum(sz));
        for i = 1:length(sz)
          H5T.insert(memtype,name(i),offset(i),doubleType);
        end
        for i = 1:length(gases)
          gas = gases{i};
          temp = struct2cell(power.gases.(gas));
          powerGas.rot = [cell2mat(temp(9)) cell2mat(temp(8)) cell2mat(temp(7))];
          powerGas.vib = [cell2mat(temp(6)) cell2mat(temp(5)) cell2mat(temp(4))];
          powerGas.ele = [cell2mat(temp(1)) cell2mat(temp(2)) cell2mat(temp(1))];
          powerGas.ion = [cell2mat(temp(10)) 0 cell2mat(temp(10))];
          powerGas.att = [cell2mat(temp(11)) 0 cell2mat(temp(11))];
          % DEBUG:
          % powerGas also includes inelastic and superelastic fields,
          % these fiels are NOT included in the initial definition!
          if output.isPulse
            reducedFieldStart = mod(output.currentJobID,output.numberOfEoverNJobs);
            if reducedFieldStart == 0
              reducedFieldStart = output.numberOfEoverNJobs;
            end
            start = [reducedFieldStart-1 0 0 extraStart i-1];
          else
            start = [output.currentJobIndeces(1)-1 0 0 extraStart i-1];
          end
          block = [1 1 3 ones(1,length(output.extraDims)) 1];
          h5_block = fliplr(block);
          memSpaceID = H5S.create_simple(length(block),h5_block,[]);
          dspaceID = H5D.get_space(dspID);
          H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",fliplr(start),[],[],h5_block);
          H5D.write(dspID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",powerGas);
        end
        H5D.close(dspID);
        H5F.close(fID);
      end

    end
    
    function saveLookUpTables(output, electronKinetics)
    % NOTE: lookUpTables are only created if output.dataFormat is 'txt'

      % name of the files containing the different lookup tables
      persistent fileName1;
      persistent fileName2;
      persistent fileName3;
      persistent fileName4;
      persistent fileName5;
      persistent folderLookUpTables;
      
      % local copies of different variables (for performance reasons)
      workCond = electronKinetics.workCond;
      power = electronKinetics.power;
      swarmParams = electronKinetics.swarmParam;
      rateCoeffAll = electronKinetics.rateCoeffAll;
      rateCoeffExtra = electronKinetics.rateCoeffExtra;
      eedf = electronKinetics.eedf;

      if isempty(folderLookUpTables) && ~isempty(output.subFolderBatches) 
        folderLookUpTables = output.subFolderBatches;
      end
      localFolderLookUpTables = output.subFolderBatches;
      
      % initialize the files in case it is needed
      if (isempty(fileName1) || ~strcmp(folderLookUpTables,localFolderLookUpTables)) 
        folderLookUpTables = output.subFolderBatches;
        % create file names
        fileName1 = [output.folder output.subFolderBatches filesep 'lookUpTableSwarm.txt'];
        fileName2 = [output.folder output.subFolderBatches filesep 'lookUpTablePower.txt'];
        fileName3 = [output.folder output.subFolderBatches filesep 'lookUpTableRateCoeff.txt'];
        % open files
        fileID1 = fopen(fileName1, 'wt');
        fileID2 = fopen(fileName2, 'wt');
        fileID3 = fopen(fileName3, 'wt');
        % write file headers
        fprintf(fileID3, [repmat('#', 1, 80) '\n# %-76s #\n'], 'ID   Description');
        strFile3 = '';
        for i = 1:length(rateCoeffAll)
          fprintf(fileID3, '# %-4d %-71s #\n', rateCoeffAll(i).collID, rateCoeffAll(i).collDescription);
          strAux = sprintf('R%d_ine(m^3s^-1)', rateCoeffAll(i).collID);
          strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          if 2 == length(rateCoeffAll(i).value)
            strAux = sprintf('R%d_sup(m^3s^-1)', rateCoeffAll(i).collID);
            strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          end
        end
        fprintf(fileID3, '#%s#\n# %-76s #\n#%s#\n# %-76s #\n', repmat(' ', 1, 78), ...
          '*** Extra rate coefficients ***', repmat(' ', 1, 78), 'ID   Description');
        for i = 1:length(rateCoeffExtra)
          fprintf(fileID3, '# %-4d %-71s #\n', rateCoeffExtra(i).collID, rateCoeffExtra(i).collDescription);
          strAux = sprintf('R%d_ine(m^3s^-1)', rateCoeffExtra(i).collID);
          strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          if 2 == length(rateCoeffExtra(i).value)
            strAux = sprintf('R%d_sup(m^3s^-1)', rateCoeffExtra(i).collID);
            strFile3 = sprintf('%s%-21s ', strFile3, strAux);
          end
        end
        fprintf(fileID3, [repmat('#', 1, 80) '\n\n']);
        if isa(electronKinetics, 'Boltzmann') 
          if electronKinetics.isTimeDependent
            fprintf(fileID1, '%-21s ', 'Time(s)');
            fprintf(fileID2, '%-21s ', 'Time(s)');
            fprintf(fileID3, '%-21s ', 'Time(s)');
            % create lookup table for the eedf
            fileName4 = [output.folder filesep output.subFolderBatches filesep 'lookUpTableEedf.txt'];
            fileID4 = fopen(fileName4, 'wt');
            % add first line with energies to eedf lookup table (eedfs will be saved as rows)
            fprintf(fileID4, '%-21.14e ', [0 electronKinetics.energyGrid.cell]);
            fprintf(fileID4, '\n');
            fclose(fileID4);
            % create lookup table for the electron density (if needed)
            if electronKinetics.eDensIsTimeDependent
             fileName5 = [output.folder filesep output.subFolderBatches filesep 'lookUpTableElectronDensity.txt'];
             fileID5 = fopen(fileName5, 'wt');
              fprintf(fileID5, '%-21s %-21s\n', 'time(s)', 'ne(m^-3)\n');
              fclose(fileID5);
            end
          end
          if output.isSimulationHF(output.currentJobID)
            fprintf(fileID1, [repmat('%-21s ', 1, 10) '\n'], 'RedField(Td)', 'RedDiff((ms)^-1)', 'RedMob((msV)^-1)', ...
              'R[RedMobHF]((msV)^-1)', 'I[RedMobHF]((msV)^-1)', 'RedDiffE(eV(ms)^-1)', 'RedMobE(eV(msV)^-1)', ...
              'MeanE(eV)', 'CharE(eV)', 'EleTemp(eV)');
          else
            fprintf(fileID1, [repmat('%-21s ', 1, 11) '\n'], 'RedField(Td)', 'RedDiff((ms)^-1)', 'RedMob((msV)^-1)', ...
              'DriftVelocity(ms^-1)', 'RedTow(m^2)', 'RedAtt(m^2)', 'RedDiffE(eV(ms)^-1)', 'RedMobE(eV(msV)^-1)', ...
              'MeanE(eV)', 'CharE(eV)', 'EleTemp(eV)');
          end
          fprintf(fileID2, '%-21s ', 'RedField(Td)');
          fprintf(fileID3, '%-21s ', 'RedField(Td)');
        else
          if output.isSimulationHF(output.currentJobID)
            fprintf(fileID1, [repmat('%-21s ', 1, 10) '\n'], 'EleTemp(eV)', 'RedField(Td)', 'RedDiff(1/(ms))', ...
              'RedMob(1/(msV))', 'R[RedMobHF](1/(msV))', 'I[RedMobHF](1/(msV))', 'RedDiffE(eV/(ms))', ...
              'RedMobE(eV/(msV))', 'MeanE(eV)', 'CharE(eV)');
          else
            fprintf(fileID1, [repmat('%-21s ', 1, 11) '\n'], 'EleTemp(eV)', 'RedField(Td)', 'RedDiff(1/(ms))', ...
            'RedMob(1/(msV))', 'DriftVelocity(m/s)', 'RedTow(m2)', 'RedAtt(m2)', 'RedDiffE(eV/(ms))', 'RedMobE(eV/(msV))', 'MeanE(eV)', ...
            'CharE(eV)');
          end
          fprintf(fileID2, '%-21s ', 'EleTemp(eV)');
          fprintf(fileID3, '%-21s ', 'EleTemp(eV)');
        end
        fprintf(fileID2, [repmat('%-21s ', 1, 21) '\n'], 'PowerField(eVm^3s^-1)', ...
          'PwrElaGain(eVm^3s^-1)', 'PwrElaLoss(eVm^3s^-1)', 'PwrElaNet(eVm^3s^-1)', 'PwrCARGain(eVm^3s^-1)', ...
          'PwrCARLoss(eVm^3s^-1)', 'PwrCARNet(eVm^3s^-1)', 'PwrEleGain(eVm^3s^-1)', 'PwrEleLoss(eVm^3s^-1)', ...
          'PwrEleNet(eVm^3s^-1)', 'PwrVibGain(eVm^3s^-1)', 'PwrVibLoss(eVm^3s^-1)', 'PwrVibNet(eVm^3s^-1)', ...
          'PwrRotGain(eVm^3s^-1)', 'PwrRotLoss(eVm^3s^-1)', 'PwrRotNet(eVm^3s^-1)', 'PwrIon(eVm^3s^-1)', ...
          'PwrAtt(eVm^3s^-1)', 'PwrGroth(eVm^3s^-1)', 'PwrBalance(eVm^3s^-1)', 'RelPwrBalance');
        fprintf(fileID3, '%s\n', strFile3);
        % close files
        fclose(fileID1);
        fclose(fileID2);
        fclose(fileID3);
      end
      
      % check if eedf lookup table needs to be saved (and append new line with data)
      if ~isempty(fileName4)
        fileID4 = fopen(fileName4, 'at');
        fprintf(fileID4, '%-21.14e ', workCond.currentTime);
        fprintf(fileID4, '%-21.14e ', eedf);
        fprintf(fileID4, '\n');
        fclose(fileID4);
      end
      % check if electron density data needs to be saved (and append new line with data)
      if ~isempty(fileName5)
        fileID5 = fopen(fileName5, 'at');
        fprintf(fileID5, '%#.14e %#.14e\n',workCond.currentTime, workCond.electronDensity);
        fclose(fileID5);
      end

      % open files
      fileID1 = fopen(fileName1, 'at');
      fileID2 = fopen(fileName2, 'at');
      fileID3 = fopen(fileName3, 'at');
      % append new lines with data
      if isa(electronKinetics, 'Boltzmann')
        if electronKinetics.isTimeDependent
          fprintf(fileID1, '%-+21.14e ', workCond.currentTime);
          fprintf(fileID2, '%-+21.14e ', workCond.currentTime);
          fprintf(fileID3, '%-+21.14e ', workCond.currentTime);
        end
        if output.isSimulationHF(output.currentJobID)
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 10) '\n'], ...
            workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
            real(swarmParams.redMobilityHF), imag(swarmParams.redMobilityHF), swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy, swarmParams.Te);
        else
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 11) '\n'], ...
            workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, swarmParams.driftVelocity, ...
            swarmParams.redTownsendCoeff, swarmParams.redAttCoeff, swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy, swarmParams.Te);
        end
        fprintf(fileID2, '%-+21.14e ', workCond.reducedField);
        fprintf(fileID3, '%-+21.14e ', workCond.reducedField);
      else
        if output.isSimulationHF(output.currentJobID)
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 10) '\n'], ...
            swarmParams.Te, workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
             real(swarmParams.redMobilityHF), imag(swarmParams.redMobilityHF), swarmParams.redDiffCoeffEnergy, ...
            swarmParams.redMobilityEnergy, swarmParams.meanEnergy, swarmParams.characEnergy);
        else
          fprintf(fileID1, [repmat('%-+21.14e ', 1, 11) '\n'], ...
            swarmParams.Te, workCond.reducedField, swarmParams.redDiffCoeff, swarmParams.redMobility, ...
            swarmParams.driftVelocity, swarmParams.redTownsendCoeff, swarmParams.redAttCoeff, ...
            swarmParams.redDiffCoeffEnergy, swarmParams.redMobilityEnergy, swarmParams.meanEnergy, ...
            swarmParams.characEnergy);
        end
        fprintf(fileID2, '%-+21.14e ', workCond.electronTemperature);
        fprintf(fileID3, '%-+21.14e ', workCond.electronTemperature);
      end
      fprintf(fileID2, [repmat('%-+21.14e ', 1, 20) '%19.14e%%\n'], power.field, ...
        power.elasticGain, power.elasticLoss, power.elasticNet, power.carGain, power.carLoss, power.carNet, ...
        power.excitationSup, power.excitationIne, power.excitationNet, power.vibrationalSup, power.vibrationalIne, ...
        power.vibrationalNet, power.rotationalSup, power.rotationalIne, power.rotationalNet, power.ionizationIne, ...
        power.attachmentIne, power.eDensGrowth, power.balance, power.relativeBalance*100);
      for i = 1:length(rateCoeffAll)
        fprintf(fileID3, '%-21.14e ', rateCoeffAll(i).value(1));
        if 2 == length(rateCoeffAll(i).value)
          fprintf(fileID3, '%-21.14e ', rateCoeffAll(i).value(2));
        end
      end
      for i = 1:length(rateCoeffExtra)
        fprintf(fileID3, '%-21.14e ', rateCoeffExtra(i).value(1));
        if 2 == length(rateCoeffExtra(i).value)
          fprintf(fileID3, '%-21.14e ', rateCoeffExtra(i).value(2));
        end
      end
      fprintf(fileID3, '\n');
      % close files
      fclose(fileID1);
      fclose(fileID2);
      fclose(fileID3);
      
    end
    
    function chemistrySolution(output, chemistry, ~) 

      % create subfolder name in case of pulsed chemistry calculations
      if chemistry.isPulsed
        output.subFolder = sprintf('%s%s%stime_%e', filesep, output.subFolderBatches, filesep, ...
          chemistry.workCond.currentTime);
      end

      if contains(output.dataFormat, 'txt')
        % create subfolder in case it is needed (when performing runs of simmulations or pulsed simulations)
        if ~isempty(output.subFolder) && (output.eedfIsToBeSaved || output.powerBalanceIsToBeSaved || ...
            output.swarmParamsIsToBeSaved || output.rateCoeffsIsToBeSaved || output.finalDensitiesIsToBeSaved || ...
            output.finalTemperaturesIsToBeSaved || output.finalParticleBalanceIsToBeSaved || ...
            output.finalThermalBalanceIsToBeSaved || output.chemParamsIsToBeSaved)
          if 7 ~= exist([output.folder filesep output.subFolder], 'file')
            mkdir([output.folder filesep output.subFolder]);
          end
        end
%       elseif strcmp(output.dataFormat, 'hdf5') % in this case we don't need extra folders
      end
              
      % save results of the last electron Kinetics solution (in case it is activated)
      if ~isempty(chemistry.electronKinetics)
        electronKinetics = chemistry.electronKinetics;
        if output.eedfIsToBeSaved
          if isa(electronKinetics, 'Boltzmann')
            output.saveEedf(electronKinetics.eedf, electronKinetics.firstAnisotropy, electronKinetics.energyGrid.cell);
          else
            output.saveEedf(electronKinetics.eedf, [], electronKinetics.energyGrid.cell);
          end
        end
        if output.powerBalanceIsToBeSaved
          output.savePower(electronKinetics.power);
        end
        if output.swarmParamsIsToBeSaved
          if contains(output.dataFormat, "hdf5")
            % Save equilibrium value for the reducedField
            fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
            dims = [1 1];
            spaceID = H5S.create_simple(2,dims,[]);
            if isa(electronKinetics, 'Boltzmann')
              dsfID = H5D.open(fID,'/electronKinetics/reducedField');
              H5D.write(dsfID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', ...
                'H5P_DEFAULT',chemistry.workCond.reducedField);
            else
              dsfID = H5D.open(fID,'/electronKinetics/electronTemperature');
              H5D.write(dsfID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', ...
                'H5P_DEFAULT',chemistry.workCond.electronTemperature);
            end
            H5S.close(spaceID);
            H5D.close(dsfID);
          end
          %
          output.saveSwarm(electronKinetics.swarmParam, chemistry.workCond.reducedField, ...
            chemistry.workCond.electronDensity);
        end
      end
      % save rate coefficients info (if selected and acording to the activated modules)
      if output.rateCoeffsIsToBeSaved
        if ~isempty(chemistry.electronKinetics)
          output.saveRateCoefficients(electronKinetics.rateCoeffAll, electronKinetics.rateCoeffExtra, ...
            chemistry.solution.reactionsInfo);
        else
          output.saveRateCoefficients([], [], chemistry.solution.reactionsInfo);
        end
      end
      % save selected results of the chemistry 
      if output.finalDensitiesIsToBeSaved
        output.saveFinalDensities(chemistry.solution.finalDischargeDensity, [chemistry.solution.reactionsInfo.netRate], ...
          chemistry.gasArray, chemistry.electronKinetics);
      end
      if output.finalTemperaturesIsToBeSaved
        output.saveFinalTemperatures(chemistry.workCond.struct);
      end
      if output.finalParticleBalanceIsToBeSaved
        output.saveFinalParticleBalance([chemistry.solution.reactionsInfo.netRate], chemistry.gasArray, ...
          chemistry.reactionArray, chemistry.workCond);
      end
      if output.finalThermalBalanceIsToBeSaved
        output.saveFinalThermalBalance(chemistry.solution.thermalModel);
      end
      if output.chemParamsIsToBeSaved
          output.saveChemParams(chemistry.targetGasPressure, chemistry.workCond.struct, chemistry.convergenceParameter, ...
              chemistry.electronKinetics, chemistry.neutralityRelErrorCurrent);
      end

      if output.chemSolutionTimeIsToBeSaved
        if chemistry.isPulsed && chemistry.workCond.currentTime == chemistry.solution.time(end)
            output.subFolder = [];
            output.saveChemSolutionTime(chemistry.solution.time, chemistry.solution.gasTemperatureTime, ...
              chemistry.solution.nearWallTemperatureTime, chemistry.solution.wallTemperatureTime, ...
              chemistry.solution.densitiesTime, chemistry.gasArray);
        elseif ~chemistry.isPulsed
          output.saveChemSolutionTime(chemistry.solution.time, chemistry.solution.gasTemperatureTime, ...
            chemistry.solution.nearWallTemperatureTime, chemistry.solution.wallTemperatureTime, ...
            chemistry.solution.densitiesTime, chemistry.gasArray);
        end
      end
    end
    
    function saveFinalDensities(output, densities, reactionRates, gasArray, electronKinetics)
    % saveFinalDensities saves the densities of all species considered in the chemistry for
    % the final discharge time of the simulation
      
      % evaluate number of gases and species
      numberOfGases = length(gasArray);
      numberOfSpecies = length(densities);
      
      % determine total gas density and relative creation-destruction rates (final time)
      gasDensities = zeros(1,numberOfGases);
      totalVolumeGasDensity = 0;
      totalSurfaceSiteDensity = 0;
      rateBalances = zeros(1,numberOfSpecies);
      for gas = gasArray
        for state = gas.stateArray
          if strcmp(state.type, 'ele') || strcmp(state.type, 'ion')
            gasDensities(gas.ID) = gasDensities(gas.ID) + densities(state.ID);
          end
          if isempty(state.childArray)
            creationRate = 0;
            for reaction = state.reactionsCreation
              for j = 1:length(reaction.productArray)
                if state.ID == reaction.productArray(j).ID
                  creationRate = creationRate + reaction.productStoiCoeff(j)*reactionRates(reaction.ID);
                  break;
                end
              end
            end
            destructionRate = 0;
            for reaction = state.reactionsDestruction
              for j = 1:length(reaction.reactantArray)
                if state.ID == reaction.reactantArray(j).ID
                  destructionRate = destructionRate + reaction.reactantStoiCoeff(j)*reactionRates(reaction.ID);
                  break;
                end
              end
            end
            rateBalances(state.ID) = (creationRate-destructionRate)/creationRate;
          end
        end
        if gas.isVolumeSpecies
          totalVolumeGasDensity = totalVolumeGasDensity + gasDensities(gas.ID);
        else
          totalSurfaceSiteDensity = totalSurfaceSiteDensity + gasDensities(gas.ID);
        end
      end


      % start writing the files or dataset
      if contains(output.dataFormat, 'txt')
        % evaluate length of the 'Species' column
        speciesColumnLength = 13;
        for gas = gasArray
          speciesColumnLength = max(speciesColumnLength, length(gas.name)+12);
          for state = gas.stateArray
            switch state.type
              case 'ele'
                speciesColumnLength = max(speciesColumnLength, length(state.name)+1);
              case 'vib'
                speciesColumnLength = max(speciesColumnLength, length(state.name)+3);
              case 'rot'
                speciesColumnLength = max(speciesColumnLength, length(state.name)+5);
              case 'ion'
                speciesColumnLength = max(speciesColumnLength, length(state.name)+1);
            end
          end
        end

        % evaluate auxiliary strings for the proper formating of the table
        auxStr1 = sprintf(' %%-%ds ',speciesColumnLength-1);
        auxStr2 = sprintf(' | %%-%ds ',speciesColumnLength-3);
        auxStr3 = sprintf(' | | %%-%ds ',speciesColumnLength-5);

        % create file name
        fileName = [output.folder output.subFolder filesep 'chemFinalDensities.txt'];

        % open file
        fileID = fopen(fileName, 'wt');

        % write volume chemistry information (final densities, final populations and final particle balances)
        fprintf(fileID, '*****************************\n');
        fprintf(fileID, '*  Chemistry (Volume phase) *\n');
        fprintf(fileID, '*****************************\n\n');
        fprintf(fileID, 'Species%s Abs.Density(m^-3)%s Population%s Balance\n%s\n', ...
          repmat(' ', 1, speciesColumnLength-7), repmat(' ', 1, 7), repmat(' ', 1, 14), repmat('-', 1, 98));
        for gas = gasArray
          if gas.isSurfaceSpecies
            continue
          end
          fprintf(fileID, '%s[%f%%]\n', gas.name, 100*gasDensities(gas.ID)/totalVolumeGasDensity);
          for eleState = gas.stateArray
            if strcmp(eleState.type, 'ele')
              if isempty(eleState.childArray)
                fprintf(fileID, [auxStr1 '%#.14e     %#.14e     %+#.14e\n'], eleState.name, densities(eleState.ID), ...
                  densities(eleState.ID)/gasDensities(gas.ID), rateBalances(eleState.ID));
              else
                fprintf(fileID, [auxStr1 '%#.14e     %#.14e\n'], eleState.name, densities(eleState.ID), ...
                  densities(eleState.ID)/gasDensities(gas.ID));
              end
              for vibState = eleState.childArray
                if isempty(vibState.childArray)
                  fprintf(fileID, [auxStr2 '| %#.14e   | %#.14e   | %+#.14e\n'], vibState.name, ...
                    densities(vibState.ID), densities(vibState.ID)/densities(eleState.ID), rateBalances(vibState.ID));
                else
                  fprintf(fileID, [auxStr2 '| %#.14e   | %#.14e\n'], vibState.name, densities(vibState.ID), ...
                    densities(vibState.ID)/densities(eleState.ID));
                end
                for rotState = vibState.childArray
                  fprintf(fileID, [auxStr3 '| | %#.14e | | %#.14e | | %+#.14e\n'], rotState.name, ...
                    densities(rotState.ID), densities(rotState.ID)/densities(vibState.ID), rateBalances(rotState.ID));
                end
              end
            end
          end
          for ionState = gas.stateArray
            if strcmp(ionState.type, 'ion')
              fprintf(fileID, [auxStr1 '%#.14e     %#.14e     %+#.14e\n'], ionState.name, densities(ionState.ID), ...
                  densities(ionState.ID)/gasDensities(gas.ID), rateBalances(ionState.ID));
            end
          end
        end
        % write surface chemistry information (final densities, final populations and final particle balances)
        if totalSurfaceSiteDensity
          fprintf(fileID, '\n*****************************\n');
          fprintf(fileID, '* Chemistry (Surface phase) *\n');
          fprintf(fileID, '*****************************\n\n');

          fprintf(fileID, 'Species%s Abs.Density(m-2)%s Population%s Balance\n%s\n', ...
            repmat(' ', 1, speciesColumnLength-7), repmat(' ', 1, 8), repmat(' ', 1, 14), repmat('-', 1, 98));
          for gas = gasArray
            if gas.isVolumeSpecies
              continue
            end
            fprintf(fileID, '%s[%f%%]\n', gas.name, 100*gasDensities(gas.ID)/totalSurfaceSiteDensity);
            for eleState = gas.stateArray
              if strcmp(eleState.type, 'ele')
                if isempty(eleState.childArray)
                  fprintf(fileID, [auxStr1 '%#.14e     %#.14e     %+#.14e\n'], eleState.name, densities(eleState.ID), ...
                    densities(eleState.ID)/gasDensities(gas.ID), rateBalances(eleState.ID));
                else
                  fprintf(fileID, [auxStr1 '%#.14e     %#.14e\n'], eleState.name, densities(eleState.ID), ...
                    densities(eleState.ID)/gasDensities(gas.ID));
                end
                for vibState = eleState.childArray
                  if isempty(vibState.childArray)
                    fprintf(fileID, [auxStr2 '| %#.14e   | %#.14e   | %+#.14e\n'], vibState.name, ...
                      densities(vibState.ID), densities(vibState.ID)/densities(eleState.ID), rateBalances(vibState.ID));
                  else
                    fprintf(fileID, [auxStr2 '| %#.14e   | %#.14e\n'], vibState.name, densities(vibState.ID), ...
                      densities(vibState.ID)/densities(eleState.ID));
                  end
                  for rotState = vibState.childArray
                    fprintf(fileID, [auxStr3 '| | %#.14e | | %#.14e | | %+#.14e\n'], rotState.name, ...
                      densities(rotState.ID), densities(rotState.ID)/densities(vibState.ID), rateBalances(rotState.ID));
                  end
                end
              end
            end
            for ionState = gas.stateArray
              if strcmp(ionState.type, 'ion')
                fprintf(fileID, [auxStr1 '%#.14e     %#.14e     %+#.14e\n'], ionState.name, densities(ionState.ID), ...
                  densities(ionState.ID)/gasDensities(gas.ID), rateBalances(ionState.ID));
              end
            end
          end
        end

        % close file
        fclose(fileID);

        % save electron kinetics populations (in case it is needed)
        if ~isempty(electronKinetics)
          % evaluate length of the 'Species' column
          speciesColumnLength = 13;
          for gas = electronKinetics.gasArray
            speciesColumnLength = max(speciesColumnLength, length(gas.name)+12);
            for state = gas.stateArray
              switch state.type
                case 'ele'
                  speciesColumnLength = max(speciesColumnLength, length(state.name)+1);
                case 'vib'
                  speciesColumnLength = max(speciesColumnLength, length(state.name)+3);
                case 'rot'
                  speciesColumnLength = max(speciesColumnLength, length(state.name)+5);
                case 'ion'
                  speciesColumnLength = max(speciesColumnLength, length(state.name)+1);
              end
            end
          end
          % evaluate auxiliary strings for the proper formating of the table
          auxStr1 = sprintf(' %%-%ds ',speciesColumnLength-1);
          auxStr2 = sprintf(' | %%-%ds ',speciesColumnLength-3);
          auxStr3 = sprintf(' | | %%-%ds ',speciesColumnLength-5);
          % create file name
          fileName = [output.folder filesep output.subFolder filesep 'electronKineticsFinalPopulations.txt'];
          % open file
          fileID = fopen(fileName, 'wt');
          % write header
          fprintf(fileID, 'Species%s Population\n%s\n', repmat(' ', 1, speciesColumnLength-7), repmat('-', 1, 98));
          % write information by gas
          for gas = electronKinetics.gasArray
            fprintf(fileID, '%s[%f%%]\n', gas.name, 100*gas.fraction);
            for eleState = gas.stateArray
              if strcmp(eleState.type, 'ele')
                fprintf(fileID, [auxStr1 '%#.14e\n'], eleState.name, eleState.population);
                for vibState = eleState.childArray
                  fprintf(fileID, [auxStr2 '| %#.14e\n'], vibState.name, vibState.population);
                  for rotState = vibState.childArray
                    fprintf(fileID, [auxStr3 '| | %#.14e\n'], rotState.name, rotState.population);
                  end
                end
              end
            end
            for ionState = gas.stateArray
              if strcmp(ionState.type, 'ion')
                fprintf(fileID, [auxStr1 '%#.14e\n'], ionState.name, ionState.population);
              end
            end
          end
          % close file
          fclose(fileID);
        end
      end
      if contains(output.dataFormat, 'hdf5')
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        % Create the base types
        strType = H5T.copy('H5T_C_S1');
        doubleType = H5T.copy('H5T_NATIVE_DOUBLE');
        H5T.set_size(strType, 'H5T_VARIABLE');
        sz(1) = H5T.get_size(strType);
        sz(2:4) = H5T.get_size(doubleType);
        offset(1) = 0;
        offset(2:4) = cumsum(sz(1:3));
        dcpl = "H5P_DEFAULT";
        gcID = H5G.open(fID,'/chemistry');
        gdID = H5G.create(gcID,'finalDensities',dcpl,dcpl,dcpl);
        if ~isempty(electronKinetics)
          % Updates the electron density (attribute)
          % - Update initial value in '/'
          attrID = H5A.open(fID,'Electron density (m-3)','H5P_DEFAULT');
          H5A.write(attrID,"H5ML_DEFAULT",electronKinetics.workCond.electronDensity);
          H5A.close(attrID);
          % - Set the value also in '/electronKinetics/swarmParameters'
          if output.swarmParamsIsToBeSaved
            sspaceID = H5S.create("H5S_SCALAR");
            acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
            dssID = H5D.open(fID,'/electronKinetics/swarmParameters');
            attrID = H5A.create(dssID,"Electron density (m-3)",doubleType,sspaceID,acpl);
            H5A.write(attrID,"H5ML_DEFAULT",electronKinetics.workCond.electronDensity);
            H5D.close(dssID);
          end 
          %
          geID = H5G.create(gcID,'electronKineticsPopulations',dcpl,dcpl,dcpl);
        end
        % process each gas in a separate dataset
        % volume species
        for gas = gasArray
          if gas.isSurfaceSpecies
            continue
          end
          % gets the max number of levels to be written
          nState = 0;
          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ele') | strcmp(gasState.type, 'ion')
              nState = nState + 1;  % we always include a sum on dens. and pop.s
              if isempty(gasState.childArray)
                continue;
              else
                for vibState = gasState.childArray
                  nState = nState + 1;
                  if isempty(vibState.childArray)
                    continue;
                  else
                    for rotState = vibState.childArray
                      nState = nState + 1;
                    end
                  end
                end
              end
            end
          end
          % create the dataset for gas.name
          dims = [nState 1];
          h5_dims = fliplr(dims);
          ctypeID = H5T.create('H5T_COMPOUND', sum(sz));
          spaceID = H5S.create_simple(2,h5_dims,[]);
          name = ["Species","Abs.Density", "Population", "Balance"];
          H5T.insert(ctypeID,name(1),offset(1),strType);
          H5T.insert(ctypeID,name(2),offset(2),doubleType);
          H5T.insert(ctypeID,name(3),offset(3),doubleType);
          H5T.insert(ctypeID,name(4),offset(4),doubleType);
          dsID = H5D.create(gdID,gas.name,ctypeID,spaceID,dcpl);
          % populates the gas.name dataset
          data = getData(gas.stateArray);
          % now we can write the dataset         
          dspaceID = H5D.get_space(dsID);
          H5D.write(dsID,ctypeID,spaceID,dspaceID,"H5P_DEFAULT",data);
          % attributes: variable and units in each column
          units = [' -  '; 'm^-3'; ' -  '; ' -  '];
          atdims = 4;
          filetype = H5T.copy('H5T_FORTRAN_S1');
          H5T.set_size(filetype, 4);
          memtype = H5T.copy('H5T_C_S1');
          H5T.set_size(memtype, 4);
          space = H5S.create_simple(1,fliplr(atdims), []);
          attr = H5A.create(dsID, 'Units', filetype, space, 'H5P_DEFAULT');
          H5A.write(attr, memtype, units');
          %
          H5A.close(attr);
          H5S.close(space);
          H5T.close(filetype);
          H5T.close(memtype);
          H5D.close(dsID);
        end
        % surface species
        if totalSurfaceSiteDensity
          for gas = gasArray
            if gas.isVolumeSpecies
              continue
            end
            % gets the max number of levels to be written
            nState = 0;
            for gasState = gas.stateArray
              if strcmp(gasState.type, 'ele') | strcmp(gasState.type, 'ion')
                nState = nState + 1;  % we always include a sum on dens. and pop.s
                if isempty(gasState.childArray)
                  continue;
                else
                  for vibState = gasState.childArray
                    nState = nState + 1;
                    if isempty(vibState.childArray)
                      continue;
                    else
                      for rotState = vibState.childArray
                        nState = nState + 1;
                      end
                    end
                  end
                end
              end
            end
            % create the dataset for gas.name
            dims = [nState 1];
            h5_dims = fliplr(dims);
            ctypeID = H5T.create('H5T_COMPOUND', sum(sz));
            spaceID = H5S.create_simple(2,h5_dims,[]);
            name = ["Species","Abs.Density", "Population", "Balance"];
            H5T.insert(ctypeID,name(1),offset(1),strType);
            H5T.insert(ctypeID,name(2),offset(2),doubleType);
            H5T.insert(ctypeID,name(3),offset(3),doubleType);
            H5T.insert(ctypeID,name(4),offset(4),doubleType);
            dsID = H5D.create(gdID,gas.name,ctypeID,spaceID,dcpl);
            % populates the gas.name dataset
            data = getData(gas.stateArray);
            % now we can write the dataset
            dspaceID = H5D.get_space(dsID);
            H5D.write(dsID,ctypeID,spaceID,dspaceID,"H5P_DEFAULT",data);
            % attributes: variable and units in each column
            units = [' -  '; 'm^-2'; ' -  '; ' -  '];
            atdims = 4;
            filetype = H5T.copy('H5T_FORTRAN_S1');
            H5T.set_size(filetype, 4);
            memtype = H5T.copy('H5T_C_S1');
            H5T.set_size(memtype, 4);
            space = H5S.create_simple(1,fliplr(atdims), []);
            attr = H5A.create(dsID, 'Units', filetype, space, 'H5P_DEFAULT');
            H5A.write(attr, memtype, units');
            %
            H5A.close(attr);
            H5S.close(space);
            H5T.close(filetype);
            H5T.close(memtype);
            %
            H5D.close(dsID);
          end
        end
        if ~isempty(electronKinetics)
          sze(1) = H5T.get_size(strType);
          sze(2) = H5T.get_size(doubleType);
          offset(1) = 0;
          offset(2) = sze(1);
          ctypeID = H5T.create('H5T_COMPOUND', sum(sze));
          name = ["Species","Population"];
          H5T.insert(ctypeID,name(1),offset(1),strType);
          H5T.insert(ctypeID,name(2),offset(2),doubleType);
          for gas = electronKinetics.gasArray
            nState = 0;
            for state = gas.stateArray
              % gets the max number of levels to be written
              if strcmp(state.type, 'ele') | strcmp(state.type, 'ion')
                nState = nState + 1;
                for vibState = state.childArray
                  nState = nState + 1;
                  for rotState = vibState.childArray
                    nState = nState + 1;
                  end
                end
              end
            end
              % create the dataset for gas.name
            h5_dims = fliplr([nState 1]);
            spaceID = H5S.create_simple(2,h5_dims,[]);
            dsID = H5D.create(geID,gas.name,ctypeID,spaceID,dcpl);
            % populates the gas.name dataset
            data = getelectronKineticsData(gas.stateArray);
            % write attribute
            attrID = H5A.create(dsID,'Concentration (%)',doubleType,sspaceID,acpl);
            H5A.write(attrID,"H5ML_DEFAULT",100*gas.fraction);
            % now we can write the dataset
            dspaceID = H5D.get_space(dsID);
            H5D.write(dsID,ctypeID,spaceID,dspaceID,"H5P_DEFAULT",data);
            %
            H5D.close(dsID);
          end
        end
        H5F.close(fID);
      end % output format selection

      % --- nested functions ---------------------------------------------------

      function data = getData(states)
        % Sort states of a given gas species by energy
        
        data.name = {};
        data.density = zeros(1,nState);
        data.population = zeros(1,nState);
        data.balance = zeros(1,nState);
        % Sort states by energy, states without indication of energy are the last
        energy = [];
        for i = 1:length(states)
          if states(i).energy >= 0
            energy(end+1) = states(i).energy;
          else
            energy(end+1) = NaN;
          end
        end
        [B, idx] = sort(energy);
        n = 0;
        for state = states(idx)
          if ~strcmp(state.type, 'ion')
            n = n+1;
            data.name{end+1} = state.name;
            data.density(n) = densities(state.ID);
            data.population(n) = densities(state.ID)/gasDensities(gas.ID);
            data.balance(n) = rateBalances(state.ID);
          end
        end
        % now the ions
        for ionState = states(idx)
          if strcmp(ionState.type, 'ion')
            n = n+1;
            % populates data
            data.name{end+1} = ionState.name;
            data.density(n) = densities(ionState.ID);
            data.population(n) = densities(ionState.ID)/gasDensities(gas.ID);
            data.balance(n) = rateBalances(ionState.ID);
          end
        end

      end   % getData

      function data = getelectronKineticsData(states)
        % Sort states of a given gas species by energy

        data.name = {};
        data.population = zeros(1,nState);
        % Sort states by energy
        energy = [];
        for i = 1:length(states)
          if states(i).energy >= 0
            energy(end+1) = states(i).energy;
          else
            energy(end+1) = NaN;
          end
        end
        [B, idx] = sort(energy);
        n = 0;
        for state = states(idx)
          if ~strcmp(state.type, 'ion')
            n = n+1;
            data.name{end+1} = state.name;
            data.population(n) = state.population;
          end
        end
        for ionState = states(idx)
          if strcmp(ionState.type, 'ion')
            n = n+1;
            data.name{end+1} = ionState.name;
            data.population(n) = ionState.population;
          end
        end
      end   % getelectronKineticsData

    end   % saveFinalDensities
    
    function saveFinalTemperatures(output, workCondStruct)

      possibleTemperatures = {'gasTemperature' 'nearWallTemperature' 'wallTemperature' 'extTemperature'};
      temperatureStr = {'Gas temperature' 'Near wall temperature' 'Wall temperature' 'External temperature'};
      temperatureStrHdf5 = {'Gas temperature (K)' 'Near wall temperature (K)' 'Wall temperature (K)' 'External temperature (K)'};

      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder output.subFolder filesep 'finalTemperatures.txt'];
        
        % open file
        fileID = fopen(fileName, 'wt');

        % save information into the file
        maxStrLength = 0;
        for idx = 1:length(possibleTemperatures)
          maxStrLength = max(maxStrLength, length(temperatureStr{idx}));
        end
        formatSpec = ['%' sprintf('%d', maxStrLength) 's = %#.14e (K)\n'];
        for idx = 1:length(possibleTemperatures)
          if ~isempty(workCondStruct.(possibleTemperatures{idx}))
            fprintf(fileID, formatSpec, temperatureStr{idx}, workCondStruct.(possibleTemperatures{idx}));
          end
        end

        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        % write temperatures as attributes to '/'
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        typeID = H5T.copy("H5T_NATIVE_DOUBLE");
        spaceID = H5S.create("H5S_SCALAR");
        acpl = H5P.create("H5P_ATTRIBUTE_CREATE");        
        for idx = 1:length(possibleTemperatures)
          if ~isempty(workCondStruct.(possibleTemperatures{idx}))
            if strcmp(possibleTemperatures{idx}, 'nearWallTemperature')
              attrID = H5A.create(fID,temperatureStrHdf5{idx},typeID,spaceID,acpl);
            else
              attrID = H5A.open(fID,temperatureStrHdf5{idx},'H5P_DEFAULT');
            end
            H5A.write(attrID,"H5ML_DEFAULT",workCondStruct.(possibleTemperatures{idx}));
            H5A.close(attrID);
          end
        end
        H5S.close(spaceID);
        H5F.close(fID);
      end
      
    end

    function saveFinalParticleBalance(output, reactionRates, gasArray, reactionArray, workCond)
      
      if contains(output.dataFormat, 'txt')

        % evaluate maximum length of reaction descriptions
        maxReactionLength = 0;
        for reaction = reactionArray
          maxReactionLength = max(maxReactionLength, length(reaction.description));
        end
        auxStr = sprintf('      %%-%ds ', maxReactionLength);

        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'finalParticleBalance.txt'];

        % open file
        fileID = fopen(fileName, 'wt');

        % save information into the file
        for gas = gasArray
          fprintf(fileID, '%s\n* Particle balance for %s species *\n%s\n\n', ...
            repmat('*', 1, 33+length(gas.name)), gas.name, repmat('*', 1, 33+length(gas.name)));
          if gas.isVolumeSpecies
            rateUnitsStr = 'm^-3s^-1';
            rateRenorm = 1;
          else
            rateUnitsStr = 'm^-2s^-1';
            rateRenorm = workCond.areaOverVolume;
          end

          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ele')
              if isempty(gasState.childArray)
                state = getRates(gasState);
                writeTxt(state);
              else
                for vibState = gasState.childArray
                  if isempty(vibState.childArray)
                      state = getRates(vibState);
                      writeTxt(state);
                  else
                    for rotState = vibState.childArray
                      state = getRates(rotState);
                      writeTxt(state);
                    end
                  end
                end
              end
            end
            fprintf(fileID, '');
          end
          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ion')
              state = getRates(gasState);
              writeTxt(state);
            end
          end
          fprintf(fileID, '\n');
        end

        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        % prepare the data required
        nSpecies = length(gasArray);
        strType = H5T.copy ('H5T_C_S1');
        H5T.set_size (strType, 'H5T_VARIABLE');
        doubleType = H5T.copy('H5T_NATIVE_DOUBLE');
        sz(1) = H5T.get_size(strType);
        sz(2:3) = H5T.get_size(doubleType);
        offset(1) = 0;
        offset(2:3) = cumsum(sz(1:2));
        dcpl = "H5P_DEFAULT";
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        gcID = H5G.open(fID,'/chemistry');
        gbID = H5G.create(gcID,'particleBalance',dcpl,dcpl,dcpl);
        for gas = gasArray
          if gas.isVolumeSpecies
            rateUnitsStr = 'm^-3/s';
            rateRenorm = 1;
          else
            rateUnitsStr = 'm^-2/s';
            rateRenorm = workCond.areaOverVolume;
          end
          % gets the max number of levels and reactions to be written
          nState = 0;
          nMaxReac = 0;
          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ele') | strcmp(gasState.type, 'ion')
              if isempty(gasState.childArray)
                nState = nState + 1;
                state = getRates(gasState);
                nMaxReac = max([nMaxReac length(state.creationRates) ...
                  length(state.destructionRates)]);
              else
                for vibState = gasState.childArray
                  if isempty(vibState.childArray)
                    nState = nState + 1;
                    state = getRates(vibState);
                    nMaxReac = max([nMaxReac length(state.creationRates) ...
                      length(state.destructionRates)]);
                  else
                    for rotState = vibState.childArray
                      nState = nState + 1;
                      state = getRates(rotState);
                      nMaxReac = max([nMaxReac length(state.creationRates) ...
                        length(state.destructionRates)]);
                    end
                  end
                end
              end
            end
          end

          % create the dataset for gas.name
          dims = [nMaxReac 1 2 nState];
          h5_dims = fliplr(dims);
          ctypeID = H5T.create ('H5T_COMPOUND', sum(sz));
          spaceID = H5S.create_simple(4,h5_dims,h5_dims);
          name = ["Reactions","Rate", "Contribution"];
          H5T.insert(ctypeID,name(1),offset(1),strType);
          H5T.insert(ctypeID,name(2),offset(2),doubleType);
          H5T.insert(ctypeID,name(3),offset(3),doubleType);
          dsID = H5D.create(gbID,gas.name,ctypeID,spaceID,dcpl);
          % populates the gas.name dataset
          stateIdx = -1;
          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ele')
              if isempty(gasState.childArray)
                % write state reactions
                stateIdx = stateIdx + 1;
                state = getRates(gasState);
                writeHdf(state);
              else
                for vibState = gasState.childArray
                  if isempty(vibState.childArray)
                      stateIdx = stateIdx + 1;
                      state = getRates(vibState);
                      writeHdf(state);
                  else
                    for rotState = vibState.childArray
                      stateIdx = stateIdx + 1;
                      state = getRates(rotState);
                      writeHdf(state);
                    end
                  end
                end
              end
            end
          end
          for gasState = gas.stateArray
            if strcmp(gasState.type, 'ion')
              stateIdx = stateIdx + 1;
              state = getRates(gasState);
              writeHdf(state);
            end
          end
          % attributes: variable and units in each column
          units = ['  -   '; rateUnitsStr; '  -   '];
          atdims = 3;
          filetype = H5T.copy('H5T_FORTRAN_S1');
          H5T.set_size(filetype, 6);
          memtype = H5T.copy('H5T_C_S1');
          H5T.set_size(memtype, 6);
          space = H5S.create_simple(1,fliplr(atdims), []);
          attr = H5A.create(dsID, 'Units', filetype, space, 'H5P_DEFAULT');
          H5A.write(attr, memtype, units');
          %
          H5A.close(attr);
          H5S.close(space);
          H5T.close(filetype);
          H5T.close(memtype);         
          %
          H5D.close(dsID);
        end
        H5G.close(gbID);
        H5G.close(gcID);
        H5F.close(fID);
      end
      
      % --- nested functions ---------------------------------------------------
      
      function state = getRates(gasState)

        state.name = gasState.name;
        state.rateUnitsStr = rateUnitsStr;
        state.rateRenorm = rateRenorm;
        % creation channels
        reactions = gasState.reactionsCreation;
        nReactions = length(reactions);
        rates = zeros(1,nReactions);
        state.creationDescription = {};
        for i = 1:nReactions
          if strcmp(gasState.type, 'ion')
            stoiCoeff = 1;
          else
            for j = 1:length(reactions(i).productArray)
              if gasState.ID == reactions(i).productArray(j).ID
                stoiCoeff = reactions(i).productStoiCoeff(j);
                break;
              end
            end
          end
          rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
          state.creationDescription{end+1} = reactions(i).description;
        end
        %
        state.creationRates = [rates(1:nReactions)];
        state.totalCreationRate = sum(rates(:));
        % destruction channels
        reactions = gasState.reactionsDestruction;
        nReactions = length(reactions);
        rates = zeros(1,nReactions);
        state.destructionDescription = {};
        for i = 1:nReactions
          if strcmp(gasState.type, 'ion')
            stoiCoeff = 1;
          else
            for j = 1:length(reactions(i).reactantArray)
              if gasState.ID == reactions(i).reactantArray(j).ID
                stoiCoeff = reactions(i).reactantStoiCoeff(j);
                break;
              end
            end
          end
          rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
          state.destructionDescription{end+1} = reactions(i).description;
        end
        %
        state.destructionRates = [rates(1:nReactions)];
        state.totalDestructionRate = sum(rates(:));

      end

      function writeTxt(state)

        fprintf(fileID, '-> Particle balance for %s:\n', state.name);
        % creation channels
        fprintf(fileID, '    * Reactions where %s is created:\n', state.name);
        fprintf(fileID, '%sRate(%s)       Contribution\n', repmat(' ', 1, maxReactionLength+7), state.rateUnitsStr);
        for i = 1:length(state.creationRates)
          fprintf(fileID, [auxStr '%#.14e %#.14e%%\n'], string(state.creationDescription(i)), ...
            state.creationRates(i), ...
            state.creationRates(i)*100/state.totalCreationRate);
        end
        fprintf(fileID, '%sTOTAL %#.14e %#.14e%%\n', repmat(' ', 1, maxReactionLength+1), state.totalCreationRate, 100);
        % destruction channels
        fprintf(fileID, '    * Reactions where %s is destroyed:\n', state.name);
        fprintf(fileID, '%sRate(%s)       Contribution\n', repmat(' ', 1, maxReactionLength+7), state.rateUnitsStr);
        for i = 1:length(state.destructionRates)
          fprintf(fileID, [auxStr '%#.14e %#.14e%%\n'], string(state.destructionDescription(i)), ...
            state.destructionRates(i), ...
            state.destructionRates(i)*100/state.totalDestructionRate);
        end
        fprintf(fileID, '%sTOTAL %#.14e %#.14e%%\n', repmat(' ', 1, maxReactionLength+1), state.totalDestructionRate, 100);
        % evaluate species balance
        fprintf(fileID, '\n    * Relative %s balance (creation-destruction)/creation: %#.14e%%\n\n', ...
          state.name, (state.totalCreationRate-state.totalDestructionRate)*100/state.totalCreationRate);
      end

      function writeHdf(state)  % sz, name, offset
        % -- Create the compound datatype for memory
        memtype = H5T.create('H5T_COMPOUND', sum(sz));
        H5T.insert(memtype,name(1),offset(1),strType);
        H5T.insert(memtype,name(2),offset(2),doubleType);
        H5T.insert(memtype,name(3),offset(3),doubleType);
        %
        % creation reactions
        data.description = state.creationDescription;
        data.rates = state.creationRates;
        data.contribution = state.creationRates*100/state.totalCreationRate;
        % write the data
        start = [0 0 0 stateIdx];
        h5_start = fliplr(start);
%         h5_block = [1 1 1 length(state.creationRates)+1];
        h5_block = [1 1 1 length(state.creationRates)];
        memSpaceID = H5S.create_simple(4,h5_block,[]);
        dspaceID = H5D.get_space(dsID);
        H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",h5_start,[],[],h5_block);
        H5D.write(dsID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);
        %
        % destruction reactions
        data.description = state.destructionDescription;
        data.rates = state.destructionRates;
        data.contribution = state.destructionRates*100/state.totalDestructionRate;
        % write the data
        start = [0 0 1 stateIdx];
        h5_start = fliplr(start);
%         h5_block = [1 1 1 length(state.destructionRates)+1];
        h5_block = [1 1 1 length(state.destructionRates)];
        memSpaceID = H5S.create_simple(4,h5_block,[]);
        dspaceID = H5D.get_space(dsID);
        H5S.select_hyperslab(dspaceID,"H5S_SELECT_SET",h5_start,[],[],h5_block);
        H5D.write(dsID,memtype,memSpaceID,dspaceID,"H5P_DEFAULT",data);      
      end


    end   % saveFinalParticleBalance

    function saveFinalThermalBalance(output, thermalModel)
      
      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder output.subFolder filesep 'finalThermalBalance.txt'];
        
        % open file
        fileID = fopen(fileName, 'wt');
        
        % save information into the file
        fprintf(fileID, '%30s = %#+.14e (eVm^-3s^-1)\n', 'Conduction', thermalModel.conduction);
        fprintf(fileID, '%30s = %#+.14e (eVm^-3s^-1)\n', 'Electron elastic collisions', thermalModel.elasticCollisions);
        fprintf(fileID, '%30s = %#+.14e (eVm^-3s^-1)\n', 'Volume source', thermalModel.volumeSource);
        fprintf(fileID, '%30s = %#+.14e (eVm^-3s^-1)\n', 'Wall source', thermalModel.wallSource);
        
        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        % write temperatures as attributes to 'chemistry' group
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        gcID = H5G.open(fID,'/chemistry');
        typeID = H5T.copy("H5T_NATIVE_DOUBLE");
        spaceID = H5S.create("H5S_SCALAR");
        acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
        attrID = H5A.create(gcID,'Thermal balance: Conduction',typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",thermalModel.conduction);
        attrID = H5A.create(gcID,'Thermal balance: Electron elastic collisions',typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",thermalModel.elasticCollisions);
        attrID = H5A.create(gcID,'Thermal balance: Volume source',typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",thermalModel.volumeSource);
        attrID = H5A.create(gcID,'Thermal balance: Wall source',typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",thermalModel.wallSource);
        %
        H5S.close(spaceID);
        H5G.close(gcID);
        H5F.close(fID);
      end
      
    end

    function saveChemSolutionTime(output, time, gasTemperatureTime, nearWallTemperatureTime, wallTemperatureTime, ...
        densitiesTime, gasArray)
    % saveChemSolutionTime saves the temporal evolution of all the variables solved in the chemistry
      
      % evaluate number of species and time steps 
      numberOfSpecies = length(densitiesTime(1,:));
      numberOfTimeSteps = length(time);

      % evaluate header of the output file - NOTE: only needed if output.dataFormat is 'txt'
      % ... but we need to know the # of columns
%       name = ["Time(s)", "GasTemp(K)"];
      name = ["Time", "GasTemp"];
      tempColumns = 1;
      if ~isempty(nearWallTemperatureTime)
%         name(end+1) = "NearWallTemp(K)";
        name(end+1) = "NearWallTemp";
        tempColumns = tempColumns+1;
      end
      if ~isempty(wallTemperatureTime)
%         name(end+1) = "WallTemp(K)";
        name(end+1) = "WallTemp";
        tempColumns = tempColumns+1;
      end
      columns = numberOfSpecies+tempColumns+1;
        
      stateIDs = [];
      for gas = gasArray
        if gas.isVolumeSpecies
          unitsStr = '(m^-3)';
        else
          unitsStr = '(m^-2)';
        end
        for state = gas.stateArray
          if strcmp(state.type, 'ele')
            for eleState = [state state.siblingArray]
%               name(end+1) = ['[' eleState.name ']' unitsStr];
              name(end+1) = eleState.name;
              stateIDs(end+1) = eleState.ID;
              if ~isempty(eleState.childArray)
                for vibState = eleState.childArray
%                   name(end+1) = ['[' vibState.name ']' unitsStr];
                  name(end+1) = vibState.name;
                  stateIDs(end+1) = vibState.ID;
                  if ~isempty(vibState.childArray)
                    for rotState = vibState.childArray
%                       name(end+1) = ['[' rotState.name ']' unitsStr];
                      name(end+1) = rotState.name;
                      stateIDs(end+1) = rotState.ID;
                    end
                  end
                end
              end
            end
            break;
          end
        end
        for state = gas.stateArray
          if strcmp(state.type, 'ion')
            for ionState = [state state.siblingArray]
%               name(end+1) = ['[' ionState.name ']' unitsStr];
              name(end+1) = ionState.name;
              stateIDs(end+1) = ionState.ID;
            end
            break;
          end
        end
      end
      
      % evaluate values of the table with the temporal evolution of the densities
      for i = columns:-1:tempColumns+2
        values(i:columns:columns*numberOfTimeSteps) = densitiesTime(:,stateIDs(i-tempColumns-1));
      end
      if ~isempty(wallTemperatureTime)
        values(4:columns:columns*numberOfTimeSteps) = wallTemperatureTime;
      end
      if ~isempty(nearWallTemperatureTime)
        values(3:columns:columns*numberOfTimeSteps) = nearWallTemperatureTime;
      end
      values(2:columns:columns*numberOfTimeSteps) = gasTemperatureTime;
      values(1:columns:columns*numberOfTimeSteps) = time;
        
      if contains(output.dataFormat, 'txt')
        
        % create file name
        fileName = [output.folder output.subFolder filesep 'chemSolutionTime.txt'];
        
        % open file
        fileID = fopen(fileName, 'wt');
        
        % save information into the file
        headerStr = '';
        for i = 1:length(name)
          headerStr = sprintf('%s %-20s', headerStr, name(i));
        end
        fprintf(fileID, '%s\n', headerStr);
        formatSpeStr = [repmat('%#.14e ', 1, columns) '\n'];
        fprintf(fileID, formatSpeStr, values);
        
        % close file
        fclose(fileID);
      end
      if contains(output.dataFormat, 'hdf5')
        
        s = struct(name(1),values(1:columns:columns*numberOfTimeSteps));
        for i = 2:columns
          s = setfield(s, 'd'+string(i), values(i:columns:columns*numberOfTimeSteps));
        end
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        typeID = H5T.copy("H5T_NATIVE_DOUBLE");
        dcpl = "H5P_DEFAULT";
        nd = length(name);
        sz(1:nd) = H5T.get_size(typeID);
        offset(1) = 0;
        offset(2:nd) = cumsum(sz(1:nd-1));
        ctypeID =  H5T.create ('H5T_COMPOUND', sum(sz));
        for i = 1:length(sz)
          H5T.insert(ctypeID,name(i),offset(i),typeID);
        end
        dims = [numberOfTimeSteps 1];
        spaceID = H5S.create_simple(2,fliplr(dims),[]);
        dscID = H5D.create(fID,'/chemistry/chemSolutionTime',ctypeID,spaceID,dcpl);       
%         H5D.write(dscID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT',values);
        dspaceID = H5D.get_space(dscID);
        H5D.write(dscID,ctypeID,dspaceID,dspaceID,'H5P_DEFAULT',s);
        % attributes
        units = ['s   '; 'K   '; 'm^-3'];
        atdims = 3;
        filetype = H5T.copy('H5T_FORTRAN_S1');
        H5T.set_size(filetype, 4);
        memtype = H5T.copy('H5T_C_S1');
        H5T.set_size(memtype, 4);
        space = H5S.create_simple(1,fliplr(atdims), []);
        attr = H5A.create(dscID, 'Units', filetype, space, 'H5P_DEFAULT');
        H5A.write(attr, memtype, units');
        %
        H5A.close(attr);
        H5S.close(space);
        H5T.close(filetype);
        H5T.close(memtype);
        H5S.close(spaceID);
        H5D.close(dscID);
        H5F.close(fID);
      end
    end

    function saveChemParams(output, targetGasPressure, workCondStruct, excitationParameter, electronKineticsStruct, ...
            neutralityRelError)
    % Saves the chemistry parameters information of the current simulation
    
      possibleParameters = {'gasPressure' 'gasTemperature' 'gasDensity' 'totalSccmInFlow' 'totalSccmOutFlow'};
      parameterStr = {'Final gas pressure' 'Final gas temperature' 'Final gas density' 'Total input flow' 'Total output flow'};
      unitStr = {'(Pa)' '(K)' '(m-3)' '(sccm)' '(sccm)'};
      parameterStrHdf5 = {'Final gas pressure (Pa)' 'Final gas temperature (K)' 'Final gas density (m-3)' ...
          'Total input flow (sccm)' 'Total output flow (sccm)'};

      if strcmp(excitationParameter, 'electronDensity')
        possibleParameters(end+1) = {'electronDensity'};
        parameterStr(end+1) = {'Target electron density'};
        unitStr(end+1) = {'(m-3)'};
        parameterStrHdf5(end+1) = {'Target electron density'};
        excParStrLength = length('Target electron density');
      elseif strcmp(excitationParameter, 'dischargeCurrent')
        possibleParameters(end+1:end+2) = {'electronDensity' 'dischargeCurrent'};
        parameterStr(end+1:end+2) = {'Final electron density' 'Target discharge current'};
        unitStr(end+1:end+2) = {'(m-3)' '(A)'};
        parameterStrHdf5(end+1:end+2) = {'Final electron density (m-3)' 'Target discharge current (A)'};
        excParStrLength = max(length('Final discharge current'), length('Target discharge current'));
        excParStrLength = max(excParStrLength, length('Final electron density'));
      elseif strcmp(excitationParameter, 'dischargePowerDensity')
        possibleParameters(end+1:end+2) = {'electronDensity' 'dischargePowerDensity'};
        parameterStr(end+1:end+2) = {'Final electron density' 'Target discharge power density'};
        unitStr(end+1:end+2) = {'(m-3)' '(W m-3)'};
        parameterStrHdf5(end+1:end+2) = {'Final electron density (m-3)' 'Target discharge power density (W m-3)'};
        excParStrLength = max(length('Final discharge power density'), length('Target discharge power density'));
        excParStrLength = max(excParStrLength, length('Final electron density'));
      end                

      if contains(output.dataFormat, 'txt')
        % create file name
        fileName = [output.folder filesep output.subFolder filesep 'chemParameters.txt'];

        % open file
        fileID = fopen(fileName, 'wt');

        % check alignment of strings
        maxStrLength = 0;
        maxStrLength = max(maxStrLength, length('Reduced electric field'));
        maxStrLength = max(maxStrLength, length('Target gas pressure'));
        maxStrLength = max(maxStrLength, excParStrLength);
        for idx = 1:length(possibleParameters)
          maxStrLength = max(maxStrLength, length(parameterStr{idx}));
        end

        formatSpec = ['%' sprintf('%d', maxStrLength) 's = %#.14e %s\n'];
        formatSpec2 = ['%' sprintf('%d', maxStrLength) 's = %#.14e %s %s %#.14e\n'];
        fprintf(fileID, formatSpec, 'Reduced electric field', workCondStruct.reducedField, '(Td)');
        if ~isempty(targetGasPressure)
          fprintf(fileID, formatSpec, 'Target gas pressure', targetGasPressure, '(Pa)');
        end          
        for idx = 1:length(possibleParameters)
          if ~isempty(workCondStruct.(possibleParameters{idx}))
            if strcmp(possibleParameters(idx), 'electronDensity')
              fprintf(fileID, formatSpec2, parameterStr{idx}, workCondStruct.(possibleParameters{idx}), unitStr{idx}, ...
                  '; Rel. error =', neutralityRelError);
            else  
              fprintf(fileID, formatSpec, parameterStr{idx}, workCondStruct.(possibleParameters{idx}), unitStr{idx});
            end
          end
        end

        if strcmp(excitationParameter, 'dischargeCurrent')
          targetDischargeCurrent = workCondStruct.dischargeCurrent;
          finalDischargeCurrent = Constant.electronCharge * workCondStruct.electronDensity * ...
              electronKineticsStruct.swarmParam.driftVelocity * pi * workCondStruct.chamberRadius^2; 
          dischargeCurrentRelError = (targetDischargeCurrent - finalDischargeCurrent) / ...
              targetDischargeCurrent;
          fprintf(fileID, formatSpec2, 'Final discharge current', finalDischargeCurrent, '(A)', '; Rel. error =', ...
              dischargeCurrentRelError);
        elseif strcmp(excitationParameter, 'dischargePowerDensity')
          targetDischargePowerDensity = workCondStruct.dischargePowerDensity;
          finalDischargePowerDensity = workCondStruct.electronDensity * workCondStruct.gasDensity * ...
              electronKineticsStruct.power.field * Constant.electronCharge;
          dischargePowerDensityRelError = (targetDischargePowerDensity - finalDischargePowerDensity) / ...
              targetDischargePowerDensity;
          fprintf(fileID, formatSpec2, 'Final discharge power density', finalDischargePowerDensity, '(W m-3)', '; Rel. error =', ...
              dischargePowerDensityRelError);
        end 

        % close file
        fclose(fileID);
      end

      if contains(output.dataFormat, 'hdf5')
        % write chemistry parameters as attributes to '/'
        fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
        gcID = H5G.open(fID,"/chemistry");
        typeID = H5T.copy("H5T_NATIVE_DOUBLE");
        spaceID = H5S.create("H5S_SCALAR");
        acpl = H5P.create("H5P_ATTRIBUTE_CREATE");
        attrID = H5A.create(gcID,"Reduced electric field (Td)",typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",workCondStruct.reducedField);
        if ~isempty(targetGasPressure)
          attrID = H5A.create(gcID,'Target gas pressure (Pa)',typeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",targetGasPressure);
        end  
        for idx = 1:length(possibleParameters)            
          if ~isempty(workCondStruct.(possibleParameters{idx}))
            if strcmp(possibleParameters{idx}, 'gasPressure')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'gasTemperature')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'gasDensity')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'totalSccmInFlow')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'totalSccmOutFlow')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'electronDensity')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'dischargeCurrent')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);
            elseif strcmp(possibleParameters{idx}, 'dischargePowerDensity')
              attrID = H5A.create(gcID,parameterStrHdf5{idx},typeID,spaceID,acpl);              
            end
            H5A.write(attrID,"H5ML_DEFAULT",workCondStruct.(possibleParameters{idx}));
            H5A.close(attrID);
          end
        end
        if strcmp(excitationParameter, 'dischargeCurrent')
          attrID = H5A.create(gcID,'Final discharge current',typeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",finalDischargeCurrent);
          attrID = H5A.create(gcID,'Discharge current rel. error',typeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",dischargeCurrentRelError);
        elseif strcmp(excitationParameter, 'dischargePowerDensity')
          attrID = H5A.create(gcID,'Final discharge power density',typeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",finalDischargePowerDensity);
          attrID = H5A.create(gcID,'Discharge power density rel. error',typeID,spaceID,acpl);
          H5A.write(attrID,"H5ML_DEFAULT",dischargePowerDensityRelError);
        end
        attrID = H5A.create(gcID,'Neutrality rel. error',typeID,spaceID,acpl);
        H5A.write(attrID,"H5ML_DEFAULT",neutralityRelError);
        H5S.close(spaceID);
        H5G.close(gcID);
        H5F.close(fID);
      end

    end
 end 

end
