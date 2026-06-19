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

    logIsToBeSaved = false;             % boolean to know if the log (as written by the CLI) must be saved
    inputsAreToBeSaved = false;         % boolean to know if the input files must be saved
    currentJobID = 1;                   % cummulative index of the job value used for data writing
    currentJobIndeces = [1];            % indeces of the job value used for data writing
    extraDims = [];                     % Number of tasks for each property, in addition to E/N, with multiple values.
    eedfIsToBeSaved = false;            % boolean to know if the eedf must be saved
    swarmParamsIsToBeSaved = false;     % boolean to know if the swarm parameters info must be saved
    rateCoeffsIsToBeSaved = false;      % boolean to know if the rate coefficients info must be saved
    powerBalanceIsToBeSaved = false;    % boolean to know if the power balance info must be saved
    lookUpTablesAreToBeSaved = false;   % boolean to know if look-up tables with results must be saved
  end
  
  methods (Access = public)
    
    function output = Output(setup)
      
      setupInfo = setup.info;
      output.dataFormat = setupInfo.output.dataFormat;

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
          % in case multiple jobs refer to a parameter different from 'reduced field'
          outputSubFolderBatches = '';
          if ~strcmp(setup.batches(iBatches).property, 'reducedField')
           % set higher-order folder for a single job or a single 'reduced field' value
            if setup.numberOfBatchTypes == 1 || isscalar(setup.info.workingConditions.reducedField)
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
            if ~strcmp(setup.batches(iBatches).property, 'reducedField')
                output.subFolderBatches = outputSubFolderBatches;
            end
        end
      end
      
      if contains(output.dataFormat, 'hdf5')
        % In this case we only need the root 'Output' folder where the hdfFile
        % is written. This file contains the data organized in groups which are 
        % the equivalent to subfoldes. The hdfFile filename is saved in output.folder.
        % The hdfFile identifier is saved in output h5fid.

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
        if isempty(setup.pulseInfo)
          for i = setup.numberOfBatchTypes:-1:1
            if strcmp(setup.batches(i).property, 'reducedField')
              numberOfEoverNJobs = setup.batches(i).jobs;
              break;
            end
          end
        else
          numberOfEoverNJobs = setup.pulseInfo.samplingPoints+1;
        end

        % Common constants for hdf5 calls
        dcpl = "H5P_DEFAULT";
        doubleType = H5T.copy("H5T_NATIVE_DOUBLE");
        intType = H5T.copy("H5T_NATIVE_INT");
        % creates the file with default library properties (overwrite, ...)
        fID = H5F.create(hdfFile, "H5F_ACC_TRUNC", dcpl, dcpl);

        % Saves working conditions as attributes if single valued
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

        % create group for electron kinetics
        if setup.enableElectronKinetics
          geid = H5G.create(fID,"electronKinetics",dcpl,dcpl,dcpl);
        end
      end

      % save what information must be saved
      dataSets = setup.info.output.dataSets;
      if ischar(dataSets)
        dataSets = {dataSets};
      end

      % save the information if the electron kinetics is HF
      if isempty(setup.pulseInfo)
        numberOfEoverNJobs = setup.batches(1).jobs;
        numberOfJobs = setup.numberOfJobs;
      else
         numberOfEoverNJobs = setup.batches(1).jobs*(setup.pulseInfo.samplingPoints+1);
         numberOfJobs = setup.numberOfJobs*(setup.pulseInfo.samplingPoints+1);         
      end
      workCond = setup.info.workingConditions;
      numberOfJobsFreq = length(workCond.excitationFrequency);
      % numberOfJobsModFreq = numberOfEoverNJobs;
      numberOfJobsModFreq = numberOfJobs/numberOfJobsFreq;      
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

      % If we study the electron kinetics, we write the reducedField or
      % electronTemperature values on geid.
      if contains(output.dataFormat, 'hdf5')
        if strcmpi(setup.info.electronKinetics.eedfType, 'boltzmann')
          % First saves the reducedField values if we don't do the chemistry.
          reducedField = setup.info.workingConditions.reducedField;
          % set dims dimensions
          if isempty(setup.pulseInfo)
            dims = [1 numberOfEoverNJobs];
          else    % reducedField(t) -> columns for t and E/N
            dims = [2 numberOfEoverNJobs];
          end
          spaceID = H5S.create_simple(2,fliplr(dims),[]);
          dsID = H5D.create(geid,'reducedField',doubleType,spaceID,dcpl);
          if isempty(setup.pulseInfo)
            H5DS.set_label(dsID,0,'E/N')
          else
            H5DS.set_label(dsID,0,'time')
            H5DS.set_label(dsID,1,'E/N')
          end
                if isempty(setup.pulseInfo)
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
        elseif strcmpi(setup.info.electronKinetics.eedfType, 'prescribedEedf')
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
          if strcmp(setup.batches(i).property, 'reducedField')
            continue    % We already wrote this dataset
          end
          dims = [1 setup.batches(i).jobs];
          spaceID = H5S.create_simple(2,fliplr(dims),[]);
          dsbID = H5D.create(geid,setup.batches(i).property,doubleType,spaceID,dcpl);
%           H5DS.set_label(dsbID,0,'E/N');   % DEBUG <- we need the quantity name!
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
%           H5DS.set_scale(dsID,setup.batches(i).property);
%           H5T.close(stypeID);
        end
      end

      % Sets the common extraDims for the datasets
      for i = 1:length(setup.batches)
        if strcmpi(setup.batches(i).property, 'reducedField')
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
              dims = [length(setup.energyGrid.cell) 1 numberOfEoverNJobs output.extraDims];
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
%               if output.isSimulationHF(output.currentJobID)
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
              dims = [numberOfEoverNJobs 1 output.extraDims];
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
              dims = [numberOfEoverNJobs 1 3 output.extraDims ngas];
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
              dims = [numberOfEoverNJobs 1 3 output.extraDims];
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
              dims = [nReactions 1 numberOfEoverNJobs output.extraDims];
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
              H5T.close(ctypeID);
              % final clean-up
              H5T.close(strType);
              H5T.close(intType);
            end
          case 'lookUpTables'
            % lookUpTables are always written as a txt file
            output.lookUpTablesAreToBeSaved = true;
        end     % switch dataSet{1}
      end     % for dataSet = dataSets

      % closes the hdf5 objects
      if contains(output.dataFormat, 'hdf5')
        H5D.close(dsID);    % We still had the scale dataset open...
        if setup.enableElectronKinetics
          H5G.close(geid);
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
      
        if setup.enableElectronKinetics
          % add listener to status messages of the electron kinetics object
          addlistener(setup.electronKinetics, 'genericStatusMessage', @output.genericStatusMessage);
        % add listener to output results when a new solution for the EEDF is found
        addlistener(setup.electronKinetics, 'obtainedNewEedf', @output.electronKineticsSolution);
      end

    end
    
  end
  
  methods (Access = private)
    
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
      
      % if output format is hdf5, save the time and reducedField values
      if contains(output.dataFormat, 'hdf5')
          data = [electronKinetics.workCond.currentTime, electronKinetics.workCond.reducedField];
          if length(data) == 2
              % write the values
              fID = H5F.open(output.h5file, "H5F_ACC_RDWR", "H5P_DEFAULT");
              dseID = H5D.open(fID,'/electronKinetics/reducedField');
              start = [0 output.currentJobIndeces(1)-1];
%               start = [0 output.currentJobID-1];
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
        if ~isempty(firstAnisotropy)
          sz(1:3) = H5T.get_size(doubleType);
          offset(2:3)=cumsum(sz(1:2));
        else
          sz(1:2) = H5T.get_size(doubleType);
          offset(2:2)=cumsum(sz(1:1));
        end
        memtype = H5T.create ('H5T_COMPOUND', sum(sz));
        data.Energy = energy;
        data.f0 = eedf;
        if ~isempty(firstAnisotropy)
          data.f1 = firstAnisotropy;
          name = ["Energy" "EEDF" "Anisotropy"];
        else
          name = ["Energy" "EEDF"];
        end
        for i = 1:length(sz)
          H5T.insert(memtype,name(i),offset(i),doubleType);
        end
        extraDims = length(output.currentJobIndeces);
        extraStart = output.currentJobIndeces - ones(1,extraDims);
        start = [0 0 extraStart];
        block = [length(energy) 1 ones(1,extraDims)];
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
%         start = [output.currentJobIndeces(1)-1 0];  % Note: location is 0-based, not 1-based!
%         h5_block = [1 1];
        extraStart = output.currentJobIndeces(2:end) - ones(1,length(output.extraDims));
        start = [output.currentJobIndeces(1)-1 0 extraStart];  % Note: location is 0-based, not 1-based!
        h5_block = [1 1 ones(1,length(output.extraDims))];
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
            extraDims = length(output.currentJobIndeces);
            extraStart = output.currentJobIndeces - ones(1,extraDims);
            start = [ratePosition 0 extraStart];
            h5_block = [1 1 ones(1,extraDims)];
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
        % H5T.close(memtype);
        % H5S.close(dspaceID);
        % H5S.close(memSpaceID);
        H5F.close(fID);
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
        start = [output.currentJobIndeces(1)-1 0 0 extraStart];
        block = [1 1 3 ones(1,length(output.extraDims))];
        h5_block = fliplr(block);
        memSpaceID = H5S.create_simple(length(block),h5_block,[]);
%         start = [output.currentJobIndeces(1)-1 0 0];
%         block = [1 1 3];
%         h5_block = fliplr(block);
%         memSpaceID = H5S.create_simple(3,h5_block,[]);
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
          start = [output.currentJobIndeces(1)-1 0 0 extraStart i-1];
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
    
 end 

end
