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

classdef CLI < handle
  %CLI Class that defines a Command Line Interface
  %   Objects of this class are the CLI of the simulation. The class has methods that displays the status/results of the
  %   simulation as it progresses in the Matlab command line
  
  properties (Access = private)
    
    setup; 
    setupFileInfoStr;
    collisionArray;
    eedfGasArray;
    chemistryGasArray;
    chemistryStateArray;
    chemistryReactionArray;
    isSimulationHF;
    
  end

  properties (Access = public)
    
    logStr = {};
    
  end
  
  methods (Access = public)
    
    function cli = CLI(setup)
      
      % display code banner (with version info)
      cli.logStr{end+1} = '******************************************************************************'; 
      cli.logStr{end+1} = '*     __    _      __    ____           __ __ ____           __  _           *';
      cli.logStr{end+1} = '*    / /   (_)____/ /_  / __ \____     / //_//  _/___  ___  / /_(_)_________ *';
      cli.logStr{end+1} = '*   / /   / / ___/ __ \/ / / / __ \   / ,<   / // __ \/ _ \/ __/ / ___/ ___/ *';
      cli.logStr{end+1} = '*  / /___/ (__  ) /_/ / /_/ / / / /  / /| |_/ // / / /  __/ /_/ / /__(__  )  *';
      cli.logStr{end+1} = '* /_____/_/____/_.___/\____/_/ /_/  /_/ |_/___/_/ /_/\___/\__/_/\___/____/   *';
      cli.logStr{end+1} = '*                                                                            *';
      cli.logStr{end+1} = '*   _          _  _____       ____ __  __         ____   __   _____          *';
      cli.logStr{end+1} = '*  | |    ___ | |/ /_ _|     / ___|  \/  | __   _|___ \ / /_ |___  |         *';
      cli.logStr{end+1} = '*  | |   / _ \|   / | |_____| |  _| |\/| | \ \ / / __) |  _ \   / /          *';
      cli.logStr{end+1} = '*  | |__| (_) | . \ | |_____| |_| | |  | |  \ V / / __/| (_) | / /           *';
      cli.logStr{end+1} = '*  |_____\___/|_|\_\___|     \____|_|  |_|   \_/ |_____|\___(_)_/            *';
      cli.logStr{end+1} = '*                                                                            *';
      cli.logStr{end+1} = '******************************************************************************';
            
      for idx = 1:length(cli.logStr)
          fprintf('%s\n', cli.logStr{idx});
      end
      
      % store handle to setup object to configure cli after setup file is parsed
      cli.setup = setup;

      % add listener to status messages of the setup object
      addlistener(setup, 'genericStatusMessage', @cli.genericStatusMessage);

    end

    function configure(cli)

      % display the setup info in the CLI
      cli.setupFileInfoStr = cli.setup.unparsedInfo;

      % evaluate flag to change the CLI in the case of HF simulations
      cli.isSimulationHF = cli.setup.workCond.reducedExcFreqSI>0;
      
      % adjust CLI to the type of simulation (ElectronKinetics only, Chemistry only or ElectronKinetics+Chemistry)
      if cli.setup.enableChemistry
        % store handle arrays to the objects used in the chemistry
        cli.chemistryGasArray = cli.setup.chemistryGasArray;
        cli.chemistryStateArray = cli.setup.chemistryStateArray;
        cli.chemistryReactionArray = cli.setup.chemistryReactionArray;
        % add listener to status messages of the chemistry
        addlistener(cli.setup.chemistry, 'genericStatusMessage', @cli.genericStatusMessage);
        % add listener to update the CLI when a new iteration of the pressure cycle is found
        addlistener(cli.setup.chemistry, 'newPressureCycleIteration', @cli.newPressureCycleIteration);
        % add listener to update the CLI when a new iteration of the neutrality cycle is found
        addlistener(cli.setup.chemistry, 'newNeutralityCycleIteration', @cli.newNeutralityCycleIteration);
        % add listener to update the CLI when a new iteration of the global cycle is found
        addlistener(cli.setup.chemistry, 'newGlobalCycleIteration', @cli.newGlobalCycleIteration);
        % store electron collisions + add electronKinetics listener (in case it is enabled)
        if cli.setup.enableElectronKinetics
          % store handle arrays to the objects used in the electron kinetics
          cli.eedfGasArray = cli.setup.electronKineticsGasArray;
          cli.collisionArray = cli.setup.electronKineticsCollisionArray;
          % add listener to status messages of the electron kinetics
          addlistener(cli.setup.electronKinetics, 'genericStatusMessage', @cli.genericStatusMessage);
        end
      else
        % store handle array for all the gases in the electron kinetics
        cli.eedfGasArray = cli.setup.electronKineticsGasArray;
        % store handle array for all the collisions in order to display their cross sections
        cli.collisionArray = cli.setup.electronKineticsCollisionArray;
        % add listener to status messages of the electron kinetics
        addlistener(cli.setup.electronKinetics, 'genericStatusMessage', @cli.genericStatusMessage);
      end

      % add listener of the working conditions object
      addlistener(cli.setup.workCond, 'genericStatusMessage', @cli.genericStatusMessage);
      
    end
    
  end
  
  methods (Access = private)

    function genericStatusMessage(cli, ~, statusEventData)
      
      str = statusEventData.message;
      if endsWith(str, '\n')
        strClean = str(1:end-2);
      else
        strClean = str;
      end
      cli.logStr{end+1} = sprintf(strClean);
      fprintf(str);

    end

    function newPressureCycleIteration(cli, chemistry, ~)

      str = sprintf('\t- New pressure cycle iteration (%d): relative error = %e', ...
        chemistry.pressureIterationCurrent, chemistry.pressureRelErrorCurrent);
      fprintf('%s\n', str);
      cli.logStr{end+1} = str; 

    end

    function newNeutralityCycleIteration(cli, chemistry, ~)

      str = sprintf('\t- New neutrality cycle iteration (%d): relative error = %e', ...
        chemistry.neutralityIterationCurrent, chemistry.neutralityRelErrorCurrent);
      fprintf('%s\n', str);
      cli.logStr{end+1} = str;

    end

    function newGlobalCycleIteration(cli, chemistry, ~)

      str = sprintf('\t- New global cycle iteration (%d): relative error = %e', ...
        chemistry.globalIterationCurrent, chemistry.globalRelErrorCurrent);
      fprintf('%s\n', str);
      cli.logStr{end+1} = str;

    end
    
  end
  
end
