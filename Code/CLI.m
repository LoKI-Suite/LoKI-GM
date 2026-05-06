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

classdef CLI < handle
  %CLI Class that defines a Command Line Interface
  %   Objects of this class are the CLI of the simulation. The class has methods that displays the status/results of the
  %   simulation as it progresses in the Matlab command line
  
  properties (Access = private)
    
    setup; 
    setupFileInfoStr;
    collisionArray;
    eedfGasArray;
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
      cli.logStr{end+1} = '*      _          _  _____      ____          ____  ____   _  ___            *';
      cli.logStr{end+1} = '*     | |    ___ | |/ /_ _|    | __ )  __   _|___ \| ___| / |/ _ \           *';
      cli.logStr{end+1} = '*     | |   / _ \|   / | |_____|  _ \  \ \ / / __) |___ \ | | | | |          *';
      cli.logStr{end+1} = '*     | |__| (_) | . \ | |_____| |_) |  \ V / / __/ ___) || | |_| |          *';
      cli.logStr{end+1} = '*     |_____\___/|_|\_\___|    |____/    \_/ |_____|____(_)_|\____/          *';
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
      % store handle array for all the gases in the electron kinetics
      cli.eedfGasArray = cli.setup.electronKineticsGasArray;
      % store handle array for all the collisions in order to display their cross sections
      cli.collisionArray = cli.setup.electronKineticsCollisionArray;
      % add listener to status messages of the electron kinetics
      addlistener(cli.setup.electronKinetics, 'genericStatusMessage', @cli.genericStatusMessage);

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
    
  end
  
end
