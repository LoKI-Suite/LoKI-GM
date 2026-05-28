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

classdef WorkingConditions < handle
  %WorkingConditions Class that contains the information about the working conditions of the simulation
  %   This class contains the information about the working conditins of the simulation. It is also the class
  %   responsible for handling changes in the working conditions by broadcasting any change in them to the
  %   corresponding objects (e.g. when the reduced field is updated the boltzmann object is notified, so it can
  %   reevaluate the field operator).
  
  properties
    gasPressure = [];
    gasTemperature = [];
    nearWallTemperature = [];
    wallTemperature = [];
    extTemperature = [];
    gasDensity = [];
    surfaceSiteDensity = [];
    electronDensity = [];
    electronTemperature = [];
    chamberLength = [];
    chamberRadius = [];
    areaOverVolume = [];
    volumeOverArea = [];
    reducedField = [];
    reducedFieldSI = [];
    excitationFrequency = [];
    reducedExcFreqSI = [];
    dischargeCurrent = [];
    dischargePowerDensity = [];
    totalSccmInFlow = [];
    totalSccmOutFlow = [];
    currentTime = [];             % only used for time-dependent calculations 
  end
  
  events
    updatedGasPressure
    updatedGasTemperature
    updatedGasDensity
    updatedElectronDensity
    updatedElectronTemperature
    updatedChamberLength
    updatedReducedField
    updatedReducedExcitationFrequency
    genericStatusMessage
  end
  
  methods (Access = public)
    
    function workCond = WorkingConditions(setup)
      
      for field = fieldnames(setup.info.workingConditions)'
        workCond.(field{1}) = setup.info.workingConditions.(field{1})(1);
      end
      if ~isempty(workCond.gasPressure)
        workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
      end
      if workCond.excitationFrequency == 0
        workCond.reducedExcFreqSI = 0;
      elseif ~isempty(workCond.gasDensity)
        workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
      end
      workCond.reducedFieldSI = workCond.reducedField*1e-21;
      if ~isempty(workCond.chamberLength) && workCond.chamberLength ~=0 && ...
          ~isempty(workCond.chamberRadius) && workCond.chamberRadius ~= 0    % regular cylinder
        workCond.areaOverVolume = 2./workCond.chamberRadius + 2./workCond.chamberLength;
        workCond.volumeOverArea = 1./workCond.areaOverVolume;
      elseif ~isempty(workCond.chamberRadius) && workCond.chamberRadius ~= 0 % infinitely long cylinder
        workCond.areaOverVolume = 2./workCond.chamberRadius;
        workCond.volumeOverArea = 1./workCond.areaOverVolume;
      elseif ~isempty(workCond.chamberLength) && workCond.chamberLength ~=0  % infinitely wide cylinder (slab)
        workCond.areaOverVolume = 2./workCond.chamberLength;
        workCond.volumeOverArea = 1./workCond.areaOverVolume;
      end
      if isfield(setup.info.workingConditions, 'totalSccmOutFlow') && ...
          ~isnumeric(setup.info.workingConditions.totalSccmOutFlow) 
        if strcmp(setup.info.workingConditions.totalSccmOutFlow,'totalSccmInFlow')
          workCond.totalSccmOutFlow = workCond.totalSccmInFlow;
        else
          workCond.totalSccmOutFlow = [];
        end  
      end
      
    end
    
    function update(workCond, propertiesToUpdate, newValues)
      
      
      if ischar(propertiesToUpdate)
        propertiesToUpdate = {propertiesToUpdate};
      end
      
      for idx = 1:length(propertiesToUpdate)
        workCond.(propertiesToUpdate{idx}) = newValues(idx);
        
        switch propertiesToUpdate{idx}
          case 'gasPressure'
            workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedGasPressure');
            notify(workCond, 'updatedGasDensity');
            notify(workCond, 'updatedReducedExcitationFrequency');
            str = sprintf('\\t- Updated gas pressure (%g Pa).\\n', newValues(idx));

          case 'gasTemperature'
            workCond.gasDensity = workCond.gasPressure/(Constant.boltzmann*workCond.gasTemperature);
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedGasTemperature'); 
            notify(workCond, 'updatedGasDensity');
            notify(workCond, 'updatedReducedExcitationFrequency');
            str = sprintf('\\t- Updated gas temperature (%g K).\\n', newValues(idx));

          case 'electronDensity'
            notify(workCond, 'updatedElectronDensity');
            str = sprintf('\\t- Updated electron density (%g m^-3).\\n', newValues(idx));
            
          case 'electronTemperature'
            notify(workCond, 'updatedElectronTemperature');
            str = sprintf('\\t- Updated electron temperature (%g eV).\\n', newValues(idx));
            
          case 'chamberLength'
            notify(workCond, 'updatedChamberLength');
            str = sprintf('\\t- Updated chamber length (%g m).\\n', newValues(idx));
            
          case 'chamberRadius'
            notify(workCond, 'updatedChamberRadius');
            str = sprintf('\\t- Updated chamber radius (%g m).\\n', newValues(idx));
            
          case 'reducedField'
            workCond.reducedFieldSI = workCond.reducedField*1e-21;
            notify(workCond, 'updatedReducedField');
            str = sprintf('\\t- Updated reduced field (%g Td).\\n', newValues(idx));
            
          case 'excitationFrequency'
            workCond.reducedExcFreqSI = workCond.excitationFrequency*2*pi/workCond.gasDensity;
            notify(workCond, 'updatedExcitationFrequency');
            str = sprintf('\\t- Updated excitation frequency (%g s^-1).\\n', newValues(idx));

        end

        % notify change of property for logging purpose 
        notify(workCond, 'genericStatusMessage', StatusEventData(str, 'status'));
      end
      
    end
    
    function workCondStruct = struct(workCond)
    % struct returns an structure with the properties of the object workCond
    
      for field = fields(workCond)'
        workCondStruct.(field{1}) = workCond.(field{1});
      end
      
    end
    
  end
  
end

