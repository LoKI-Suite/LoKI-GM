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

classdef GUI < handle
  %GUI Class that defines a Graphical User Interface
  %   Objects of this class are the GUI of the simulation. The class has methods that create and update the
  %   information that appears in the different panels and elements of the GUI
  
  properties (Access = private)
    
    handle;
    cli;
    eedfGasArray;
    collisionArray;
    chemistryGasArray;
    chemistryStateArray;
    chemistryReactionArray;
    solutions = struct.empty;
    refreshFrequency;
    evolvingParameter;
    isSingleEvolvingParameter = false;
    evolvingParameter2;
    evolvingParameterPopUpMenuStr;
    evolvingParameterPopUpMenuStr2;
    evolvingParameterSwarm;    
    isSimulationHF = false;
    
    setupPanel;
    setupTabGroup;
    setupFileTab;
    setupFileInfo;

    resultsGraphsPanel;
    resultsGraphsTabGroup;
    eedfTab;
    eedfPlot;
    eedfHandles;
    eedfLegend;
    eedfPopUpMenu;
    eedfClearButton;
    anisotropyCheckBox;
    eedfPlottedSolutions;
    eedfLogScaleCheckBox;
    chemistryVsTimeTab;
    chemistryVsTimePlot;
    chemistryVsTimeLegend;
    chemistryVsTimeRunPopUpMenu;
    chemistryVsTimeGasPopUpMenu;
    chemistryVsTimeStatePopUpMenu;
    chemistryVsTimeLogScaleCheckBoxX;
    chemistryVsTimeLogScaleCheckBoxY;
    chemistryVsTimeClearButton;
    isGasTemperaturePlotted = false;
    crossSectionTab;
    crossSectionPlot;
    crossSectionLegend;
    crossSectionPopUpMenu;
    crossSectionClearButton;
    inputRateCoeffTab;
    inputRateCoeffList;
    inputRateCoeffInfo1;
    inputRateCoeffInfo2;
    inputRateCoeffInfo3;
    inputRateCoeffInfo4;
    inputRateCoeffInfo5;
    swarmLegend
    redDiffTab;
    redDiffLogScaleCheckBoxX;
    redDiffLogScaleCheckBoxY;
    redDiffPlot;
    redDiffLegend;
    redMobTab;
    redMobLogScaleCheckBoxX;
    redMobLogScaleCheckBoxY;
    redMobPlot;
    redMobLegend;
    redDiffEnergyTab;
    redDiffEnergyLogScaleCheckBoxX;
    redDiffEnergyLogScaleCheckBoxY;
    redDiffEnergyPlot;
    redDiffEnergyLegend;
    redMobEnergyTab;
    redMobEnergyLogScaleCheckBoxX;
    redMobEnergyLogScaleCheckBoxY;
    redMobEnergyPlot;
    redMobEnergyLegend;
    energyTab;
    energyLogScaleCheckBoxX;
    energyLogScaleCheckBoxY;
    energyPlot;
    TeLegend;    
    redTownsendTab;
    redTownsendLogScaleCheckBoxX;
    redTownsendLogScaleCheckBoxY;
    redTownsendPlot;
    redTownsendLegend;    
    redAttachmentTab;
    redAttachmentLogScaleCheckBoxX;
    redAttachmentLogScaleCheckBoxY;
    redAttachmentPlot;
    redAttachmentLegend;
    powerTab;
    powerLogScaleCheckBoxX;
    powerPlot;
    powerLegend;
    powerFieldColor = [0 0 0];
    powerElasticColor = [1 0 0];
    powerCARColor = [0 1 0];
    powerRotColor = [30 110 50]/255;
    powerVibColor = [0 0 1];
    powerEleColor = [200 0 255]/255;
    powerIonColor = [180 50 50]/255;
    powerAttColor = [0 220 220]/255;
    powerGrowthColor = [200 200 50]/255;
    
    resultsTextPanel;
    resultsTextSubPanel;
    resultsTextPopUpMenu;
    resultsTextTabGroup;
    powerBalanceTab;
    powerBalanceInfo;
    swarmParametersTab;
    swarmParametersInfo;
    electronImpactRateCoeffOutputTab;
    electronImpactRateCoeffOutputInfo;
    finalDensitiesTab;
    finalDensitiesInfo;
    finalTemperaturesTab;
    finalTemperaturesInfo;
    chemParametersTab;
    chemParametersInfo;    
    finalBalanceTab;
    finalBalanceInfo;

    statusPanel;
    statusTabGroup;
    logTab;
    logInfo;
    chemistryIterationsTab;
    globalErrorInfo;
    globalErrorData;
    globalErrorPlot;
    neutralityErrorInfo;
    neutralityErrorData;
    neutralityErrorPlot;
    pressureErrorInfo;
    pressureErrorData;
    pressureErrorPlot;
    
  end
  
  methods (Access = public)
    
    function gui = GUI(setup)
      
      % store handle to the CLI
      gui.cli = setup.cli;

      % add listener to status messages of the setup object
      addlistener(setup, 'genericStatusMessage', @gui.genericStatusMessage);

      % add listener of the working conditions object
      addlistener(setup.workCond, 'genericStatusMessage', @gui.genericStatusMessage);

      % store refresh frequency of the GUI
      if isfield(setup.info.gui, 'refreshFrequency')
        gui.refreshFrequency = setup.info.gui.refreshFrequency; 
      else
        gui.refreshFrequency = 1;
      end
      
      % create window of the GUI
      gui.createWindow();
      
      % display the setup info in the GUI
      gui.setupFileInfo.String = setup.unparsedInfo;
      
      % adjust GUI to the type of simulation (ElectronKinetics only, Chemistry only or ElectronKinetics+Chemistry)
      if setup.enableChemistry
        % store handle array to the objects used in the chemistry
        gui.chemistryGasArray = setup.chemistryGasArray;
        gui.chemistryStateArray = setup.chemistryStateArray;
        gui.chemistryReactionArray = setup.chemistryReactionArray;
        % add listener to status messages of the chemistry object
        addlistener(setup.chemistry, 'genericStatusMessage', @gui.genericStatusMessage);
        % add listener to update the GUI when a new iteration of the pressure cycle is found
        addlistener(setup.chemistry, 'newPressureCycleIteration', @gui.newPressureCycleIteration);
        % add listener to update the GUI when a new iteration of the neutrality cycle is found
        addlistener(setup.chemistry, 'newNeutralityCycleIteration', @gui.newNeutralityCycleIteration);
        % add listener to update the GUI when a new iteration of the global cycle is found
        addlistener(setup.chemistry, 'newGlobalCycleIteration', @gui.newGlobalCycleIteration);
        % add listener to update the GUI when a new solution for the chemistry is found
        addlistener(setup.chemistry, 'obtainedNewChemistrySolution', @gui.newChemistrySolution);
        % create chemistry related tabs
        gui.createChemistryVsTimeTab(~isempty(setup.chemistry.surfacePhaseSpeciesIDs));
        gui.createChemistryIterationsTab();
        gui.createFinalDensitiesTab();
        gui.createFinalTemperaturesTab();
        gui.createChemParametersTab();
        gui.createFinalBalanceTab();
        gui.createInputRateCoeffTab();
        % select by default this tab
        gui.statusTabGroup.SelectedTab = gui.chemistryIterationsTab;
        % store electron collisions + create electronKinetics tabs (in case it is enabled)
        if setup.enableElectronKinetics
          % store handle array for all the gases in the electron kinetics
          gui.eedfGasArray = setup.electronKineticsGasArray;
          % store handle array for all the collisions in order to display their cross sections
          gui.collisionArray = setup.electronKineticsCollisionArray;
          % add listener to status messages of the electron kinetics object
          addlistener(setup.electronKinetics, 'genericStatusMessage', @gui.genericStatusMessage);
          % evaluate flag to change the GUI in the case of HF simulations (initialized as false)
          if setup.workCond.reducedExcFreqSI>0
            gui.isSimulationHF = true;
          end
          % create electronKinetics related tabs
          gui.createEedfTab();
          gui.createCrossSectionTab();
          gui.createPowerBalanceTab();
          gui.createSwarmParametersTab();
          gui.createElectronImpactRateCoeffOutputTab();
        end
      else
        % store handle array for all the gases in the electron kinetics
        gui.eedfGasArray = setup.electronKineticsGasArray;
        % store handle array for all the collisions in order to display their cross sections
        gui.collisionArray = setup.electronKineticsCollisionArray;
        % add listener to status messages of the electron kinetics object
        addlistener(setup.electronKinetics, 'genericStatusMessage', @gui.genericStatusMessage);
        % add listener to update the GUI when a new solution for the EEDF is found
        addlistener(setup.electronKinetics, 'obtainedNewEedf', @gui.newEedf);
        % evaluate flag to change the GUI in the case of HF simulations (initialized as false)
        if setup.workCond.reducedExcFreqSI>0
          gui.isSimulationHF = true;
        end
        % create electronKinetics related tabs
        switch class(setup.electronKinetics)
          case 'Boltzmann'
            if setup.electronKinetics.isTimeDependent
              xLabelText = 'Time (s)';
              gui.evolvingParameter = 'currentTime';
              if isscalar(setup.info.workingConditions.gasTemperature)
                 gui.evolvingParameterPopUpMenuStr = 't = %9.3e (s)';
              else
                 gui.evolvingParameter2 = 'gasTemperature';
                 gui.evolvingParameterPopUpMenuStr2 = 't = %9.3e (s); Tg = %9.3e (K)';
                 gui.evolvingParameterSwarm = 'Tg = %9.3e (K)';
              end                  
            else
              xLabelText = 'Reduced Field (Td)';
              gui.evolvingParameter = 'reducedField';
              if isscalar(setup.info.workingConditions.reducedField)
                 gui.isSingleEvolvingParameter = 'true'; 
              end
              if isscalar(setup.info.workingConditions.gasTemperature)
                gui.evolvingParameterPopUpMenuStr = 'E/N = %9.3e (Td)';
              else  
                gui.evolvingParameter2 = 'gasTemperature';
                gui.evolvingParameterPopUpMenuStr2 = 'E/N = %9.3e (Td); Tg = %9.3e (K)';
                gui.evolvingParameterSwarm = 'Tg = %9.3e (K)';
              end
            end
          case 'PrescribedEedf'
            xLabelText = 'Electron Temperature (eV)';
            gui.evolvingParameter = 'electronTemperature';
            if isscalar(setup.info.workingConditions.gasTemperature)
              gui.evolvingParameterPopUpMenuStr = 'Te = %9.3e (eV)';
            else  
              gui.evolvingParameter2 = 'gasTemperature';
              gui.evolvingParameterPopUpMenuStr2 = 'Te = %9.3e (eV); Tg = %9.3e (K)';
              gui.evolvingParameterSwarm = 'Tg = %9.3e (K)';
            end
        end
        gui.createEedfTab();
        gui.createRedDiffTab(xLabelText);
        gui.createRedMobTab(xLabelText);
        gui.createRedDiffEnergyTab(xLabelText);
        gui.createRedMobEnergyTab(xLabelText);
        gui.createEnergyTab(xLabelText);
        if ~gui.isSimulationHF
          gui.createRedTownsendTab(xLabelText);
          gui.createRedAttachmentTab(xLabelText);
        end
        gui.createPowerTab(xLabelText);
        gui.createCrossSectionTab();
        gui.createPowerBalanceTab();
        gui.createSwarmParametersTab();
        gui.createElectronImpactRateCoeffOutputTab();
      end
      
      % display the gui
      drawnow;
      
    end
    
  end
  
  methods (Access = private)
    
    function createWindow(gui)
      
      % create figure for GUI
      screenSize = get(groot,'ScreenSize');
      gui.handle = figure('name', 'LoKI Simulation Tool', 'OuterPosition', [1 45 screenSize(3) screenSize(4)-45], ...
        'MenuBar', 'none', 'NumberTitle', 'off');

      % create results (graph) panel
      gui.resultsGraphsPanel = uipanel('Parent', gui.handle, 'FontSize', 12, 'FontWeight', 'Bold', 'Title', ...
        'Results (graphical)', 'Position', [0.01 0.51 0.48 0.48]);
      gui.resultsGraphsTabGroup = uitabgroup('Parent', gui.resultsGraphsPanel);
      
      % create results (text) panel and pop up menu for selecting results
      gui.resultsTextPanel = uipanel('Parent', gui.handle, 'FontSize', 12, 'FontWeight', 'Bold', 'Title', ...
        'Results (text)', 'Position', [0.01 0.01 0.48 0.48]);
      gui.resultsTextPopUpMenu = uicontrol('Parent', gui.resultsTextPanel, 'Style', 'popupmenu', 'Units', ...
        'normalized', 'Position', [0.01 0.95 0.3 0.05], 'String', {''}, 'Callback', @gui.resultsTextPopUpMenuHandler);
      gui.resultsTextSubPanel = uipanel('Parent', gui.resultsTextPanel, 'Units', 'normalized', 'BorderType', ...
        'none', 'Position', [0.0 0.0 1.0 0.95]);
      gui.resultsTextTabGroup = uitabgroup('Parent', gui.resultsTextSubPanel);
      
      % create status panel
      gui.statusPanel = uipanel('Parent', gui.handle, 'FontSize', 12, 'FontWeight', 'Bold', 'Title', 'Status', ...
        'Position', [0.51 0.01 0.48 0.48]);
      gui.statusTabGroup = uitabgroup('Parent', gui.statusPanel);
      gui.logTab = uitab('Parent', gui.statusTabGroup, 'Title', 'Simulation Log');
      gui.logInfo = uicontrol('Parent', gui.logTab, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', 'FontName', 'Monospaced', ...
        'Fontsize', 10, 'HorizontalAlignment', 'left', 'String', gui.cli.logStr);
      
      % create setup panel
      gui.setupPanel = uipanel('Parent', gui.handle, 'FontSize', 12, 'FontWeight', 'Bold', 'Title', 'Setup', ...
        'Position', [0.51 0.51 0.48 0.48]);
      gui.setupTabGroup = uitabgroup('Parent', gui.setupPanel);
      gui.setupFileTab = uitab('Parent', gui.setupTabGroup, 'Title', 'Setup file');
      gui.setupFileInfo = uicontrol('Parent', gui.setupFileTab, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', 'FontName', 'Monospaced', ...
        'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function createEedfTab(gui)
      
      gui.eedfTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'EEDF');
      gui.eedfPopUpMenu = uicontrol('Parent', gui.eedfTab, 'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [0.1 0.9 0.3 0.05], 'String', {'All'}, 'Callback', @gui.eedfPopUpMenuHandler);
      gui.eedfClearButton = uicontrol('Parent', gui.eedfTab, 'Style', 'pushbutton', 'Units', 'normalized', ...
        'Position', [0.75 0.9 0.15 0.05], 'String', 'Clear Graph', 'Callback', @gui.clearEedfPlot);
      gui.anisotropyCheckBox = uicontrol('Parent', gui.eedfTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.45 0.9 0.3 0.05], 'String', 'Show First Anisotropy (dashed)', 'Value', 1, 'Callback', ...
        @gui.toggleEedfAnisotropy);
      gui.eedfLogScaleCheckBox = uicontrol('Parent', gui.eedfTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.45 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 1, 'Callback', @gui.changeEedfScale);
      gui.eedfPlot = axes('Parent', gui.eedfTab, 'Units', 'normalized', 'OuterPosition', [0 0 1 0.9], ...
        'Box', 'on', 'YScale', 'log'); 
      xlabel('Energy (eV)');
      ylabel('Distribution Function (eV^{-3/2})');
      hold on;
      
    end
        
    function createChemistryVsTimeTab(gui, isSurfaceKineticsEnabled)
      
      gui.chemistryVsTimeTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Chemistry Vs Time');
      gui.chemistryVsTimeRunPopUpMenu = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'popupmenu', 'Units', ...
        'normalized', 'Position', [0.1 0.9 0.1 0.05], 'String', {'All'}, 'Callback', ...
        @gui.chemistryVsTimeRunPopUpMenuHandler);
      gui.chemistryVsTimeGasPopUpMenu = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'popupmenu', 'Units', ...
        'normalized', 'Position', [0.21 0.9 0.1 0.05], 'String', {gui.chemistryGasArray.name 'All'}, 'Value', ...
        length(gui.chemistryGasArray)+1, 'Callback', @gui.chemistryVsTimeGasPopUpMenuHandler);
      gui.chemistryVsTimeStatePopUpMenu = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'popupmenu', 'Units', ...
        'normalized', 'Position', [0.32 0.9 0.1 0.05], 'String', {gui.chemistryStateArray.name 'All'}, 'Value', ...
        length(gui.chemistryStateArray)+1, 'Callback', @gui.chemistryVsTimeStatePopUpMenuHandler);
      gui.chemistryVsTimeLogScaleCheckBoxX = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'checkbox', ...
        'Units', 'normalized', 'Position', [0.45 0.9 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, ....
        'Callback', @gui.changeChemistryVsTimeXScale);
      gui.chemistryVsTimeLogScaleCheckBoxY = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'checkbox', ...
        'Units', 'normalized', 'Position', [0.45 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 1, ...
        'Callback', @gui.changeChemistryVsTimeYScale);
      gui.chemistryVsTimeClearButton = uicontrol('Parent', gui.chemistryVsTimeTab, 'Style', 'pushbutton', 'Units', ...
        'normalized', 'Position', [0.75 0.9 0.15 0.05], 'String', 'Clear Graph', 'Callback', ...
        @gui.clearChemistryVsTimePlot);
      gui.chemistryVsTimePlot = axes('Parent', gui.chemistryVsTimeTab, 'Units', 'normalized', ...
        'OuterPosition', [0 0 1 0.9], 'Box', 'on', 'XScale', 'log', 'YScale', 'log'); 
      xlabel('Time (s)');
      yyaxis(gui.chemistryVsTimePlot, 'left');
      if isSurfaceKineticsEnabled
        ylabel('Volume Density (m^{-3}) / Surface Density (m^{-2})');
      else
        ylabel('Volume Density (m^{-3})');
      end
      yyaxis(gui.chemistryVsTimePlot, 'right');
      ylabel('Temperature (K)');
      hold on;
      
    end
    
    function createRedDiffTab(gui, xLabelText)
      
      gui.redDiffTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Electron Reduced Diffusion');
      gui.redDiffLogScaleCheckBoxX = uicontrol('Parent', gui.redDiffTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', @gui.changeRedDiffXScale);
      gui.redDiffLogScaleCheckBoxY = uicontrol('Parent', gui.redDiffTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 0, 'Callback', @gui.changeRedDiffYScale);
      gui.redDiffPlot = axes('Parent', gui.redDiffTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Reduced diffusion (1/(ms))');
      hold on;
      
    end
      
    function createRedMobTab(gui, xLabelText)
      
      gui.redMobTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Electron Reduced Mobility');
      if gui.isSimulationHF
        uicontrol('Parent', gui.redMobTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.20 0.85 0.15 0.05], ...
          'HorizontalAlignment', 'left', 'ForegroundColor', 'black', 'String', 'DC reduced mobility');
        uicontrol('Parent', gui.redMobTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.36 0.85 0.15 0.05], ...
          'HorizontalAlignment', 'left', 'ForegroundColor', 'red', 'String', 'Re[HF reduced mobility]');
        uicontrol('Parent', gui.redMobTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.52 0.85 0.15 0.05], ...
          'HorizontalAlignment', 'left', 'ForegroundColor', 'blue', 'String', '-Im[HF reduced mobility]');
      end
      gui.redMobLogScaleCheckBoxX = uicontrol('Parent', gui.redMobTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', @gui.changeRedMobXScale);
      gui.redMobLogScaleCheckBoxY = uicontrol('Parent', gui.redMobTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 0, 'Callback', @gui.changeRedMobYScale);
      gui.redMobPlot = axes('Parent', gui.redMobTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Reduced mobility (1/(msV))');
      hold on;
      
    end
    
    function createRedDiffEnergyTab(gui, xLabelText)
      
      gui.redDiffEnergyTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Electron Reduced Energy Diffusion');
      gui.redDiffEnergyLogScaleCheckBoxX = uicontrol('Parent', gui.redDiffEnergyTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedDiffEnergyXScale);
      gui.redDiffEnergyLogScaleCheckBoxY = uicontrol('Parent', gui.redDiffEnergyTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 0, 'Callback', ...
        @gui.changeRedDiffEnergyYScale);
      gui.redDiffEnergyPlot = axes('Parent', gui.redDiffEnergyTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Reduced energy diffusion (eV/(ms))');
      hold on;
      
    end
      
    function createRedMobEnergyTab(gui, xLabelText)
      
      gui.redMobEnergyTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Electron Reduced Energy Mobility');
      gui.redMobEnergyLogScaleCheckBoxX = uicontrol('Parent', gui.redMobEnergyTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedMobEnergyXScale);
      gui.redMobEnergyLogScaleCheckBoxY = uicontrol('Parent', gui.redMobEnergyTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 0, 'Callback', ...
        @gui.changeRedMobEnergyYScale);
      gui.redMobEnergyPlot = axes('Parent', gui.redMobEnergyTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Reduced energy mobility (eV/(msV))');
      hold on;
      
    end
      
    function createEnergyTab(gui, xLabelText)
      
      gui.energyTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Electron Energies');
      uicontrol('Parent', gui.energyTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.15 0.85 0.15 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', 'red', 'String', 'Electron Temperature');
      uicontrol('Parent', gui.energyTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.31 0.85 0.15 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', 'blue', 'String', 'Characteristic Energy');
      gui.energyLogScaleCheckBoxX = uicontrol('Parent', gui.energyTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', @gui.changeEnergyXScale);
      gui.energyLogScaleCheckBoxY = uicontrol('Parent', gui.energyTab, 'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 0, 'Callback', @gui.changeEnergyYScale);
      gui.energyPlot = axes('Parent', gui.energyTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log');
      xlabel(xLabelText);
      ylabel('Energy (eV)');
      hold on;
      
    end
      
    function createRedTownsendTab(gui,xLabelText)
      
      gui.redTownsendTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Townsend Coefficient');
      gui.redTownsendLogScaleCheckBoxX = uicontrol('Parent', gui.redTownsendTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedTownsendXScale);
      gui.redTownsendLogScaleCheckBoxY = uicontrol('Parent', gui.redTownsendTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedTownsendYScale);
      gui.redTownsendPlot = axes('Parent', gui.redTownsendTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log', 'YScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Red. Townsend Coeff. (m^2)');
      hold on;
      
    end
    
    function createRedAttachmentTab(gui, xLabelText)
      
      gui.redAttachmentTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Attachment Coefficient');
      gui.redAttachmentLogScaleCheckBoxX = uicontrol('Parent', gui.redAttachmentTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedAttachmentXScale);
      gui.redAttachmentLogScaleCheckBoxY = uicontrol('Parent', gui.redAttachmentTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.85 0.3 0.05], 'String', 'Y axis logscale', 'Value', 1, 'Callback', ...
        @gui.changeRedAttachmentYScale);
      gui.redAttachmentPlot = axes('Parent', gui.redAttachmentTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log', 'YScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Red. Attach. Coeff. (m^2)');
      hold on;
      
    end
    
    function createPowerTab(gui, xLabelText)
      
      gui.powerTab = uitab('Parent', gui.resultsGraphsTabGroup, 'Title', 'Power Channels');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.15 0.85 0.07 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerFieldColor, 'String', 'Field');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.22 0.85 0.07 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerElasticColor, 'String', 'Elastic');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.29 0.85 0.08 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerCARColor, 'String', 'CC-CAR');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.37 0.85 0.08 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerRotColor, 'String', 'Rotations');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.45 0.85 0.1 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerVibColor, 'String', 'Vibrations');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.55 0.85 0.1 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerEleColor, 'String', 'Electronic');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.65 0.85 0.1 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerIonColor, 'String', 'Ionization');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.75 0.85 0.1 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerAttColor, 'String', 'Attachment');
      uicontrol('Parent', gui.powerTab, 'Style', 'text', 'Units', 'normalized', 'Position', [0.85 0.85 0.1 0.05], ...
        'HorizontalAlignment', 'left', 'ForegroundColor', gui.powerGrowthColor, 'String', 'Growth');
      gui.powerLogScaleCheckBoxX = uicontrol('Parent', gui.powerTab, 'Style', 'checkbox', 'Units', ...
        'normalized', 'Position', [0.75 0.90 0.3 0.05], 'String', 'X axis logscale', 'Value', 1, 'Callback', ...
        @gui.changePowerXScale);
      gui.powerPlot = axes('Parent', gui.powerTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'XScale', 'log'); 
      xlabel(xLabelText);
      ylabel('Normalized power');
      hold on;
      
    end
    
    function createCrossSectionTab(gui)
      
      gui.crossSectionTab = uitab('Parent', gui.setupTabGroup, 'Title', 'Cross Sections');
      values = cell(1,length(gui.collisionArray)+1);
      for idx = 1:length(gui.collisionArray)
        values{idx} = sprintf('%s', gui.collisionArray(idx).description);
      end
      values(end) = {'All'};
      gui.crossSectionPopUpMenu = uicontrol('Parent', gui.crossSectionTab, 'Style', 'popupmenu', 'Units', ...
        'normalized', 'Position', [0.1 0.9 0.6 0.05], 'String', values, 'Callback', @gui.crossSectionPopUpMenuHandler);
      gui.crossSectionClearButton = uicontrol('Parent', gui.crossSectionTab, 'Style', 'pushbutton', 'Units', ...
        'normalized', 'Position', [0.75 0.9 0.15 0.05], 'String', 'Clear Graph', 'Callback', @gui.clearCrossSectionPlot);
      gui.crossSectionPlot = axes('Parent', gui.crossSectionTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 1 0.9], 'Box', 'on', 'YScale', 'log', 'XScale', 'log'); 
      xlabel('Energy (eV)');
      ylabel('CrossSection (m^2)');
      hold on;
      
    end

    function createInputRateCoeffTab(gui)
      
      gui.inputRateCoeffTab = uitab('Parent', gui.setupTabGroup, 'Title', 'Rate Coefficients');
      values = cell(1,length(gui.chemistryReactionArray));
      for idx = 1:length(gui.chemistryReactionArray)
        values{idx} = sprintf('%s', gui.chemistryReactionArray(idx).description);
      end
      gui.inputRateCoeffList = uicontrol('Parent', gui.inputRateCoeffTab,'Style', 'listbox', 'Units', 'normalized', ...
        'Position', [0.01 0.4 0.98 0.58], 'String', values, 'Callback', @gui.inputRateCoeffListHandler);
      uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', 'characters', 'Position', ...
        [2 9 12 2], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Description:');
      uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', 'characters', 'Position', ...
        [2 7 6 2], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Type:');
      uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', 'characters', 'Position', ...
        [2 5 12 2], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Parameters:');
      uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', 'characters', 'Position', ...
        [2 3 10 2], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Enthalpy:');
      uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', 'characters', 'Position', ...
        [2 1 27 2], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Present in electron kinetics:');
      gui.inputRateCoeffInfo1 = uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', ...
        'characters', 'Position', [14 9 50 2], 'HorizontalAlignment', 'left', 'String', 'N/A');
      gui.inputRateCoeffInfo2 = uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', ...
        'characters', 'Position', [8 7 50 2], 'HorizontalAlignment', 'left', 'String', 'N/A');
      gui.inputRateCoeffInfo3 = uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', ...
        'characters', 'Position', [14 5 50 2], 'HorizontalAlignment', 'left', 'String', 'N/A');
      gui.inputRateCoeffInfo4 = uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', ...
        'characters', 'Position', [12 3 50 2], 'HorizontalAlignment', 'left', 'String', 'N/A');
      gui.inputRateCoeffInfo5 = uicontrol('Parent', gui.inputRateCoeffTab, 'Style', 'text', 'Units', ...
        'characters', 'Position', [29 1 50 2], 'HorizontalAlignment', 'left', 'String', 'N/A');
      
    end
    
    function createPowerBalanceTab(gui)
      
      gui.powerBalanceTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Power Balance');
      gui.powerBalanceInfo = uicontrol('Parent', gui.powerBalanceTab, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', 'FontName', 'Monospaced', ...
        'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function createSwarmParametersTab(gui)
      
      gui.swarmParametersTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Swarm Parameters');
      gui.swarmParametersInfo = uicontrol('Parent', gui.swarmParametersTab, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function createElectronImpactRateCoeffOutputTab(gui)
      
      gui.electronImpactRateCoeffOutputTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Rate Coefficients');
      gui.electronImpactRateCoeffOutputInfo = uicontrol('Parent', gui.electronImpactRateCoeffOutputTab, ...
        'Style', 'edit', 'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function createChemistryIterationsTab(gui)
      
      gui.chemistryIterationsTab = uitab('Parent', gui.statusTabGroup, 'Title', 'Chemistry Iterations');
      uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', 'normalized', 'Position', ...
        [0.03 0.9 0.3 0.05], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Global cycle:');
      gui.globalErrorInfo = uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', ...
        'normalized', 'Position', [0.03 0.8 0.3 0.1], 'HorizontalAlignment', 'left', 'String', ...
        {'iteration: N/A' 'error: N/A'});
      gui.globalErrorPlot = axes('Parent', gui.chemistryIterationsTab, 'Units', 'normalized', 'OuterPosition', ...
        [0 0 0.33 0.7], 'Box', 'on', 'YScale', 'log'); 
      xlabel('Global Iteration');
      hold on;
      uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', 'normalized', 'Position', ...
        [0.36 0.9 0.3 0.05], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Quasineutrality cycle:');
      gui.neutralityErrorInfo = uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', ...
        'normalized', 'Position', [0.36 0.7 0.3 0.2], 'HorizontalAlignment', 'left', 'String', ...
        {'iteration: N/A' 'error: N/A' 'exc. param.: N/A'});
      gui.neutralityErrorPlot = axes('Parent', gui.chemistryIterationsTab, 'Units', 'normalized', 'OuterPosition', ...
        [0.33 0 0.33 0.7], 'Box', 'on', 'YScale', 'log'); 
      xlabel('Neutrality Iteration');
      hold on;
      uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', 'normalized', 'Position', ...
        [0.69 0.9 0.3 0.05], 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'String', 'Pressure cycle:');
      gui.pressureErrorInfo = uicontrol('Parent', gui.chemistryIterationsTab, 'Style', 'text', 'Units', ...
        'normalized', 'Position', [0.69 0.8 0.3 0.1], 'HorizontalAlignment', 'left', 'String', ...
        {'iteration: N/A' 'error: N/A'});
      gui.pressureErrorPlot = axes('Parent', gui.chemistryIterationsTab, 'Units', 'normalized', 'OuterPosition', ...
        [0.66 0 0.33 0.7], 'Box', 'on', 'YScale', 'log');  
      xlabel('Pressure Iteration');
      hold on;
      
    end
    
    function createFinalDensitiesTab(gui)
      
      gui.finalDensitiesTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Final densities');
      gui.finalDensitiesInfo = uicontrol('Parent', gui.finalDensitiesTab, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end

    function createFinalTemperaturesTab(gui)
      
      gui.finalTemperaturesTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Final temperatures');
      gui.finalTemperaturesInfo = uicontrol('Parent', gui.finalTemperaturesTab, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function createChemParametersTab(gui)
      
      gui.chemParametersTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Chem Parameters');
      gui.chemParametersInfo = uicontrol('Parent', gui.chemParametersTab, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end

    
    function createFinalBalanceTab(gui)
      
      gui.finalBalanceTab = uitab('Parent', gui.resultsTextTabGroup, 'Title', 'Final balance');
      gui.finalBalanceInfo = uicontrol('Parent', gui.finalBalanceTab, 'Style', 'edit', ...
        'Units', 'normalized', 'Position', [0.01 0.01 0.98 0.98], 'Max', 2, 'Enable', 'inactive', ...
        'FontName', 'Monospaced', 'Fontsize', 10, 'HorizontalAlignment', 'left');
      
    end
    
    function newEedf(gui, electronKinetics, ~)
      
      % evaluate new solution ID
      newSolutionID = length(gui.solutions)+1;
      
      % save solutions for later use on the gui
      gui.solutions(newSolutionID).eedf = electronKinetics.eedf;
      if isa(electronKinetics, 'Boltzmann')
        gui.solutions(newSolutionID).firstAnisotropy = electronKinetics.firstAnisotropy;
      end
      gui.solutions(newSolutionID).energyValues = electronKinetics.energyGrid.cell;
      gui.solutions(newSolutionID).power = electronKinetics.power;
      gui.solutions(newSolutionID).swarmParam = electronKinetics.swarmParam;
      gui.solutions(newSolutionID).rateCoeffAll = electronKinetics.rateCoeffAll;
      gui.solutions(newSolutionID).rateCoeffExtra = electronKinetics.rateCoeffExtra;
      gui.solutions(newSolutionID).workCond = electronKinetics.workCond.struct;
      
      % add new entry to eedfPopUpMenu
      if isempty(gui.evolvingParameter2)
        newString = sprintf(gui.evolvingParameterPopUpMenuStr, electronKinetics.workCond.(gui.evolvingParameter));
      else
        newString = sprintf(gui.evolvingParameterPopUpMenuStr2, electronKinetics.workCond.(gui.evolvingParameter),...
          electronKinetics.workCond.(gui.evolvingParameter2));
      end  
      if length(gui.eedfPopUpMenu.String) == 1
        newString = [newString; gui.eedfPopUpMenu.String];
      else
        newString = [gui.eedfPopUpMenu.String(1:end-1); newString; gui.eedfPopUpMenu.String(end)];        
      end
      set(gui.eedfPopUpMenu, 'String', newString);

      % add new entry to swarmLegend
      if ~isempty(gui.evolvingParameter2)
        newStringSwarm = sprintf(gui.evolvingParameterSwarm, electronKinetics.workCond.(gui.evolvingParameter2)); 
        if ~isempty(gui.swarmLegend)
          newStringSwarm = [gui.swarmLegend(:,:); newStringSwarm];
        end  
      gui.swarmLegend = newStringSwarm;
      end

      % add new entry to resultsPopUpMenu
      set(gui.resultsTextPopUpMenu, 'String', newString(1:end-1));
      
      % update graphs and results panels with the new solution
      gui.addEedfPlot(newSolutionID, 0);
      gui.updateSwarmParamGraphs(gui.evolvingParameter);
      gui.updatePowerGraphs(gui.evolvingParameter);
      gui.updatePowerBalanceInfo(newSolutionID);
      gui.updateSwarmParamInfo(newSolutionID);
      gui.updateRateCoeffInfo(newSolutionID);
      
      % refresh gui
      if mod(newSolutionID, gui.refreshFrequency) == 0
        drawnow;
      end
      
    end
    
    function addEedfPlot(gui, solutionID, includeLegend)
      
      % add legend for the new plot
      gui.eedfLegend{end+1} = gui.eedfPopUpMenu.String{solutionID};

      % add solutionID to vector of plotted solutions (needed for anisotropy toggle interactivity)
      gui.eedfPlottedSolutions(end+1) = solutionID;

      % add new plot
      if isfield(gui.solutions(solutionID), 'firstAnisotropy') && get(gui.anisotropyCheckBox, 'Value')
        gui.eedfHandles(end+1) = plot(gui.eedfPlot, gui.solutions(solutionID).energyValues, ...
          gui.solutions(solutionID).eedf, '-');
        if gui.eedfPlot.ColorOrderIndex == 1
          gui.eedfPlot.ColorOrderIndex = length(gui.eedfPlot.ColorOrder(:,1));
        else
        gui.eedfPlot.ColorOrderIndex = gui.eedfPlot.ColorOrderIndex-1;
        end
        plot(gui.eedfPlot, gui.solutions(solutionID).energyValues, gui.solutions(solutionID).firstAnisotropy, '--');
      else
        gui.eedfHandles(end+1) = plot(gui.eedfPlot, gui.solutions(solutionID).energyValues, ...
          gui.solutions(solutionID).eedf, '-');
      end
      
      % add legend
      if includeLegend
        legend(gui.eedfHandles, gui.eedfLegend);
      end
      
    end
    
    function updateSwarmParamGraphs(gui, evolvingParameter)
      
      numberOfSolutions = length(gui.solutions);
      inputParamValues = zeros(1,numberOfSolutions);
      redDiff = zeros(1,numberOfSolutions);
      redMob = zeros(1,numberOfSolutions);
      redMobHF = zeros(1,numberOfSolutions);
      redDiffEnergy = zeros(1,numberOfSolutions);
      redMobEnergy = zeros(1,numberOfSolutions);
      Te = zeros(1,numberOfSolutions);
      charE = zeros(1,numberOfSolutions);
      redTown = zeros(1,numberOfSolutions);
      redAtt = zeros(1,numberOfSolutions);
      
      idx = numberOfSolutions;
      inputParamValues(idx) = gui.solutions(idx).workCond.(evolvingParameter);
      redDiff(idx) = gui.solutions(idx).swarmParam.redDiffCoeff;
      redMob(idx) = gui.solutions(idx).swarmParam.redMobility;
      if gui.isSimulationHF
        redMobHF(idx) = gui.solutions(idx).swarmParam.redMobilityHF;
      end
      redDiffEnergy(idx) = gui.solutions(idx).swarmParam.redDiffCoeffEnergy;
      redMobEnergy(idx) = gui.solutions(idx).swarmParam.redMobilityEnergy;
      Te(idx) = gui.solutions(idx).swarmParam.Te;
      charE(idx) = gui.solutions(idx).swarmParam.characEnergy;
      redTown(idx) = gui.solutions(idx).swarmParam.redTownsendCoeff;
      redAtt(idx) = gui.solutions(idx).swarmParam.redAttCoeff;

      % Select legends to write in cases of evolvingParameter2 solution
      % (e.g., avoiding repetitions of Tg-legends for different E/N or t values)
      if ~isempty(gui.swarmLegend)
        for idxLegend = 1:idx
          if strcmp(gui.swarmLegend(idxLegend,:),gui.swarmLegend(idx,:))
            idxCurrent = 1;
            for idx2 = 1:idx-1
                idxCurrent = idxCurrent + ~strcmp(gui.swarmLegend(idx2,:),gui.swarmLegend(idx2+1,:));
            end
            if idxLegend == idx
                if gui.isSimulationHF
                    gui.redMobLegend{3*idx-2} = gui.swarmLegend(idx,:);
                    gui.redMobLegend{3*idx-1} = gui.swarmLegend(idx,:);
                    gui.redMobLegend{3*idx} = gui.swarmLegend(idx,:);                 
                else     
                    gui.redMobLegend{idx} = gui.swarmLegend(idx,:);
                    gui.redTownsendLegend{idx} = gui.swarmLegend(idx,:);
                    gui.redAttachmentLegend{idx} = gui.swarmLegend(idx,:);
                end
                gui.redDiffLegend{idx} = gui.swarmLegend(idx,:);
                gui.redDiffEnergyLegend{idx} = gui.swarmLegend(idx,:);
                gui.redMobEnergyLegend{idx} = gui.swarmLegend(idx,:);
                gui.TeLegend{2*idx-1} = gui.swarmLegend(idx,:);
                gui.TeLegend{2*idx} = gui.swarmLegend(idx,:);
            else
                if gui.isSimulationHF
                    gui.redMobLegend{3*idx-2} = "";
                    gui.redMobLegend{3*idx-1} = "";
                    gui.redMobLegend{3*idx} = "";
                else     
                    gui.redMobLegend{idx} = "";
                    gui.redTownsendLegend{idx} = "";
                    gui .redAttachmentLegend{idx} = "";
                end
                gui.redDiffLegend{idx} = "";
                gui.redDiffEnergyLegend{idx} = "";
                gui.redMobEnergyLegend{idx} = "";
                gui.TeLegend{2*idx-1} = "";
                gui.TeLegend{2*idx} = "";
            end
            break
          end
        end   

        % Generate a different marker for each evolvingParameter2 solution 
        % (e.g., for different Tg)
        markers = {'o', '+', '*', '.', 'x', '_', '|', 'square', 'diamond', '^', 'v', '>', '<', 'pentagram', 'hexagram'};
        % Circular selection of markers
        % Example:  getprop(markers, 3) = '*'
        getFirst = @(v)v{1}; 
        getprop = @(options, idxCurrent)getFirst(circshift(options,-idxCurrent+1));
        markerChosen =  getprop(markers,idxCurrent);

      else

        % Consider the same marker for all solutions 
        markerChosen =  'o';

      end

      % Plot the swarm parameters with different makers for each evolvingParameter2 solution
      plot(gui.redDiffPlot, inputParamValues(idx), redDiff(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redDiffplot');

      if gui.isSimulationHF
        plot(gui.redMobPlot, inputParamValues(idx), redMob(idx), 'Color', 'black', 'LineStyle', 'none', ...
            'Marker', markerChosen, 'Tag', 'redMobplot');
        plot(gui.redMobPlot, inputParamValues(idx), real(redMobHF(idx)), 'Color', 'red', 'LineStyle', 'none', ...
            'Marker', markerChosen, 'Tag', 'redMobplot');
        plot(gui.redMobPlot, inputParamValues(idx), -imag(redMobHF(idx)), 'Color', 'blue', 'LineStyle', 'none', ...
            'Marker', markerChosen, 'Tag', 'redMobplot');
      else
        plot(gui.redMobPlot, inputParamValues(idx), redMob(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redMobplot');
        plot(gui.redTownsendPlot, inputParamValues(idx), redTown(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redTownsendplot');
        plot(gui.redAttachmentPlot, inputParamValues(idx), redAtt(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redAttachmentplot');
      end

      plot(gui.redDiffEnergyPlot, inputParamValues(idx), redDiffEnergy(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redDiffEnergyplot');
      plot(gui.redMobEnergyPlot, inputParamValues(idx), redMobEnergy(idx), 'Color', 'black', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'redMobEnergyplot');
      plot(gui.energyPlot, inputParamValues(idx), Te(idx), 'Color', 'red', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'meanEplot');
      plot(gui.energyPlot, inputParamValues(idx), charE(idx), 'Color', 'blue', 'LineStyle', 'none', ...
          'Marker', markerChosen, 'Tag', 'meanEplot');

      % Update legends in cases of evolvingParameter2 solution  
      if ~isempty(gui.swarmLegend)
        legend(gui.redDiffPlot, gui.redDiffLegend);
        legend(gui.redMobPlot, gui.redMobLegend);
        legend(gui.redMobEnergyPlot, gui.redMobEnergyLegend);
        legend(gui.redDiffEnergyPlot, gui.redDiffEnergyLegend);
        legend(gui.energyPlot, gui.TeLegend);
        if ~gui.isSimulationHF
            legend(gui.redTownsendPlot, gui.redTownsendLegend);
            legend(gui.redAttachmentPlot, gui.redAttachmentLegend);
        end    
      end

    end
    
    function changeRedDiffXScale(gui, ~, ~)
    % changeRedDiffXScale is the callback function of the checkbox "redDiffLogScaleCheckBoxX", it sets the x axis of the 
    % reduced diffusion plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redDiffLogScaleCheckBoxX, 'Value')
        set(gui.redDiffPlot, 'XScale', 'log');
      else
        set(gui.redDiffPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedDiffYScale(gui, ~, ~)
    % changeRedDiffYScale is the callback function of the checkbox "redDiffLogScaleCheckBoxY", it sets the y axis of the
    % reduced diffusion plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redDiffLogScaleCheckBoxY, 'Value')
        set(gui.redDiffPlot, 'YScale', 'log');
      else
        set(gui.redDiffPlot, 'YScale', 'linear');
      end

    end

    function changeRedMobXScale(gui, ~, ~)
    % changeRedMobXScale is the callback function of the checkbox "redMobLogScaleCheckBoxX", it sets the x axis of the 
    % reduced mobility plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redMobLogScaleCheckBoxX, 'Value')
        set(gui.redMobPlot, 'XScale', 'log');
      else
        set(gui.redMobPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedMobYScale(gui, ~, ~)
    % changeRedMobYScale is the callback function of the checkbox "redMobLogScaleCheckBoxY", it sets the y axis of the
    % reduced mobility plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redMobLogScaleCheckBoxY, 'Value')
        set(gui.redMobPlot, 'YScale', 'log');
      else
        set(gui.redMobPlot, 'YScale', 'linear');
      end

    end

    function changeRedDiffEnergyXScale(gui, ~, ~)
    % changeRedDiffEnergyXScale is the callback function of the checkbox "redDiffEnergyLogScaleCheckBoxX", it sets the 
    % x axis of the reduced diffusion energy plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redDiffEnergyLogScaleCheckBoxX, 'Value')
        set(gui.redDiffEnergyPlot, 'XScale', 'log');
      else
        set(gui.redDiffEnergyPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedDiffEnergyYScale(gui, ~, ~)
    % changeRedDiffEnergyYScale is the callback function of the checkbox "redDiffEnergyLogScaleCheckBoxY", it sets the 
    % y axis of the reduced diffusion energy plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redDiffEnergyLogScaleCheckBoxY, 'Value')
        set(gui.redDiffEnergyPlot, 'YScale', 'log');
      else
        set(gui.redDiffEnergyPlot, 'YScale', 'linear');
      end

    end

    function changeRedMobEnergyXScale(gui, ~, ~)
    % changeRedMobEnergyXScale is the callback function of the checkbox "redMobEnergyLogScaleCheckBoxX", it sets the 
    % x axis of the reduced mobility energy plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redMobEnergyLogScaleCheckBoxX, 'Value')
        set(gui.redMobEnergyPlot, 'XScale', 'log');
      else
        set(gui.redMobEnergyPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedMobEnergyYScale(gui, ~, ~)
    % changeRedMobEnergyYScale is the callback function of the checkbox "redMobEnergyLogScaleCheckBoxY", it sets the 
    % y axis of the reduced mobility energy plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redMobEnergyLogScaleCheckBoxY, 'Value')
        set(gui.redMobEnergyPlot, 'YScale', 'log');
      else
        set(gui.redMobEnergyPlot, 'YScale', 'linear');
      end

    end

    function changeEnergyXScale(gui, ~, ~)
    % changeEnergyXScale is the callback function of the checkbox "energyLogScaleCheckBoxX", it sets the x axis of the 
    % energy plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.energyLogScaleCheckBoxX, 'Value')
        set(gui.energyPlot, 'XScale', 'log');
      else
        set(gui.energyPlot, 'XScale', 'linear');
      end
      
    end

    function changeEnergyYScale(gui, ~, ~)
    % changeEnergyYScale is the callback function of the checkbox "energyLogScaleCheckBoxY", it sets the y axis of the 
    % energy plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.energyLogScaleCheckBoxY, 'Value')
        set(gui.energyPlot, 'YScale', 'log');
      else
        set(gui.energyPlot, 'YScale', 'linear');
      end
      
    end

    function changeRedTownsendXScale(gui, ~, ~)
    % changeRedTownsendXScale is the callback function of the checkbox "redTownsendLogScaleCheckBoxX", it sets the x 
    % axis of the reduced Townsend plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redTownsendLogScaleCheckBoxX, 'Value')
        set(gui.redTownsendPlot, 'XScale', 'log');
      else
        set(gui.redTownsendPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedTownsendYScale(gui, ~, ~)
    % changeRedTownsendYScale is the callback function of the checkbox "redTownsendLogScaleCheckBoxY", it sets the y 
    % axis of the reduced Townsend plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redTownsendLogScaleCheckBoxY, 'Value')
        set(gui.redTownsendPlot, 'YScale', 'log');
      else
        set(gui.redTownsendPlot, 'YScale', 'linear');
      end

    end

    function changeRedAttachmentXScale(gui, ~, ~)
    % changeRedAttachmentXScale is the callback function of the checkbox "redAttachmentLogScaleCheckBoxX", it sets the x 
    % axis of the reduced Townsend plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.redAttachmentLogScaleCheckBoxX, 'Value')
        set(gui.redAttachmentPlot, 'XScale', 'log');
      else
        set(gui.redAttachmentPlot, 'XScale', 'linear');
      end
      
    end

    function changeRedAttachmentYScale(gui, ~, ~)
    % changeRedAttachmentYScale is the callback function of the checkbox "redAttachmentLogScaleCheckBoxY", it sets the y 
    % axis of the reduced Townsend plot as linear or logscale acording to the value of the checkbox.

      if get(gui.redAttachmentLogScaleCheckBoxY, 'Value')
        set(gui.redAttachmentPlot, 'YScale', 'log');
      else
        set(gui.redAttachmentPlot, 'YScale', 'linear');
      end

    end

    function changePowerXScale(gui, ~, ~)
    % changePowerXScale is the callback function of the checkbox "powerLogScaleCheckBoxX", it sets the x 
    % axis of the power plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.powerLogScaleCheckBoxX, 'Value')
        set(gui.powerPlot, 'XScale', 'log');
      else
        set(gui.powerPlot, 'XScale', 'linear');
      end
      
    end

    function updatePowerGraphs(gui, evolvingParameter)
      
      numberOfSolutions = length(gui.solutions);
      inputParamValues = zeros(1,numberOfSolutions);
      powerField = zeros(1,numberOfSolutions);
      powerElasticGain = zeros(1,numberOfSolutions);
      powerElasticLoss = zeros(1,numberOfSolutions);
      powerCARGain = zeros(1,numberOfSolutions);
      powerCARLoss = zeros(1,numberOfSolutions);
      powerRotationalGain = zeros(1,numberOfSolutions);
      powerRotationalLoss = zeros(1,numberOfSolutions);
      powerVibrationalGain = zeros(1,numberOfSolutions);
      powerVibrationalLoss = zeros(1,numberOfSolutions);
      powerElectronicGain = zeros(1,numberOfSolutions);
      powerElectronicLoss = zeros(1,numberOfSolutions);
      powerIonization = zeros(1,numberOfSolutions);
      powerAttachment = zeros(1,numberOfSolutions);
      powerGrowth = zeros(1,numberOfSolutions);
      powerRef = zeros(1,numberOfSolutions);
      
      for idx = 1:numberOfSolutions
        inputParamValues(idx) = gui.solutions(idx).workCond.(evolvingParameter);
        powerField(idx) = gui.solutions(idx).power.field;
        powerElasticGain(idx) = gui.solutions(idx).power.elasticGain;
        powerElasticLoss(idx) = gui.solutions(idx).power.elasticLoss;
        powerCARGain(idx) = gui.solutions(idx).power.carGain;
        powerCARLoss(idx) = gui.solutions(idx).power.carLoss;
        powerRotationalGain(idx) = gui.solutions(idx).power.rotationalSup;
        powerRotationalLoss(idx) = gui.solutions(idx).power.rotationalIne;
        powerVibrationalGain(idx) = gui.solutions(idx).power.vibrationalSup;
        powerVibrationalLoss(idx) = gui.solutions(idx).power.vibrationalIne;
        powerElectronicGain(idx) = gui.solutions(idx).power.excitationSup;
        powerElectronicLoss(idx) = gui.solutions(idx).power.excitationIne;
        powerIonization(idx) = gui.solutions(idx).power.ionizationIne;
        powerAttachment(idx) = gui.solutions(idx).power.attachmentIne;
        powerGrowth(idx) = gui.solutions(idx).power.eDensGrowth;
        powerRef(idx) = gui.solutions(idx).power.reference;

        % Normalize swarmLegend
        swarmLegendx = strtrim(string(gui.swarmLegend));        
      end 

      % Generate a different lineStyle or marker for each evolvingParameter2 solution 
      % (e.g., for different Tg)
      lineStyles = {'-', '--', ':', '-+', '-*', '-.', '-x', '-_', '-|', ...
          '-square', '-diamond', '-^', '-v', '->', '-<','-pentagram', '-hexagram'};
      markers = {'o', '+', '*', 	'.', 'x', '_', '|', 'square', 'diamond', ...
          '^', 'v', '>', '<', 'pentagram', 'hexagram'};
      % Circular selection of lineStyles or markers
      % Example:  getprop(lineStyles, 3) = ':'
      getFirst = @(v)v{1}; 
      getprop = @(options, id)getFirst(circshift(options,-id +1));

      % Legend groups
      % to enable selecting legends to write in cases of evolvingParameter2 solution
      % (e.g., avoiding repetitions of Tg-legends for different power channels, and E/N or t values)
      legendLabels = unique(swarmLegendx, 'stable');
      nLegends = numel(legendLabels);
      % Particular case for a single Tg-value 
      if nLegends == 0
        nLegends = 1;
      end  
      
      hField = gobjects(nLegends,1);
      
      for k = 1:nLegends
        % indices belonging to this legend group for multiple Tg-values
        if ~isempty(swarmLegendx)
          idxGroup = swarmLegendx == legendLabels(k);
        % indices belonging to this legend group for a single Tg-value
        else
          idxGroup = 1:numberOfSolutions;
        end  
        if ~any(idxGroup)
           continue
        end

        idLineStyle = k;

        % Single E/N calculation - use markers instead of lines in plot
        if gui.isSingleEvolvingParameter
            markerChosen = getprop(markers, idLineStyle);

            % hField (only legend-visible)
            hField(k) = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerField(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerFieldColor, 'LineStyle', 'none', 'Marker', markerChosen);

            % Other handles (hidden from legend)
            hElasticGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElasticGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerElasticColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hElasticLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElasticLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerElasticColor, 'LineStyle', 'none', 'Marker', markerChosen);       
            hCARGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerCARGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerCARColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hCARLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerCARLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerCARColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hRotationalGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerRotationalGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerRotColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hRotationalLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerRotationalLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerRotColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hVibrationalGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerVibrationalGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerVibColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hVibrationalLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerVibrationalLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerVibColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hElectronicGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElectronicGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerEleColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hElectronicLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElectronicLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerEleColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hIonization = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerIonization(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerIonColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hAttachment = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerAttachment(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerAttColor, 'LineStyle', 'none', 'Marker', markerChosen);
            hGrowth = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerGrowth(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerGrowthColor, 'LineStyle', 'none', 'Marker', markerChosen);

        % Multiple E/N calculations - use lines in plot
        else    
            lineStyleChosen = getprop(lineStyles, idLineStyle);

            % hField (only legend-visible)
            hField(k) = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerField(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerFieldColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);

            % Other handles (hidden from legend)
            hElasticGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElasticGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerElasticColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hElasticLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElasticLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerElasticColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);       
            hCARGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerCARGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerCARColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hCARLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerCARLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerCARColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hRotationalGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerRotationalGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerRotColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hRotationalLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerRotationalLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerRotColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hVibrationalGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerVibrationalGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerVibColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hVibrationalLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerVibrationalLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerVibColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hElectronicGain = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElectronicGain(idxGroup)./ powerRef(idxGroup), ...
                'Color', gui.powerEleColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hElectronicLoss = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerElectronicLoss(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerEleColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hIonization = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerIonization(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerIonColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hAttachment = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerAttachment(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerAttColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2);
            hGrowth = plot(gui.powerPlot, ...
                inputParamValues(idxGroup), powerGrowth(idxGroup)./powerRef(idxGroup), ...
                'Color', gui.powerGrowthColor, 'LineStyle', lineStyleChosen, 'LineWidth', 2); 

        end    
        
        hElasticGain.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hElasticLoss.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hCARGain.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hCARLoss.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hRotationalGain.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hRotationalLoss.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hVibrationalGain.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hVibrationalLoss.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hElectronicGain.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hElectronicLoss.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hIonization.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hAttachment.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hGrowth.Annotation.LegendInformation.IconDisplayStyle = 'off';
      end

      % --- Create legend (only valid handles)
      if ~isempty(swarmLegendx)
        valid = isgraphics(hField);
        legend(gui.powerPlot, hField(valid), legendLabels(valid))
      end  

    end
    
    function addCrossSectionPlot(gui, collID)
      
      % add legend for the new plot
      gui.crossSectionLegend{end+1} = gui.collisionArray(collID).description;
      
      % add new plot
      loglog(gui.crossSectionPlot, gui.collisionArray(collID).rawCrossSection(1,:), ...
        gui.collisionArray(collID).rawCrossSection(2,:), '-', 'Tag', 'eedfplot');
      
      % add legend
      legend(gui.crossSectionPlot, gui.crossSectionLegend);
      
    end
    
    function clearEedfPlot(gui, ~, ~)

      % clear plot
      cla(gui.eedfPlot);
      % clear legend
      legend(gui.eedfPlot, 'off');
      gui.eedfLegend = cell.empty;
      gui.eedfHandles = [];
      % clear vector of plotted solutions
      gui.eedfPlottedSolutions = [];
      
      % refresh gui
      drawnow;
      
    end

    function toggleEedfAnisotropy(gui, ~, ~)
    % toogleEedfAnisotropy is the callback function of the checkbox "anisotropyCheckBox", it adds or removes the first anisotropy
    % from the eedf plot according to the value of the checkbox.

      % local save of currently plotted solution(s)
      currentlyPlottedSolutions = gui.eedfPlottedSolutions;

      % check if legend is active
      if length(gui.eedfPlot.Legend)>0
        includeLegend = true;
      else
        includeLegend = false;
      end

      % clear eedf plot
      gui.clearEedfPlot;

      % plot selected solution(s)
      for solutionID = currentlyPlottedSolutions
        gui.addEedfPlot(solutionID, includeLegend);
      end
      
      % refresh gui
      drawnow;
      
    end
    
    function changeEedfScale(gui, ~, ~)
    % changeEeedfScale is the callback function of the checkbox "eedfLogScaleCheckbox", it sets the y axis of the 
    % eedf plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.eedfLogScaleCheckBox, 'Value')
        set(gui.eedfPlot, 'YScale', 'log');
      else
        set(gui.eedfPlot, 'YScale', 'linear');
      end
      
    end

    function changeChemistryVsTimeXScale(gui, ~, ~)
    % changeChemistryVsTimeXScale is the callback function of the checkbox "chemistryVsTimeLogScaleCheckBoxX", it sets 
    % the x axis of the chemistryVsTime plot as linear or logscale acording to the value of the checkbox.
      
      if get(gui.chemistryVsTimeLogScaleCheckBoxX, 'Value')
        set(gui.chemistryVsTimePlot, 'XScale', 'log');
      else
        set(gui.chemistryVsTimePlot, 'XScale', 'linear');
      end
      
    end

    function changeChemistryVsTimeYScale(gui, ~, ~)
    % changeChemistryVsTimeYScale is the callback function of the checkbox "chemistryVsTimeLogScaleCheckBoxY", it sets 
    % the y axis of the chemistryVsTime plot as linear or logscale acording to the value of the checkbox.
      
      yyaxis(gui.chemistryVsTimePlot, 'left');
      if get(gui.chemistryVsTimeLogScaleCheckBoxY, 'Value')
        set(gui.chemistryVsTimePlot, 'YScale', 'log');
      else
        set(gui.chemistryVsTimePlot, 'YScale', 'linear');
      end
      
    end

    function eedfPopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution(s) to plot
      solutionIDArray = gui.eedfPopUpMenu.Value;
      if solutionIDArray == length(gui.solutions)+1
        gui.clearEedfPlot;
        solutionIDArray = 1:length(gui.solutions);
      end
    
      % plot selected solution(s)
      for solutionID = solutionIDArray
        gui.addEedfPlot(solutionID, 1);
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function clearCrossSectionPlot(gui, ~, ~)

      % clear plot
      cla(gui.crossSectionPlot);
      % clear legend
      legend(gui.crossSectionPlot, 'off');
      gui.crossSectionLegend = cell.empty;
      
      % refresh gui
      drawnow;
      
    end
    
    function crossSectionPopUpMenuHandler(gui, ~, ~)
    
      % evaluate ID of the cross section(s) to plot
      collIDArray = gui.crossSectionPopUpMenu.Value;
      if collIDArray == length(gui.collisionArray)+1
        gui.clearCrossSectionPlot;
        collIDArray = 1:length(gui.collisionArray);
      end
    
      % plot selected cross section(s)
      for collID = collIDArray
        gui.addCrossSectionPlot(collID);
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function inputRateCoeffListHandler(gui, ~, ~)

      % evaluate ID of the reaction to display its information
      reactionID = gui.inputRateCoeffList.Value;
      set(gui.inputRateCoeffInfo1, 'String', gui.chemistryReactionArray(reactionID).description);
      set(gui.inputRateCoeffInfo2, 'String', gui.chemistryReactionArray(reactionID).type);
      parametersStr = '';
      for idx = 1:length(gui.chemistryReactionArray(reactionID).rateCoeffParams)
        parameter = gui.chemistryReactionArray(reactionID).rateCoeffParams{idx};
        if isnumeric(parameter)
          parametersStr = sprintf('%s%g, ', parametersStr, parameter);
        elseif ischar(parameter)
          parametersStr = sprintf('%s%s, ', parametersStr, parameter);
        elseif islogical(parameter)
          if parameter
            logicalStr = 'true';
          else
            logicalStr = 'false';
          end
          parametersStr = sprintf('%s%s, ', parametersStr, logicalStr);
        end
      end
      if length(parametersStr) == 2
        parametersStr = 'N/A';
      else
        parametersStr = parametersStr(1:end-2);
      end
      set(gui.inputRateCoeffInfo3, 'String', parametersStr);
      set(gui.inputRateCoeffInfo4, 'String', sprintf('%g eV', gui.chemistryReactionArray(reactionID).enthalpy));
      if isempty(gui.chemistryReactionArray(reactionID).eedfEquivalent)
        eedfEquivalentStr = 'No';
      else
        if gui.chemistryReactionArray(reactionID).eedfEquivalent.isExtra
          eedfEquivalentStr = 'Yes (Extra)';
        else
          eedfEquivalentStr = 'Yes';
        end
      end
      set(gui.inputRateCoeffInfo5, 'String', eedfEquivalentStr);

    end

    function resultsTextPopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution to show
      solutionID = gui.resultsTextPopUpMenu.Value;
    
      % show selected solution(s)
      if isfield(gui.solutions(solutionID), 'densitiesTime')
        gui.updateFinalDensitiesInfo(solutionID);
      end
      if isfield(gui.solutions(solutionID), 'gasTemperatureTime')
        gui.updateFinalTemperaturesInfo(solutionID);
      end
      if isfield(gui.solutions(solutionID), 'reactionRates')
        gui.updateFinalBalanceInfo(solutionID);
      end
      if isfield(gui.solutions(solutionID), 'power')
        gui.updatePowerBalanceInfo(solutionID);
      end
      if isfield(gui.solutions(solutionID), 'swarmParam')
        gui.updateSwarmParamInfo(solutionID);
      end
      if isfield(gui.solutions(solutionID), 'rateCoeffAll')
        gui.updateRateCoeffInfo(solutionID);
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function updatePowerBalanceInfo(gui, solutionID)
      
      % save local copy of the solution
      power = gui.solutions(solutionID).power;
      % create information to display
      gases = fields(power.gases);
      powerStr = cell(1,44+length(gases)*20);
      powerStr{1} = sprintf('                               Field = %#+.3e (eVm^3s^-1)', power.field);
      powerStr{2} = sprintf('           Elastic collisions (gain) = %#+.3e (eVm^3s^-1)', power.elasticGain);
      powerStr{3} = sprintf('           Elastic collisions (loss) = %#+.3e (eVm^3s^-1)', power.elasticLoss);
      powerStr{4} = sprintf('                          CAR (gain) = %#+.3e (eVm^3s^-1)', power.carGain);
      powerStr{5} = sprintf('                          CAR (loss) = %#+.3e (eVm^3s^-1)', power.carLoss);
      powerStr{6} = sprintf('     Excitation inelastic collisions = %#+.3e (eVm^3s^-1)', power.excitationIne);
      powerStr{7} = sprintf('  Excitation superelastic collisions = %#+.3e (eVm^3s^-1)', power.excitationSup);
      powerStr{8} = sprintf('    Vibrational inelastic collisions = %#+.3e (eVm^3s^-1)', power.vibrationalIne);
      powerStr{9} = sprintf(' Vibrational superelastic collisions = %#+.3e (eVm^3s^-1)', power.vibrationalSup);
      powerStr{10} = sprintf('     Rotational inelastic collisions = %#+.3e (eVm^3s^-1)', power.rotationalIne);
      powerStr{11} = sprintf('  Rotational superelastic collisions = %#+.3e (eVm^3s^-1)', power.rotationalSup);
      powerStr{12} = sprintf('               Ionization collisions = %#+.3e (eVm^3s^-1)', power.ionizationIne);
      powerStr{13} = sprintf('               Attachment collisions = %#+.3e (eVm^3s^-1)', power.attachmentIne);
      powerStr{14} = sprintf('             Electron density growth = %#+.3e (eVm^3s^-1) +', power.eDensGrowth);
      powerStr{15} = [' ' repmat('-', 1, 69)];
      powerStr{16} = sprintf('                       Power Balance = %#+.3e (eVm^3s^-1)', power.balance);
      powerStr{17} = sprintf('              Relative Power Balance = % #.3e%%', power.relativeBalance*100);
      powerStr{18} = '';
      powerStr{19} = '';
      powerStr{20} = sprintf('           Elastic collisions (gain) = %#+.3e (eVm^3s^-1)', power.elasticGain);
      powerStr{21} = sprintf('           Elastic collisions (loss) = %#+.3e (eVm^3s^-1) +', power.elasticLoss);
      powerStr{22} = [' ' repmat('-', 1, 69)];
      powerStr{23} = sprintf('           Elastic collisions (net) = %#+.3e (eVm^3s^-1)', power.elasticNet);
      powerStr{24} = '';
      powerStr{25} = sprintf('                          CAR (gain) = %#+.3e (eVm^3s^-1)', power.carGain);
      powerStr{26} = sprintf('                          CAR (loss) = %#+.3e (eVm^3s^-1) +', power.carLoss);
      powerStr{27} = [' ' repmat('-', 1, 69)];
      powerStr{28} = sprintf('                           CAR (net) = %#+.3e (eVm^3s^-1)', power.carNet);
      powerStr{29} = '';
      powerStr{30} = sprintf('     Excitation inelastic collisions = %#+.3e (eVm^3s^-1)', power.excitationIne);
      powerStr{31} = sprintf('  Excitation superelastic collisions = %#+.3e (eVm^3s^-1) +', power.excitationSup);
      powerStr{32} = [' ' repmat('-', 1, 69)];
      powerStr{33} = sprintf('         Excitation collisions (net) = %#+.3e (eVm^3s^-1)', power.excitationNet);
      powerStr{34} = '';
      powerStr{35} = sprintf('    Vibrational inelastic collisions = %#+.3e (eVm^3s^-1)', power.vibrationalIne);
      powerStr{36} = sprintf(' Vibrational superelastic collisions = %#+.3e (eVm^3s^-1) +', power.vibrationalSup);
      powerStr{37} = [' ' repmat('-', 1, 69)];
      powerStr{38} = sprintf('        Vibrational collisions (net) = %#+.3e (eVm^3s^-1)', power.vibrationalNet);
      powerStr{39} = '';
      powerStr{40} = sprintf('     Rotational inelastic collisions = %#+.3e (eVm^3s^-1)', power.rotationalIne);
      powerStr{41} = sprintf('  Rotational superelastic collisions = %#+.3e (eVm^3s^-1) +', power.rotationalSup);
      powerStr{42} = [' ' repmat('-', 1, 69)];
      powerStr{43} = sprintf('         Rotational collisions (net) = %#+.3e (eVm^3s^-1)', power.rotationalNet);
      
      %power balance by gases
      powerByGas = power.gases;
      index = 44;
      for i = 1:length(gases)
        gas = gases{i};
        powerStr{index} = '';
        powerStr{index+1} = [repmat('*', 1, 35) ' ' gas ' ' repmat('*', 1, 37-length(gas))];
        powerStr{index+2} = '';
        powerStr{index+3} = sprintf('     Excitation inelastic collisions = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).excitationIne);
        powerStr{index+4} = sprintf('  Excitation superelastic collisions = %#+.3e (eVm^3s^-1) +', ...
          powerByGas.(gas).excitationSup);
        powerStr{index+5} = [' ' repmat('-', 1, 69)];
        powerStr{index+6} = sprintf('         Excitation collisions (net) = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).excitationNet);
        powerStr{index+7} = '';
        powerStr{index+8} = sprintf('    Vibrational inelastic collisions = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).vibrationalIne);
        powerStr{index+9} = sprintf(' Vibrational superelastic collisions = %#+.3e (eVm^3s^-1) +', ...
          powerByGas.(gas).vibrationalSup);
        powerStr{index+10} = [' ' repmat('-', 1, 69)];
        powerStr{index+11} = sprintf('        Vibrational collisions (net) = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).vibrationalNet);
        powerStr{index+12} = '';
        powerStr{index+13} = sprintf('     Rotational inelastic collisions = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).rotationalIne);
        powerStr{index+14} = sprintf('  Rotational superelastic collisions = %#+.3e (eVm^3s^-1) +', ...
          powerByGas.(gas).rotationalSup);
        powerStr{index+15} = [' ' repmat('-', 1, 69)];
        powerStr{index+16} = sprintf('         Rotational collisions (net) = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).rotationalNet);
        powerStr{index+17} = '';
        powerStr{index+18} = sprintf('               Ionization collisions = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).ionizationIne);
        powerStr{index+19} = sprintf('               Attachment collisions = %#+.3e (eVm^3s^-1)', ...
          powerByGas.(gas).attachmentIne);
        index = index+20;
      end

      %update the powerBalanceInfo object (uicontrol object)
      set(gui.powerBalanceInfo, 'String', powerStr);
      
    end
    
    function updateSwarmParamInfo(gui, solutionID)
    
      % save local copy of the solution
      electronDensity = gui.solutions(solutionID).workCond.electronDensity;
      swarmParam = gui.solutions(solutionID).swarmParam;
      reducedField = gui.solutions(solutionID).workCond.reducedField;
      % create information to display
      swarmStr = cell(0);
      swarmStr{end+1} = sprintf('                     Electron density = %#.3e (m^-3)', electronDensity);      swarmStr{end+1} = sprintf('               Reduced electric field = %#.3e (Td)', reducedField);
      swarmStr{end+1} = sprintf('        Reduced diffusion coefficient = %#.3e ((ms)^-1)', swarmParam.redDiffCoeff);
      swarmStr{end+1} = sprintf('                     Reduced mobility = %#.3e ((msV)^-1)', swarmParam.redMobility);
      if gui.isSimulationHF
        swarmStr{end+1} = sprintf('                  Reduced mobility HF = %#.3e%+#.3ei ((msV)^-1)', ...
          real(swarmParam.redMobilityHF), imag(swarmParam.redMobilityHF));
      else
        swarmStr{end+1} = sprintf('                       Drift velocity = %#.3e (ms^-1)', swarmParam.driftVelocity);
        swarmStr{end+1} = sprintf('         Reduced Townsend coefficient = %#.3e (m^2)', swarmParam.redTownsendCoeff);
        swarmStr{end+1} = sprintf('       Reduced attachment coefficient = %#.3e (m^2)', swarmParam.redAttCoeff);
      end
      swarmStr{end+1} = sprintf(' Reduced energy diffusion coefficient = %#.3e (eV(ms)^-1)', swarmParam.redDiffCoeffEnergy);
      swarmStr{end+1} = sprintf('              Reduced energy mobility = %#.3e (eV(msV)^-1)', swarmParam.redMobilityEnergy);
      swarmStr{end+1} = sprintf('                          Mean energy = %#.3e (eV)', swarmParam.meanEnergy);
      swarmStr{end+1} = sprintf('                Characteristic energy = %#.3e (eV)', swarmParam.characEnergy);
      swarmStr{end+1} = sprintf('                 Electron temperature = %#.3e (eV)', swarmParam.Te);
      
      %update the transportParametersInfo object (uicontrol object)
      set(gui.swarmParametersInfo, 'String', swarmStr);
    
    end
    
    function updateRateCoeffInfo(gui, solutionID)
    
      % save local copy of the solution
      rateCoeffAll = gui.solutions(solutionID).rateCoeffAll;
      rateCoeffExtra = gui.solutions(solutionID).rateCoeffExtra;
      % evaluate number of collisions (regular + extra)
      numberOfCollisions = length(rateCoeffAll);
      numberOfExtraCollisions = length(rateCoeffExtra);
      
      % create information to display
      numberOfCells = numberOfCollisions+numberOfExtraCollisions;
      if ~isempty(rateCoeffAll)
        numberOfCells = numberOfCells + 9;
      end
      if ~isempty(rateCoeffExtra)
        numberOfCells = numberOfCells + 9;
      end
      rateCoeffStr = cell(1,numberOfCells);
      if ~isempty(rateCoeffAll)
        rateCoeffStr{1} = repmat('*', 1,38);
        rateCoeffStr{2} = '*    e-Kinetics Rate Coefficients    *';
        rateCoeffStr{3} = repmat('*', 1,38);
        rateCoeffStr{4} = '';
        rateCoeffStr{5} = 'ID   Inel.     Superel.  Description';
        rateCoeffStr{6} = '     (m^3s^-1) (m^3s^-1)';
        rateCoeffStr{7} = repmat('-', 1,80);
        for idx = 1:numberOfCollisions
          if length(rateCoeffAll(idx).value) == 1
            rateCoeffStr{7+idx} = sprintf('%4d %9.3e (N/A)     %s', rateCoeffAll(idx).collID, ...
              rateCoeffAll(idx).value, rateCoeffAll(idx).collDescription);
          else
            rateCoeffStr{7+idx} = sprintf('%4d %9.3e %9.3e %s', rateCoeffAll(idx).collID, ...
              rateCoeffAll(idx).value(1), rateCoeffAll(idx).value(2), rateCoeffAll(idx).collDescription);
          end
        end
        rateCoeffStr{8+numberOfCollisions} = repmat('-', 1,80);
        rateCoeffStr{9+numberOfCollisions} = '';
      end
      if ~isempty(rateCoeffExtra)
        if isempty(rateCoeffAll)
          initialIdx = 0;
        else
          initialIdx = 9+numberOfCollisions;
        end
        rateCoeffStr{initialIdx+1} = repmat('*', 1,38);
        rateCoeffStr{initialIdx+2} = '* e-Kinetics Extra Rate Coefficients *';
        rateCoeffStr{initialIdx+3} = repmat('*', 1,38);
        rateCoeffStr{initialIdx+4} = '';
        rateCoeffStr{initialIdx+5} = 'ID   Inel.     Superel.  Description';
        rateCoeffStr{initialIdx+6} = '     (m^3s^-1) (m^3s^-1)';
        rateCoeffStr{initialIdx+7} = repmat('-', 1,80);
        for idx = 1:numberOfExtraCollisions
          if length(rateCoeffExtra(idx).value) == 1
            rateCoeffStr{initialIdx+7+idx} = sprintf('%4d %9.3e (N/A)     %s', rateCoeffExtra(idx).collID, ...
              rateCoeffExtra(idx).value, rateCoeffExtra(idx).collDescription);
          else
            rateCoeffStr{initialIdx+7+idx} = sprintf('%4d %9.3e %9.3e %s', rateCoeffExtra(idx).collID, ...
              rateCoeffExtra(idx).value(1), rateCoeffExtra(idx).value(2), rateCoeffExtra(idx).collDescription);
          end
        end
        rateCoeffStr{initialIdx+8+numberOfExtraCollisions} = repmat('-', 1,80);
        rateCoeffStr{initialIdx+9+numberOfExtraCollisions} = '';
      end
      
      %update the powerBalanceInfo object (uicontrol object)
      set(gui.electronImpactRateCoeffOutputInfo, 'String', rateCoeffStr);
    
    end
    
    function newPressureCycleIteration(gui, chemistry, ~)
      
      % add data of the new iteration to the array pressureErrorData
      if isempty(gui.pressureErrorData)
        data = [chemistry.pressureIterationCurrent chemistry.pressureRelErrorCurrent];
      else
        data = [gui.pressureErrorData; chemistry.pressureIterationCurrent chemistry.pressureRelErrorCurrent];
      end
      gui.pressureErrorData = data;
      
      % Separate positive and negative errors
      pos_indices = data(:,2) >= 0;
      neg_indices = data(:,2) < 0;
      
      % Plot positive errors in blue
      h1 = plot(gui.pressureErrorPlot, data(pos_indices,1), abs(data(pos_indices,2)), 'o', ...
        'MarkerEdgeColor', 'blue', 'DisplayName', 'Positive Errors');
      hold(gui.pressureErrorPlot, 'on');
      
      % Plot negative errors in red
      h2 = plot(gui.pressureErrorPlot, data(neg_indices,1), abs(data(neg_indices,2)), 'o', ...
        'MarkerEdgeColor', 'red', 'DisplayName', 'Negative Errors');
      
      % Plot all errors with a connecting line in black
      plot(gui.pressureErrorPlot, data(:,1), abs(data(:,2)), '-', 'Color', 'black');
      
      % Add legend
      legend(gui.pressureErrorPlot, [h1, h2], 'Location', 'northeast');

      % Update text in the chemistry iterations tab
      pressureErrorInfoStr = cell(1,2);
      pressureErrorInfoStr{1} = sprintf('iteration: %d', data(end,1));
      pressureErrorInfoStr{2} = sprintf('error: %g', data(end,2));
      set(gui.pressureErrorInfo, 'String', pressureErrorInfoStr);

      % add info to log tab
      gui.logInfo.String{end+1} = sprintf('\t- New pressure cycle iteration (%d): relative error = %e', ...
        data(end,1), data(end,2));
      
      % refresh gui
      if mod(length(data(:,1)), gui.refreshFrequency) == 0
        drawnow;
      end
      
    end
    
    function newNeutralityCycleIteration(gui, chemistry, ~)
      
      % add data of the new iteration to the array neutralityErrorData
      if isempty(gui.neutralityErrorData)
        data = [chemistry.neutralityIterationCurrent chemistry.neutralityRelErrorCurrent];
      else
        data = [gui.neutralityErrorData; chemistry.neutralityIterationCurrent chemistry.neutralityRelErrorCurrent];
      end
      gui.neutralityErrorData = data;
      
      % Separate positive and negative errors
      pos_indices = data(:,2) >= 0;
      neg_indices = data(:,2) < 0;
      
      % Plot positive errors in blue
      h1 = plot(gui.neutralityErrorPlot, data(pos_indices,1), abs(data(pos_indices,2)), 'o', ...
        'MarkerEdgeColor', 'blue', 'DisplayName', 'Positive Errors');
      hold(gui.neutralityErrorPlot, 'on');
      
      % Plot negative errors in red
      h2 = plot(gui.neutralityErrorPlot, data(neg_indices,1), abs(data(neg_indices,2)), 'o', ...
        'MarkerEdgeColor', 'red', 'DisplayName', 'Negative Errors');
      
      % Plot all errors with a connecting line in black
      plot(gui.neutralityErrorPlot, data(:,1), abs(data(:,2)), '-', 'Color', 'black');

      % Add legend
      legend(gui.neutralityErrorPlot, [h1, h2], 'Location', 'northeast');

      % Update text in the chemistry iterations tab
      neutralityErrorInfoStr = cell(1,2);
      neutralityErrorInfoStr{1} = sprintf('iteration: %d', data(end,1));
      neutralityErrorInfoStr{2} = sprintf('error: %g', data(end,2));
      switch class(chemistry.electronKinetics)
        case 'Boltzmann'
          neutralityErrorInfoStr{3} = sprintf('E/N: %g(Td)', chemistry.workCond.reducedField);
        case 'PrescribedEedf'
          neutralityErrorInfoStr{3} = sprintf('Te: %g(eV)', chemistry.workCond.electronTemperature);
      end
      set(gui.neutralityErrorInfo, 'String', neutralityErrorInfoStr);
      
      % add info to log tab
      gui.logInfo.String{end+1} = sprintf('\t- New neutrality cycle iteration (%d): relative error = %e', ...
        data(end,1), data(end,2));

      % refresh gui
      drawnow;
      
    end
    
    function newGlobalCycleIteration(gui, chemistry, ~)
      
      % add data of the new iteration to the array neutralityErrorData
      if isempty(gui.globalErrorData)
        data = [chemistry.globalIterationCurrent chemistry.globalRelErrorCurrent];
      else
        data = [gui.globalErrorData; chemistry.globalIterationCurrent chemistry.globalRelErrorCurrent];
      end
      gui.globalErrorData = data;
      
      % Separate positive and negative errors
      pos_indices = data(:,2) >= 0;
      neg_indices = data(:,2) < 0;
      
      % Plot positive errors in blue
      h1 = plot(gui.globalErrorPlot, data(pos_indices,1), abs(data(pos_indices,2)), 'o', ...
        'MarkerEdgeColor', 'blue', 'DisplayName', 'Positive Errors');
      hold(gui.globalErrorPlot, 'on');
      
      % Plot negative errors in red
      h2 = plot(gui.globalErrorPlot, data(neg_indices,1), abs(data(neg_indices,2)), 'o', ...
        'MarkerEdgeColor', 'red', 'DisplayName', 'Negative Errors');
      
      % Plot all errors with a connecting line in black
      plot(gui.globalErrorPlot, data(:,1), abs(data(:,2)), '-', 'Color', 'black');
      
      % Add legend
      legend(gui.globalErrorPlot, [h1, h2], 'Location', 'northeast');
      
      % update text in the chemistry iterations tab
      globalErrorInfoStr = cell(1,2);
      globalErrorInfoStr{1} = sprintf('iteration: %d', data(end,1));
      globalErrorInfoStr{2} = sprintf('error: %g', data(end,2));
      set(gui.globalErrorInfo, 'String', globalErrorInfoStr);
      
      % add info to log tab
      gui.logInfo.String{end+1} = sprintf('\t- New global cycle iteration (%d): relative error = %e', ...
        data(end,1), data(end,2));

      % refresh gui
      drawnow;
      
    end
    
    function newChemistrySolution(gui, chemistry, ~)
      
      % evaluate new solution ID
      newSolutionID = length(gui.solutions)+1;
      
      % save chemistry solutions for later use on the gui
      gui.solutions(newSolutionID).time = chemistry.solution.time;
      gui.solutions(newSolutionID).finalDischargeDensity = chemistry.solution.finalDischargeDensity;
      gui.solutions(newSolutionID).densitiesTime = chemistry.solution.densitiesTime;
      gui.solutions(newSolutionID).gasTemperatureTime = chemistry.solution.gasTemperatureTime;
      
      if isfield(chemistry.solution, 'nearWallTemperatureTime')
        gui.solutions(newSolutionID).nearWallTemperatureTime = chemistry.solution.nearWallTemperatureTime;
      end
      if isfield(chemistry.solution, 'wallTemperatureTime')
        gui.solutions(newSolutionID).wallTemperatureTime = chemistry.solution.wallTemperatureTime;
      end
      gui.solutions(newSolutionID).reactionsInfo = chemistry.solution.reactionsInfo;
      gui.solutions(newSolutionID).workCond = chemistry.workCond.struct;
      gui.solutions(newSolutionID).targetGasPressure = chemistry.targetGasPressure;
      gui.solutions(newSolutionID).excitationParameter = chemistry.convergenceParameter;
      gui.solutions(newSolutionID).neutralityRelError = chemistry.neutralityRelErrorCurrent;

      % save associated electron Kinetic solution for later use on the gui (in case it is activated)
      if ~isempty(chemistry.electronKinetics)
        switch class(chemistry.electronKinetics)
          case 'Boltzmann'
            gui.solutions(newSolutionID).reducedField = chemistry.workCond.reducedField;
            gui.solutions(newSolutionID).firstAnisotropy = chemistry.electronKinetics.firstAnisotropy;
          case 'PrescribedEedf'
            gui.solutions(newSolutionID).electronTemperature = chemistry.workCond.electronTemperature;
        end
        gui.solutions(newSolutionID).eedf = chemistry.electronKinetics.eedf;
        gui.solutions(newSolutionID).energyValues = chemistry.electronKinetics.energyGrid.cell;
        gui.solutions(newSolutionID).power = chemistry.electronKinetics.power;
        gui.solutions(newSolutionID).swarmParam = chemistry.electronKinetics.swarmParam;
        gui.solutions(newSolutionID).rateCoeffAll = chemistry.electronKinetics.rateCoeffAll;
        gui.solutions(newSolutionID).rateCoeffExtra = chemistry.electronKinetics.rateCoeffExtra;
      end
      
      % add new entry to heavySpeciesDensitiesPopUpMenu
      if length(gui.chemistryVsTimeRunPopUpMenu.String) == 1
        newString = [sprintf('Run %d', newSolutionID); gui.chemistryVsTimeRunPopUpMenu.String];
        set(gui.chemistryVsTimeRunPopUpMenu, 'String', newString);
        % add new entry to eedfPopUpMenu (in case electron kinetics is activated)
        if ~isempty(chemistry.electronKinetics)
          set(gui.eedfPopUpMenu, 'String', newString);
        end
      end
      
      % add new entry to resultsTextPopUpMenu
      if chemistry.isPulsed
        newString = [sprintf('t = %g s', chemistry.workCond.currentTime); gui.resultsTextPopUpMenu.String];
        set(gui.resultsTextPopUpMenu, 'String', newString);
      else
        set(gui.resultsTextPopUpMenu, 'String', gui.chemistryVsTimeRunPopUpMenu.String(1:end-1));
      end
      
      % update graphical and text result panels with the new solution
      if chemistry.isPulsed && chemistry.workCond.currentTime == chemistry.solution.time(end)
        gui.addChemistryPlots(newSolutionID,1:length(gui.chemistryGasArray),1:length(gui.chemistryStateArray));
      elseif ~chemistry.isPulsed
        gui.addChemistryPlots(newSolutionID,1:length(gui.chemistryGasArray),1:length(gui.chemistryStateArray));
      end
      gui.updateFinalDensitiesInfo(newSolutionID);
      gui.updateFinalTemperaturesInfo(newSolutionID);
      gui.updateChemParametersInfo(newSolutionID);
      gui.updateFinalBalanceInfo(newSolutionID);
      if ~isempty(chemistry.electronKinetics)
        gui.addEedfPlot(newSolutionID, 1);
        gui.updatePowerBalanceInfo(newSolutionID);
        gui.updateSwarmParamInfo(newSolutionID);
        gui.updateRateCoeffInfo(newSolutionID);
      end
      
      % refresh gui
      drawnow;
      
    end
    
    function addChemistryPlots(gui, solutionID, gasIDArray, stateIDArray)
      
      % add new plot
      time = gui.solutions(solutionID).time;
      densitiesTime = gui.solutions(solutionID).densitiesTime;
      gasTemperatureTime = gui.solutions(solutionID).gasTemperatureTime;
      if isfield(gui.solutions(solutionID), 'nearWallTemperatureTime')
        nearWallTemperatureTime = gui.solutions(solutionID).nearWallTemperatureTime;
      else
        nearWallTemperatureTime = [];
      end
      if isfield(gui.solutions(solutionID), 'wallTemperatureTime')
        wallTemperatureTime = gui.solutions(solutionID).wallTemperatureTime;
      else
        wallTemperatureTime = [];
      end
      
      % organise states as a family tree
      familyIDs = [];
      for gas = gui.chemistryGasArray
        if any(gas.ID == gasIDArray)
          for eleState = gas.stateArray
            if strcmp(eleState.type, 'ele')
              if any(eleState.ID == stateIDArray)
                familyIDs(end+1) = eleState.ID;
                for vibState = eleState.childArray
                  familyIDs(end+1) = vibState.ID;
                  for rotState = vibState.childArray
                    familyIDs(end+1) = rotState.ID;
                  end
                end
              else
                for vibState = eleState.childArray
                  if any(vibState.ID == stateIDArray)
                    familyIDs(end+1) = vibState.ID;
                    for rotState = vibState.childArray
                      familyIDs(end+1) = rotState.ID;
                    end
                  else
                    for rotState = vibState.childArray
                      if any(rotState.ID == stateIDArray)
                        familyIDs(end+1) = rotState.ID;
                      end
                    end
                  end
                end
              end
            end
          end
          for ionState = gas.stateArray
            if strcmp(ionState.type, 'ion') && any(ionState.ID == stateIDArray)
              familyIDs(end+1) = ionState.ID;
            end
          end
        end
      end
      
      % plot results
      yyaxis(gui.chemistryVsTimePlot, 'left');
      for i = 1:length(familyIDs)
        plot(gui.chemistryVsTimePlot, time, densitiesTime(:,familyIDs(i)), '-', 'Color', [rand rand rand], 'Tag', ...
          'chemistryVsTimePlot');
      end
      
      % plot gas temperature in case it is not plotted already 
      if ~gui.isGasTemperaturePlotted
        yyaxis(gui.chemistryVsTimePlot, 'right')
        plot(gui.chemistryVsTimePlot, time, gasTemperatureTime, '--', 'Color', 'r', 'LineWidth', 2);
        if ~isempty(nearWallTemperatureTime)
          plot(gui.chemistryVsTimePlot, time, nearWallTemperatureTime, '--', 'Color', 'b', 'LineWidth', 2);
        end
        if ~isempty(wallTemperatureTime)
          plot(gui.chemistryVsTimePlot, time, wallTemperatureTime, '--', 'Color', 'k', 'LineWidth', 2);
        end
      end
      
      % add info to legend
      if ~gui.isGasTemperaturePlotted
        legendStr = gui.chemistryVsTimeLegend;
        for i = 1:length(familyIDs)
          legendStr{end+1} = gui.chemistryStateArray(familyIDs(i)).name;
        end
        legendStr{end+1} = 'Gas temperature';
        if ~isempty(nearWallTemperatureTime)
          legendStr{end+1} = 'Near wall temperature';
        end
        if ~isempty(wallTemperatureTime)
          legendStr{end+1} = 'Wall temperature';
        end
        % set flag indicating if the gas temperature is plotted to true
        gui.isGasTemperaturePlotted = true;
      else
        if isempty(wallTemperatureTime)
          if isempty(nearWallTemperatureTime)
            legendStr = {gui.chemistryVsTimeLegend{1:end-1}};
          else
            legendStr = {gui.chemistryVsTimeLegend{1:end-2}};
          end
        else
          legendStr = {gui.chemistryVsTimeLegend{1:end-3}};
        end
        for i = 1:length(familyIDs)
          legendStr{end+1} = gui.chemistryStateArray(familyIDs(i)).name;
        end
        legendStr{end+1} = 'Gas temperature';
        if ~isempty(nearWallTemperatureTime)
          legendStr{end+1} = 'Near wall temperature';
        end
        if ~isempty(wallTemperatureTime)
          legendStr{end+1} = 'Wall temperature';
        end
      end
      
      % draw and save legend info
      legend(gui.chemistryVsTimePlot, legendStr);
      gui.chemistryVsTimeLegend = legendStr;
      
    end
    
    function chemistryVsTimeRunPopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution(s) to plot
      solutionIDArray = gui.chemistryVsTimeRunPopUpMenu.Value;
      if solutionIDArray == length(gui.solutions)+1
        gui.clearChemistryVsTimePlot();
        solutionIDArray = 1:length(gui.solutions);
      end
    
      % plot selected solution(s)
      for solutionID = solutionIDArray
        gui.addChemistryPlots(solutionID,1:length(gui.chemistryGasArray),1:length(gui.chemistryStateArray));
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function chemistryVsTimeGasPopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution(s) to plot
      solutionIDArray = gui.chemistryVsTimeRunPopUpMenu.Value;
      if solutionIDArray == length(gui.solutions)+1
        gui.clearChemistryVsTimePlot();
        solutionIDArray = 1:length(gui.solutions);
      end
      
      % evaluate gas(s) to plot
      gasIDArray = gui.chemistryVsTimeGasPopUpMenu.Value;
      if gasIDArray == length(gui.chemistryGasArray)+1
        gui.clearChemistryVsTimePlot();
        gasIDArray = 1:length(gui.chemistryGasArray);
      end
    
      % plot selected solution(s)
      for solutionID = solutionIDArray
        gui.addChemistryPlots(solutionID, gasIDArray, 1:length(gui.chemistryStateArray));
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function chemistryVsTimeStatePopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution(s) to plot
      solutionIDArray = gui.chemistryVsTimeRunPopUpMenu.Value;
      if solutionIDArray == length(gui.solutions)+1
        gui.clearChemistryVsTimePlot();
        solutionIDArray = 1:length(gui.solutions);
      end
      
      % evaluate state(s) to plot
      stateIDArray = gui.chemistryVsTimeStatePopUpMenu.Value;
      if stateIDArray == length(gui.chemistryStateArray)+1
        gui.clearChemistryVsTimePlot();
        stateIDArray = 1:length(gui.chemistryStateArray);
      end
    
      % plot selected solution(s)
      for solutionID = solutionIDArray
        gui.addChemistryPlots(solutionID, 1:length(gui.chemistryGasArray), stateIDArray);
      end
      
      % refresh gui
      drawnow;
    
    end
    
    function clearChemistryVsTimePlot(gui, ~, ~)

      % clear plot
      yyaxis(gui.chemistryVsTimePlot, 'left');
      cla(gui.chemistryVsTimePlot);
      yyaxis(gui.chemistryVsTimePlot, 'right');
      cla(gui.chemistryVsTimePlot);
      gui.isGasTemperaturePlotted = false;
      % clear legend
      legend(gui.chemistryVsTimePlot, 'off');
      gui.chemistryVsTimeLegend = cell.empty;
      
      % refresh gui
      drawnow;
      
    end
    
    function updateFinalDensitiesInfo(gui, solutionID)
      
      % save local copy of different variables
      finalDensities = gui.solutions(solutionID).finalDischargeDensity;
      reactionRates = [gui.solutions(solutionID).reactionsInfo.netRate];
      % evaluate number of gases and species
      numberOfGases = length(gui.chemistryGasArray);
      numberOfSpecies = length(finalDensities);
      % determine total gas density and creation/destruction rates (final time)
      gasDensities = zeros(1,numberOfGases);
      totalVolumeGasDensity = 0;
      totalSurfaceSiteDensity = 0;
      rateBalances = zeros(1,numberOfSpecies);
      for gas = gui.chemistryGasArray
        for state = gas.stateArray
          if strcmp(state.type, 'ele') || strcmp(state.type, 'ion')
            gasDensities(gas.ID) = gasDensities(gas.ID) + finalDensities(state.ID);
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
      % evaluate length of the 'Species' column
      speciesColumnLength = 13;
      for gas = gui.chemistryGasArray
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
      % write volume chemistry information (final densities, final populations and final particle balances)
      finalDensitiesStr = cell.empty;
      finalDensitiesStr{end+1} = sprintf('*****************************');
      finalDensitiesStr{end+1} = sprintf('*  Chemistry (Volume phase) *');
      finalDensitiesStr{end+1} = sprintf('*****************************');
      finalDensitiesStr{end+1} = sprintf('Species%s Abs. Density  Population    Balance', ...
        repmat(' ', 1, speciesColumnLength-7));
      finalDensitiesStr{end+1} = sprintf('%s (m^-3)', repmat(' ', 1, speciesColumnLength));
      finalDensitiesStr{end+1} = repmat('-', 1, speciesColumnLength+43);
      for gas = gui.chemistryGasArray
        if gas.isSurfaceSpecies
          continue
        end
        finalDensitiesStr{end+1} = sprintf('%s[%f%%]', gas.name, 100*gasDensities(gas.ID)/totalVolumeGasDensity);
        for eleState = gas.stateArray
          if strcmp(eleState.type, 'ele')
            if isempty(eleState.childArray)
              finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e     %+9.3e'], eleState.name, ...
                finalDensities(eleState.ID), finalDensities(eleState.ID)/gasDensities(gas.ID), ...
                rateBalances(eleState.ID));
            else
              finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e'], eleState.name, ...
                finalDensities(eleState.ID), finalDensities(eleState.ID)/gasDensities(gas.ID));
            end
            for vibState = eleState.childArray
              if isempty(vibState.childArray)
                finalDensitiesStr{end+1} = sprintf([auxStr2 '| %9.3e   | %9.3e   | %+9.3e'], vibState.name, ...
                  finalDensities(vibState.ID), finalDensities(vibState.ID)/finalDensities(eleState.ID), ...
                  rateBalances(vibState.ID));
              else
                finalDensitiesStr{end+1} = sprintf([auxStr2 '| %9.3e   | %9.3e'], vibState.name, ...
                  finalDensities(vibState.ID), finalDensities(vibState.ID)/finalDensities(eleState.ID));
              end
              for rotState = vibState.childArray
                finalDensitiesStr{end+1} = sprintf([auxStr3 '| | %9.3e | | %9.3e | | %+9.3e'], rotState.name, ...
                  finalDensities(rotState.ID), finalDensities(rotState.ID)/finalDensities(vibState.ID), ...
                  rateBalances(rotState.ID));
              end
            end
          end
        end
        for ionState = gas.stateArray
          if strcmp(ionState.type, 'ion')
            finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e     %+9.3e'], ionState.name, ...
              finalDensities(ionState.ID), finalDensities(ionState.ID)/gasDensities(gas.ID), rateBalances(ionState.ID));
          end
        end
      end
      % write surface chemistry information (final densities, final populations and final particle balances)
      if totalSurfaceSiteDensity ~= 0
        finalDensitiesStr{end+2} = sprintf('*****************************');
        finalDensitiesStr{end+1} = sprintf('* Chemistry (Surface phase) *');
        finalDensitiesStr{end+1} = sprintf('*****************************');
        finalDensitiesStr{end+1} = sprintf('Species%s Abs. Density  Population    Balance', ...
          repmat(' ', 1, speciesColumnLength-7));
        finalDensitiesStr{end+1} = sprintf('%s (m^-2)', repmat(' ', 1, speciesColumnLength));
        finalDensitiesStr{end+1} = repmat('-', 1, speciesColumnLength+43);
        for gas = gui.chemistryGasArray
          if gas.isVolumeSpecies
            continue
          end
          finalDensitiesStr{end+1} = sprintf('%s[%f%%]', gas.name, 100*gasDensities(gas.ID)/totalSurfaceSiteDensity);
          for eleState = gas.stateArray
            if strcmp(eleState.type, 'ele')
              if isempty(eleState.childArray)
                finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e     %+9.3e'], eleState.name, ...
                  finalDensities(eleState.ID), finalDensities(eleState.ID)/gasDensities(gas.ID), ...
                  rateBalances(eleState.ID));
              else
                finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e'], eleState.name, ...
                  finalDensities(eleState.ID), finalDensities(eleState.ID)/gasDensities(gas.ID));
              end
              for vibState = eleState.childArray
                if isempty(vibState.childArray)
                  finalDensitiesStr{end+1} = sprintf([auxStr2 '| %9.3e   | %9.3e   | %+9.3e'], vibState.name, ...
                    finalDensities(vibState.ID), finalDensities(vibState.ID)/finalDensities(eleState.ID), ...
                    rateBalances(vibState.ID));
                else
                  finalDensitiesStr{end+1} = sprintf([auxStr2 '| %9.3e   | %9.3e'], vibState.name, ...
                    finalDensities(vibState.ID), finalDensities(vibState.ID)/finalDensities(eleState.ID));
                end
                for rotState = vibState.childArray
                  finalDensitiesStr{end+1} = sprintf([auxStr3 '| | %9.3e | | %9.3e | | %+9.3e'], rotState.name, ...
                    finalDensities(rotState.ID), finalDensities(rotState.ID)/finalDensities(vibState.ID), ...
                    rateBalances(rotState.ID));
                end
              end
            end
          end
          for ionState = gas.stateArray
            if strcmp(ionState.type, 'ion')
              finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e     %9.3e     %+9.3e'], ionState.name, ...
                finalDensities(ionState.ID), finalDensities(ionState.ID)/gasDensities(gas.ID), rateBalances(ionState.ID));
            end
          end
        end
      end
      % write eedf information (final populations)
      if isfield(gui.solutions(solutionID), 'eedf')
        % evaluate length of the 'Species' column
        speciesColumnLength = 13;
        for gas = gui.eedfGasArray
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
        finalDensitiesStr{end+2} = sprintf('*********************');
        finalDensitiesStr{end+1} = sprintf('* Electron Kinetics *');
        finalDensitiesStr{end+1} = sprintf('*********************');
        finalDensitiesStr{end+1} = sprintf('Species%s Population', repmat(' ', 1, speciesColumnLength-7));
        finalDensitiesStr{end+1} = repmat('-', 1, speciesColumnLength+14);
        for gas = gui.eedfGasArray
          finalDensitiesStr{end+1} = sprintf('%s[%f%%]', gas.name, 100*gas.fraction);
          for eleState = gas.stateArray
            if strcmp(eleState.type, 'ele')
              finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e      %9.3e'], eleState.name, eleState.population);
              if ~isempty(eleState.childArray)
                for vibState = eleState.childArray
                  finalDensitiesStr{end+1} = sprintf([auxStr2 '| %9.3e'], vibState.name, vibState.population);
                  if ~isempty(vibState.childArray)
                    for rotState = vibState.childArray
                      finalDensitiesStr{end+1} = sprintf([auxStr3 '| | %9.3e'], rotState.name, rotState.population);
                    end
                  end
                end
              end
            end
          end
          for ionState = gas.stateArray
            if strcmp(ionState.type, 'ion')
              finalDensitiesStr{end+1} = sprintf([auxStr1 '%9.3e'], ionState.name, ionState.population);
            end
          end
        end
      end
      
      % update data in the GUI
      set(gui.finalDensitiesInfo, 'String', finalDensitiesStr);
      
    end
    
    function updateFinalTemperaturesInfo(gui, solutionID)
      
      % write temperature information in cell array
      finalTemperaturesStr = cell.empty;
      finalTemperaturesStr{end+1} = sprintf('      Gas temperature = %9.3e (K)', ...
        gui.solutions(solutionID).workCond.gasTemperature);
      if ~isempty(gui.solutions(solutionID).workCond.nearWallTemperature)
        finalTemperaturesStr{end+1} = sprintf('Near wall temperature = %9.3e (K)', ...
          gui.solutions(solutionID).workCond.nearWallTemperature);
      end
      if ~isempty(gui.solutions(solutionID).workCond.wallTemperature)
        finalTemperaturesStr{end+1} = sprintf('     Wall temperature = %9.3e (K)', ...
          gui.solutions(solutionID).workCond.wallTemperature);
      end
      if ~isempty(gui.solutions(solutionID).workCond.extTemperature)
        finalTemperaturesStr{end+1} = sprintf(' External temperature = %9.3e (K)', ...
          gui.solutions(solutionID).workCond.extTemperature);
      end
      
      % update data in the GUI
      set(gui.finalTemperaturesInfo, 'String', finalTemperaturesStr);
      
    end

    function updateChemParametersInfo(gui, solutionID)

      % write chemistry parameters information in cell array
      chemParametersStr = cell.empty;
      chemParametersStr{end+1} = sprintf('          Reduced electric field = %9.3e (Td)', ...
        gui.solutions(solutionID).workCond.reducedField);
      chemParametersStr{end+1} = sprintf('             Target gas pressure = %9.3e (Pa)', ...
        gui.solutions(solutionID).targetGasPressure);
      chemParametersStr{end+1} = sprintf('              Final gas pressure = %9.3e (Pa)', ...
        gui.solutions(solutionID).workCond.gasPressure);
      chemParametersStr{end+1} = sprintf('           Final gas temperature = %9.3e (K)', ...
        gui.solutions(solutionID).workCond.gasTemperature);
      chemParametersStr{end+1} = sprintf('               Final gas density = %9.3e (m-3)', ...
        gui.solutions(solutionID).workCond.gasDensity);
      chemParametersStr{end+1} = sprintf('                Total input flow = %9.3e (sccm)', ...
        gui.solutions(solutionID).workCond.totalSccmInFlow);
      chemParametersStr{end+1} = sprintf('               Total output flow = %9.3e (sccm)', ...
        gui.solutions(solutionID).workCond.totalSccmOutFlow);      

      % save local copy of the solution
      excitationParameter = gui.solutions(solutionID).excitationParameter;

      % continue writing chemistry parameters information in cell array
      if strcmp(excitationParameter, 'electronDensity')
        chemParametersStr{end+1} = sprintf('         Target electron density = %9.3e (m-3)', ...
          gui.solutions(solutionID).workCond.electronDensity);
        chemParametersStr{end+1} = sprintf('           Neutrality rel. error = %9.3e', ...
          gui.solutions(solutionID).neutralityRelError);
      elseif strcmp(excitationParameter, 'dischargeCurrent')
        chemParametersStr{end+1} = sprintf('          Final electron density = %9.3e (m-3)', ...
          gui.solutions(solutionID).workCond.electronDensity);
        chemParametersStr{end+1} = sprintf('           Neutrality rel. error = %9.3e', ...
          gui.solutions(solutionID).neutralityRelError);
        % local calculation and writing of discharge current solution and relative error 
        targetDischargeCurrent = gui.solutions(solutionID).workCond.dischargeCurrent;
        finalDischargeCurrent = Constant.electronCharge * gui.solutions(solutionID).workCond.electronDensity * ...
            gui.solutions(solutionID).swarmParam.driftVelocity * pi * ...
            gui.solutions(solutionID).workCond.chamberRadius^2; 
        dischargeCurrentRelError = (targetDischargeCurrent - finalDischargeCurrent) / targetDischargeCurrent;
        chemParametersStr{end+1} = sprintf('        Target discharge current = %9.3e (A)', targetDischargeCurrent);
        chemParametersStr{end+1} = sprintf('         Final discharge current = %9.3e (A)', finalDischargeCurrent);
        chemParametersStr{end+1} = sprintf('    Discharge current rel. error = %9.3e', dischargeCurrentRelError);
      elseif strcmp(excitationParameter, 'dischargePowerDensity')
        chemParametersStr{end+1} = sprintf('          Final electron density = %9.3e (m-3)', ...
          gui.solutions(solutionID).workCond.electronDensity);
        chemParametersStr{end+1} = sprintf('           Neutrality rel. error = %9.3e', ...
          gui.solutions(solutionID).neutralityRelError);
        % local calculation and writing of discharge power density solution and relative error 
        targetDischargePowerDensity = gui.solutions(solutionID).workCond.dischargePowerDensity;
        finalDischargePowerDensity = gui.solutions(solutionID).workCond.electronDensity * ...
            gui.solutions(solutionID).workCond.gasDensity * ...
              gui.solutions(solutionID).power.field * Constant.electronCharge;
        dischargePowerDensityRelError = (targetDischargePowerDensity - finalDischargePowerDensity) / ...
              targetDischargePowerDensity;
        chemParametersStr{end+1} = sprintf('  Target discharge power density = %9.3e (W m-3)', targetDischargePowerDensity);
        chemParametersStr{end+1} = sprintf('   Final discharge power density = %9.3e (W m-3)', finalDischargePowerDensity);
        chemParametersStr{end+1} = sprintf('Discharge power dens. rel. error = %9.3e', dischargePowerDensityRelError);
      end  
                
      % update data in the GUI
      set(gui.chemParametersInfo, 'String', chemParametersStr);
      
    end

    function updateFinalBalanceInfo(gui, solutionID)
      
      % evaluate maximum length of reaction descriptions
      maxReactionLength = 0;
      for reaction = gui.chemistryReactionArray
        maxReactionLength = max(maxReactionLength, length(reaction.description));
      end
      auxStr = sprintf('      %%-%ds ', maxReactionLength);
      % save local copy of the solution
      reactionRates = [gui.solutions(solutionID).reactionsInfo.netRate];
      % write detailed particle balance
      finalBalanceStr = cell.empty;
      for gas = gui.chemistryGasArray
        finalBalanceStr{end+1} = repmat('*', 1, 33+length(gas.name));
        finalBalanceStr{end+1} = sprintf('* Particle balance for %s species *', gas.name);
        finalBalanceStr{end+1} = repmat('*', 1, 33+length(gas.name));
        if gas.isVolumeSpecies
          rateUnitsStr = 'm^-3s^-1';
          rateRenorm = 1;
        else
          rateUnitsStr = 'm^-2s^-1';
          rateRenorm = gui.solutions(solutionID).workCond.areaOverVolume;
        end
        for eleState = gas.stateArray
          if strcmp(eleState.type, 'ele')
            if isempty(eleState.childArray)
              finalBalanceStr{end+1} = sprintf('-> Particle balance for %s:', eleState.name);
              % evaluate creation channels 
              finalBalanceStr{end+1} = sprintf('    * Reactions where %s is created:', eleState.name);
              finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
                rateUnitsStr);
              reactions = eleState.reactionsCreation;
              rates = zeros(1,length(reactions));
              for i = 1:length(reactions)
                for j = 1:length(reactions(i).productArray)
                  if eleState.ID == reactions(i).productArray(j).ID
                    stoiCoeff = reactions(i).productStoiCoeff(j);
                    break;
                  end
                end
                rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
              end
              totalCreationRate = sum(rates(:));
              for i = 1:length(reactions)
                finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                  rates(i), rates(i)*100/totalCreationRate);
              end
              finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', repmat(' ', 1, 1+maxReactionLength), ...
                totalCreationRate, 100);
              % evaluate destruction channels
              finalBalanceStr{end+1} = sprintf('    * Reactions where %s is destroyed:', eleState.name);
              finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
                rateUnitsStr);
              reactions = eleState.reactionsDestruction;
              rates = zeros(1,length(reactions));
              for i = 1:length(reactions)
                for j = 1:length(reactions(i).reactantArray)
                  if eleState.ID == reactions(i).reactantArray(j).ID
                    stoiCoeff = reactions(i).reactantStoiCoeff(j);
                    break;
                  end
                end
                rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
              end
              totalDestructionRate = sum(rates(:));
              for i = 1:length(reactions)
                finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                  rates(i), rates(i)*100/totalDestructionRate);
              end
              finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', repmat(' ', 1, 1+maxReactionLength), ...
                totalDestructionRate, 100);
              % evaluate species balance
              finalBalanceStr{end+1} = sprintf('');
              finalBalanceStr{end+1} = sprintf('    * Relative %s balance (creation-destruction)/creation: %9.3e%%', ...
                eleState.name, (totalCreationRate-totalDestructionRate)*100/totalCreationRate);
              finalBalanceStr{end+1} = sprintf('');
            else
              for vibState = eleState.childArray
                if isempty(vibState.childArray)
                  finalBalanceStr{end+1} = sprintf('-> Particle balance for %s:', vibState.name);
                  % evaluate creation channels
                  finalBalanceStr{end+1} = sprintf('    * Reactions where %s is created:', vibState.name);
                  finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
                    rateUnitsStr);
                  reactions = vibState.reactionsCreation;
                  rates = zeros(1,length(reactions));
                  for i = 1:length(reactions)
                    for j = 1:length(reactions(i).productArray)
                      if vibState.ID == reactions(i).productArray(j).ID
                        stoiCoeff = reactions(i).productStoiCoeff(j);
                        break;
                      end
                    end
                    rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
                  end
                  totalCreationRate = sum(rates(:));
                  for i = 1:length(reactions)
                    finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                      rates(i), rates(i)*100/totalCreationRate);
                  end
                  finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', ...
                    repmat(' ', 1, 1+maxReactionLength), totalCreationRate, 100);
                  % evaluate destruction channels
                  finalBalanceStr{end+1} = sprintf('    * Reactions where %s is destroyed:', vibState.name);
                  finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
                    rateUnitsStr);
                  reactions = vibState.reactionsDestruction;
                  rates = zeros(1,length(reactions));
                  for i = 1:length(reactions)
                    for j = 1:length(reactions(i).reactantArray)
                      if vibState.ID == reactions(i).reactantArray(j).ID
                        stoiCoeff = reactions(i).reactantStoiCoeff(j);
                        break;
                      end
                    end
                    rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
                  end
                  totalDestructionRate = sum(rates(:));
                  for i = 1:length(reactions)
                    finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                      rates(i), rates(i)*100/totalDestructionRate);
                  end
                  finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', ...
                    repmat(' ', 1, 1+maxReactionLength), totalDestructionRate, 100);
                  % evaluate species balance
                  finalBalanceStr{end+1} = sprintf('');
                  finalBalanceStr{end+1} = sprintf(['    * Relative %s balance (creation-destruction)/creation: ' ...
                    '%9.3e%%'], vibState.name, (totalCreationRate-totalDestructionRate)*100/totalCreationRate);
                  finalBalanceStr{end+1} = sprintf('');
                else
                  for rotState = vibState.childArray
                    finalBalanceStr{end+1} = sprintf('-> Particle balance for %s:', rotState.name);
                    % evaluate creation channels
                    finalBalanceStr{end+1} = sprintf('    * Reactions where %s is created:', rotState.name);
                    finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', ...
                      repmat(' ', 1, maxReactionLength+7), rateUnitsStr);
                    reactions = rotState.reactionsCreation;
                    rates = zeros(1,length(reactions));
                    for i = 1:length(reactions)
                      for j = 1:length(reactions(i).productArray)
                        if rotState.ID == reactions(i).productArray(j).ID
                          stoiCoeff = reactions(i).productStoiCoeff(j);
                          break;
                        end
                      end
                      rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
                    end
                    totalCreationRate = sum(rates(:));
                    for i = 1:length(reactions)
                      finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                        rates(i), rates(i)*100/totalCreationRate);
                    end
                    finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', ...
                      repmat(' ', 1, 1+maxReactionLength), totalCreationRate, 100);
                    % evaluate destruction channels
                    finalBalanceStr{end+1} = sprintf('    * Reactions where %s is destroyed:', rotState.name);
                    finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', ...
                      repmat(' ', 1, maxReactionLength+7), rateUnitsStr);
                    reactions = rotState.reactionsDestruction;
                    rates = zeros(1,length(reactions));
                    for i = 1:length(reactions)
                      for j = 1:length(reactions(i).reactantArray)
                        if rotState.ID == reactions(i).reactantArray(j).ID
                          stoiCoeff = reactions(i).reactantStoiCoeff(j);
                          break;
                        end
                      end
                      rates(i) = stoiCoeff*reactionRates(reactions(i).ID)/rateRenorm;
                    end
                    totalDestructionRate = sum(rates(:));
                    for i = 1:length(reactions)
                      finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                        rates(i), rates(i)*100/totalDestructionRate);
                    end
                    finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', ...
                      repmat(' ', 1, 1+maxReactionLength), totalDestructionRate, 100);
                    % evaluate species balance
                    finalBalanceStr{end+1} = sprintf('');
                    finalBalanceStr{end+1} = sprintf(['    * Relative %s balance (creation-destruction)/creation: ' ...
                      '%9.3e%%'], rotState.name, (totalCreationRate-totalDestructionRate)*100/totalCreationRate);
                    finalBalanceStr{end+1} = sprintf('');
                  end
                end
              end
            end
          end
        end
        for ionState = gas.stateArray
          if strcmp(ionState.type, 'ion')
            % escribir creacion-destruccion
            finalBalanceStr{end+1} = sprintf('-> Particle balance for %s:', ionState.name);
            % evaluate creation channels
            finalBalanceStr{end+1} = sprintf('    * Reactions where %s is created:', ionState.name);
            finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
              rateUnitsStr);
            reactions = ionState.reactionsCreation;
            rates = zeros(1,length(reactions));
            for i = 1:length(reactions)
              rates(i) = reactionRates(reactions(i).ID)/rateRenorm;
            end
            totalCreationRate = sum(rates(:));
            for i = 1:length(reactions)
              finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                rates(i), rates(i)*100/totalCreationRate);
            end
            finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', repmat(' ', 1, 1+maxReactionLength), ...
              totalCreationRate, 100);
            % evaluate destruction channels
            finalBalanceStr{end+1} = sprintf('    * Reactions where %s is destroyed:', ionState.name);
            finalBalanceStr{end+1} = sprintf('%sRate(%s) Contribution', repmat(' ', 1, maxReactionLength+7), ...
              rateUnitsStr);
            reactions = ionState.reactionsDestruction;
            rates = zeros(1,length(reactions));
            for i = 1:length(reactions)
              rates(i) = reactionRates(reactions(i).ID)/rateRenorm;
            end
            totalDestructionRate = sum(rates(:));
            for i = 1:length(reactions)
              finalBalanceStr{end+1} = sprintf([auxStr '%9.3e      %9.3e%%'], reactions(i).description, ...
                rates(i), rates(i)*100/totalDestructionRate);
            end
            finalBalanceStr{end+1} = sprintf('%sTOTAL %9.3e      %9.3e%%', repmat(' ', 1, 1+maxReactionLength), ...
              totalDestructionRate, 100);
            % evaluate species balance
            finalBalanceStr{end+1} = sprintf('');
            finalBalanceStr{end+1} = sprintf('    * Relative %s balance (creation-destruction)/creation: %9.3e%%', ...
              ionState.name, (totalCreationRate-totalDestructionRate)*100/totalCreationRate);
            finalBalanceStr{end+1} = sprintf('');
          end
        end
        finalBalanceStr{end+1} = '';
      end
      
      set(gui.finalBalanceInfo, 'String', finalBalanceStr);
      
    end
    
    function genericStatusMessage(gui, ~, statusEventData)

      str = statusEventData.message;
      if endsWith(str, '\n')
        str = str(1:end-2);
      end
      gui.logInfo.String{end+1} = sprintf(str);

    end

  end
  
end
