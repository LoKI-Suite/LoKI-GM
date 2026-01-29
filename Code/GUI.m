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

classdef GUI < handle
  %GUI Class that defines a Graphical User Interface
  %   Objects of this class are the GUI of the simulation. The class has methods that create and update the
  %   information that appears in the different panels and elements of the GUI
  
  properties (Access = private)
    
    handle;
    cli;
    eedfGasArray;
    collisionArray;
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

    statusPanel;
    statusTabGroup;
    logTab;
    logInfo;
    
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
            gui.evolvingParameterPopUpMenuStr = 'Te = %9.3e (eV)';
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
      
      % display the gui
      drawnow;
      
    end
    
  end
  
  methods (Access = private)
    
    function createWindow(gui)
      
      % create figure for GUI
      screenSize = get(groot,'ScreenSize');
      gui.handle = figure('name', 'LoKI Simulation Tool', 'OuterPosition', [0 35 screenSize(3) screenSize(4)-35], ...
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
    
    function resultsTextPopUpMenuHandler(gui, ~, ~)
    
      % evaluate solution to show
      solutionID = gui.resultsTextPopUpMenu.Value;
    
      % show selected solution(s)
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
      swarmStr{end+1} = sprintf('               Reduced electric field = %#.3e (Td)', reducedField);
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
    
    function genericStatusMessage(gui, ~, statusEventData)

      str = statusEventData.message;
      if endsWith(str, '\n')
        str = str(1:end-2);
      end
      gui.logInfo.String{end+1} = sprintf(str);

    end

  end
  
end
