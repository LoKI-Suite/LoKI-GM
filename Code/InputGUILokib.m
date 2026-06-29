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

classdef InputGUILokib < handle
    %InputGUILokib Class that provides a user-friendly interface for LoKI setup

    properties
        Setup;      % Struct holding the simulation setup data
        Fig;        % Main figure handle
        UIControls; % Struct to hold handles to UI controls for easy access/update
        IsApplyingDefaultSetup = false; % Guard against recursive GUI updates while resetting defaults
    end

    methods
        function gui = InputGUILokib(inputFile)
            % Constructor
            close all force;
            delete(findall(groot, 'Type', 'figure'));
            clearvars -except gui inputFile;

            % Show loading animation with progress bar while initializing
            % Center window on screen
            screenSize = get(0, 'ScreenSize');
            windowWidth = 400;
            windowHeight = 100;
            windowX = (screenSize(3) - windowWidth) / 2;
            windowY = (screenSize(4) - windowHeight) / 2;
            
            loadingFig = uifigure('Name', 'LoKI-B', 'WindowStyle', 'modal', ...
                'Position', [windowX, windowY, windowWidth, windowHeight], 'Resize', 'off', ...
                'Color', [0.94 0.94 0.94]);
            loadingGrid = uigridlayout(loadingFig, [2, 1]);
            loadingGrid.RowHeight = {'fit', 14}; % thinner progress bar
            loadingGrid.Padding = [30 20 30 16];
            loadingGrid.RowSpacing = 12;
            
            loadingLabel = uilabel(loadingGrid, 'Text', 'LoKI-B is starting...', ...
                'FontSize', 14, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
            loadingLabel.Layout.Row = 1;

            % Progress bar (rounded) using uihtml when available; fallback to panels otherwise
            useHtmlProgress = exist('uihtml', 'file') == 2;
            progressHTML = [];
            progressContainer = [];
            progressFill = [];
            containerWidth = windowWidth - 60; % fallback estimate

            if useHtmlProgress
                progressHTML = uihtml(loadingGrid);
                progressHTML.Layout.Row = 2;
                % Initial render
                progressHTML.HTMLSource = localProgressHtml(0);
                drawnow;
            else
                % Fallback: simple panels (no rounded corners available)
                progressContainer = uipanel(loadingGrid, 'BackgroundColor', [0.85 0.85 0.85], ...
                    'BorderType', 'none');
                progressContainer.Layout.Row = 2;
                progressFill = uipanel(progressContainer, 'BackgroundColor', [0.18 0.70 0.25], ...
                    'BorderType', 'none', 'Position', [0 0 0 14]);
                drawnow;
                try
                    containerWidth = progressContainer.Position(3);
                    if containerWidth == 0
                        containerWidth = windowWidth - 60;
                    end
                catch
                end
            end

            % Animate progress bar early so user sees movement immediately
            for p = 1:20
                localSetProgress(p);
                pause(0.015);
            end
            
            if nargin < 1 || isempty(inputFile)
                gui.initializeDefaultSetup(); % Initialize with defaults if no file provided
            else
                % Placeholder: Implement loading from an existing file
                try
                    % gui.Setup = gui.parseInputFile(inputFile); % You'd need to implement parseInputFile
                    gui.initializeDefaultSetup(); % For now, still use defaults
                    fprintf('Note: Loading from file not yet implemented. Using default setup.\n');
                catch ME
                    warning('Error loading input file: %s. Using default setup.', ME.message);
                    gui.initializeDefaultSetup();
                end
            end
            gui.UIControls = struct(); % Initialize empty struct for control handles
            
            % Update progress: 30% - Setup initialized
            localSetProgress(30);
            
            gui.createGUI();
            
            % Update progress: 80% - GUI created
            localSetProgress(80);
            
            gui.populateGUIFromSetup(); % Populate GUI fields with Setup data
            % Ensure run button reflects current electron kinetics state
            try
                gui.updateRunButtonState();
            catch
            end
            % Ensure run button reflects current electron kinetics state
            try
                gui.updateRunButtonState();
            catch
            end
            
            % Update progress: 100% - Complete
            localSetProgress(100);
            pause(0.1); % Brief pause to show 100%

            % Show only after everything is laid out to avoid "narrow then expand" flicker
            drawnow;
            gui.handleFigureSizeChanged();
            gui.Fig.Visible = 'on';
            drawnow; % ensure grid-managed panels have final pixel sizes
            gui.layoutTotalSccmOutFlowOverlay();
            
            % Close loading window
            close(loadingFig);

            function localSetProgress(pct)
                % pct in [0,100]
                try
                    pct = max(0, min(100, pct));
                catch
                end
                try
                    if useHtmlProgress && ~isempty(progressHTML) && isvalid(progressHTML)
                        progressHTML.HTMLSource = localProgressHtml(pct);
                        drawnow;
                        return;
                    end
                catch
                end
                % Fallback (panel fill)
                try
                    if ~isempty(progressFill) && isvalid(progressFill)
                        progressFill.Position = [0 0 (pct/100)*containerWidth 14];
                        drawnow;
                    end
                catch
                end
            end

            function html = localProgressHtml(pct)
                % Thin rounded progress bar (CSS). uihtml fills its layout cell.
                % Use fixed height with border-radius for rounded ends.
                html = sprintf([ ...
                    '<div style="width:100%%;height:8px;background:rgb(217,217,217);border-radius:999px;overflow:hidden;">' ...
                    '<div style="width:%.1f%%;height:100%%;background:rgb(46,179,67);border-radius:999px;"></div>' ...
                    '</div>'], pct);
            end
        end

        function initializeDefaultSetup(gui)
            % Initialize the Setup struct with default values (mirroring the input file)
            gui.Setup = struct();

            % Working Conditions [cite: 1, 2, 3, 4]
            gui.Setup.workingConditions.reducedField = 'logspace(-3,3,100)'; % Store as string initially
            gui.Setup.workingConditions.electronTemperature = 'linspace(0.03, 5, 100)'; % Store as string
            gui.Setup.workingConditions.excitationFrequency = 0;
            gui.Setup.workingConditions.gasPressure = 133.32;
            gui.Setup.workingConditions.gasTemperature = 300;
            gui.Setup.workingConditions.wallTemperature = 390;
            gui.Setup.workingConditions.extTemperature = 300;
            gui.Setup.workingConditions.surfaceSiteDensity = 1e19;
            gui.Setup.workingConditions.electronDensity = 1e19;
            gui.Setup.workingConditions.chamberLength = 1.0;
            gui.Setup.workingConditions.chamberRadius = 1.0;
            gui.Setup.workingConditions.totalSccmInFlow = 1;
            gui.Setup.workingConditions.totalSccmOutFlow = 'ensureIsobaric'; % Can be a number or model string
            gui.Setup.workingConditions.dischargeCurrent = 1e-2; % Optional
            gui.Setup.workingConditions.dischargePowerDensity = 1e5; % Optional

            % Electron Kinetics [cite: 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
            gui.Setup.electronKinetics.isOn = true;
            gui.Setup.electronKinetics.eedfType = 'boltzmann'; % 'boltzmann' or 'prescribedEedf'
            gui.Setup.electronKinetics.shapeParameter = 1; % Only used for prescribedEedf
            gui.Setup.electronKinetics.ionizationOperatorType = 'usingSDCS'; % 'conservative', 'oneTakesAll', 'usingSDCS'
            gui.Setup.electronKinetics.growthModelType = 'temporal'; % 'temporal' or 'spatial'
            gui.Setup.electronKinetics.includeEECollisions = false;
            gui.Setup.electronKinetics.LXCatFiles = {'Nitrogen/N2_LXCat.txt', 'Nitrogen/N2_rot_LXCat.txt'}; % Cell array of strings
            gui.Setup.electronKinetics.LXCatExtraFiles = {}; % Optional
            gui.Setup.electronKinetics.effectiveCrossSectionPopulations = {}; % Optional
            gui.Setup.electronKinetics.CARgases = {}; % Optional
            
            % Gas Properties [cite: 5] - Inside electronKinetics
            gui.Setup.electronKinetics.gasProperties.mass = {'Databases/masses.txt'};
            gui.Setup.electronKinetics.gasProperties.fraction = {'N2 = 1'}; % Cell array of strings
            gui.Setup.electronKinetics.gasProperties.harmonicFrequency = {'Databases/harmonicFrequencies.txt'};
            gui.Setup.electronKinetics.gasProperties.anharmonicFrequency = {'Databases/anharmonicFrequencies.txt'};
            gui.Setup.electronKinetics.gasProperties.rotationalConstant = {'Databases/rotationalConstants.txt'};
            gui.Setup.electronKinetics.gasProperties.electricQuadrupoleMoment = {'Databases/quadrupoleMoment.txt'};
            gui.Setup.electronKinetics.gasProperties.OPBParameter = {'Databases/OPBParameter.txt'};
            
            % State Properties [cite: 5]
            gui.Setup.electronKinetics.stateProperties.energy = {'N2(X,v=*) = harmonicOscillatorEnergy', 'N2(X,v=0,J=*) = rigidRotorEnergy'};
            gui.Setup.electronKinetics.stateProperties.statisticalWeight = {'N2(X,v=*) = 1.0', 'N2(X,v=0,J=*) = rotationalDegeneracy_N2'};
            gui.Setup.electronKinetics.stateProperties.population = {'N2(X) = 1.0', 'Nitrogen/N2_vibpop.txt', 'N2(X,v=0,J=*) = boltzmannPopulation@gasTemperature'};
            
            % Numerics
            gui.Setup.electronKinetics.numerics.energyGrid.maxEnergy = 1; % (use 18-20 for time-dependent simulations)
            gui.Setup.electronKinetics.numerics.energyGrid.cellNumber = 1000; % (use 1800-2000 for time-dependent simulations)
            gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.isOn = false;
            gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay = 20;
            gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay = 25;
            gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.updateFactor = 0.05;
            gui.Setup.electronKinetics.numerics.maxPowerBalanceRelError = 1e-9;
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.algorithm = 'mixingDirectSolutions';
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.mixingParameter = 0.7;
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError = 1e-9;
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol = 1e-300;
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol = 1e-6;
            gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep = 1e-7;

            % GUI [cite: 13, 14]
            gui.Setup.gui.isOn = true;
            gui.Setup.gui.refreshFrequency = 1;

            % Output [cite: 14]
            gui.Setup.output.isOn = false;
            gui.Setup.output.dataFormat = 'hdf5+txt'; % 'txt', 'hdf5', 'hdf5+txt'
            gui.Setup.output.folder = 'simulation_1';
            gui.Setup.output.dataSets = {'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTables'}; % Cell array
        end

        function createGUI(gui)
            % Create the main figure
            iconPath = fullfile('figs', 'icon.png'); % Relative path to icon in figs folder
            if ~isfile(iconPath)
                warning('Icon file not found: %s. Using default icon.', iconPath);
                iconPath = ''; % Use default if not found
            end

            figWidth = 1000;
            figHeight = 750;
            monitorPositions = get(groot, 'MonitorPositions');
            screenSize = monitorPositions(1, :);
            figX = screenSize(1) + (screenSize(3) - figWidth) / 2;
            figY = screenSize(2) + (screenSize(4) - figHeight) / 2;

            gui.Fig = uifigure('Name', 'LoKI-B', ...
                'Position', [figX, figY, figWidth, figHeight], ...
                'NumberTitle', 'off', ...
                'Resize', 'on', ... % Allow resizing
                'Icon', iconPath, ...
                'WindowStyle', 'normal', ... % Keep window visible
                'WindowState', 'maximized', ... % Start maximized to avoid hidden title bar on different screens
                'Visible', 'off'); % Prevent initial layout "snap" while building/populating UI

            % Allow SizeChangedFcn to execute (MATLAB suppresses it when AutoResizeChildren is on)
            gui.Fig.AutoResizeChildren = 'off';

            % Handle responsive UI adjustments (font autoscaling, etc.)
            gui.Fig.SizeChangedFcn = @(src, evt) gui.handleFigureSizeChanged();

            % Main grid layout: two columns side by side
            mainGrid = uigridlayout(gui.Fig, [2, 2]);
            mainGrid.RowHeight = {'1x', 'fit'}; % Content takes most space, buttons at bottom
            mainGrid.ColumnWidth = {'1x', '1.5x'}; % Left (Working Conditions) narrower than right tabs
            mainGrid.ColumnSpacing = 10;
            mainGrid.RowSpacing = 10;
            mainGrid.Padding = [10 10 10 10];

            % --- Left Panel: Working Conditions (always visible) ---
            leftPanel = uipanel(mainGrid, 'BorderType', 'none');
            leftPanel.Layout.Row = 1;
            leftPanel.Layout.Column = 1;

            % Tab-style container inside left panel for visual homogeneity
            leftGrid = uigridlayout(leftPanel, [1, 1]);
            leftGrid.RowHeight = {'1x'};
            leftGrid.ColumnWidth = {'1x'};
            leftGrid.Padding = [0 0 0 0];

            leftTabGroup = uitabgroup(leftGrid);
            leftWorkingTab = uitab(leftTabGroup, 'Title', 'Working Conditions', 'Scrollable', 'on');

            gui.createWorkingConditionsPanel(leftWorkingTab);

            % --- Right Panel: Other tabs ---
            rightPanel = uipanel(mainGrid, 'BorderType', 'none');
            rightPanel.Layout.Row = 1;
            rightPanel.Layout.Column = 2;

            % Grid layout inside right panel to fill space
            rightGrid = uigridlayout(rightPanel, [1, 1]);
            rightGrid.RowHeight = {'1x'};
            rightGrid.ColumnWidth = {'1x'};
            rightGrid.Padding = [0 0 0 0];

            % Create tabs in right panel
            tabGroup = uitabgroup(rightGrid);
            gui.UIControls.tabs.tabGroup = tabGroup;

            % --- Electron Kinetics Tab ---
            kineticsTab = uitab(tabGroup, 'Title', 'Electron Kinetics', 'Scrollable', 'on');
            gui.UIControls.tabs.kineticsTab = kineticsTab;
            gui.createElectronKineticsPanel(kineticsTab);

            % --- Gas Properties Tab ---
            gasPropsTab = uitab(tabGroup, 'Title', 'Gas Properties', 'Scrollable', 'on');
            gui.createGasPropertiesPanel(gasPropsTab);

            % --- State Properties Tab ---
            statePropsTab = uitab(tabGroup, 'Title', 'State Properties', 'Scrollable', 'on');
            gui.createStatePropertiesPanel(statePropsTab);

            % --- Numerics Tab ---
            numericsTab = uitab(tabGroup, 'Title', 'Numerics', 'Scrollable', 'on');
            gui.createNumericsPanel(numericsTab);

            % --- Output Tab ---
            outputTab = uitab(tabGroup, 'Title', 'Output', 'Scrollable', 'on');
            gui.createOutputPanel(outputTab);

            % Enforce required LXCat files when leaving Electron Kinetics tab
            tabGroup.SelectionChangedFcn = @(src, evt) gui.handleTabChange(evt);

            % --- Button Panel (spans both columns) ---
            buttonPanel = uipanel(mainGrid, 'BorderType', 'none');
            buttonPanel.Layout.Row = 2;
            buttonPanel.Layout.Column = [1, 2]; % Span both columns
            buttonGrid = uigridlayout(buttonPanel, [1, 7]);
            % Add some padding/spacing if needed
            buttonGrid.ColumnWidth = {'fit', 'fit', '1x', 'fit', '1x', 'fit', 'fit'}; % Left buttons, centered help, right buttons
            buttonGrid.Padding = [10 10 10 10];
            buttonGrid.ColumnSpacing = 10;

            % Clear Button
            clearBtn = uibutton(buttonGrid, 'Text', 'Clear Setup', ...
                'ButtonPushedFcn', @gui.clearAllFields);
            clearBtn.Layout.Column = 1;

            % Load Button
            loadBtn = uibutton(buttonGrid, 'Text', 'Load Settings', ...
                'ButtonPushedFcn', @gui.loadSettings);
            loadBtn.Layout.Column = 2;

            % Help Button
            helpIcon = fullfile(matlabroot, 'toolbox', 'matlab', 'datamanager', ...
                '+datamanager', '+basicfit', '+icons', 'help.png');
            if ~isfile(helpIcon)
                helpIcon = '';
            end
            helpBtn = uibutton(buttonGrid, 'Text', 'Help', ...
                'Icon', helpIcon, ...
                'IconAlignment', 'left', ...
                'ButtonPushedFcn', @gui.openSetupFileHelp);
            helpBtn.Layout.Column = 4;

            % Save Button
            saveBtn = uibutton(buttonGrid, 'Text', 'Save Input File', ...
                'ButtonPushedFcn', @gui.saveInputFile);
            saveBtn.Layout.Column = 6;
            gui.UIControls.saveInputButton = saveBtn;

            % Run Button
            runBtn = uibutton(buttonGrid, 'Text', 'Generate & Run', ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @gui.runSimulation);
            runBtn.Layout.Column = 7;
            runBtn.BackgroundColor = [0.18 0.70 0.25];
            runBtn.FontColor = [1 1 1];
            gui.UIControls.runButton = runBtn;
        end

        function createWorkingConditionsPanel(gui, parent)
            % Use grid layout for better alignment and resizing
            grid = uigridlayout(parent, [16, 2]); % 15 data rows + 1 info section row
            % Fixed label column prevents the input column from collapsing to ~0 in narrow windows
            % (which was causing fields like totalSccmOutFlow to "disappear").
            grid.ColumnWidth = {220, '1x'}; % Label/checkbox, Edit field
            grid.RowHeight = [{42}, {42}, repmat({22}, 1, 13), {'1x'}]; % Required-condition rows + 13 fixed rows + flexible info row
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 5;
            grid.ColumnSpacing = 10;

            row = 1;
            % Reduced Field [cite: 1]
            reducedFieldLabel = uilabel(grid, 'Text', 'Reduced Field (Td):');
            reducedFieldLabel.Layout.Row = row;
            reducedFieldLabel.Layout.Column = 1;
            reducedFieldGrid = uigridlayout(grid, [2, 1]);
            reducedFieldGrid.RowHeight = {22, 14};
            reducedFieldGrid.ColumnWidth = {'1x'};
            reducedFieldGrid.Padding = [0 0 0 0];
            reducedFieldGrid.RowSpacing = 2;
            reducedFieldGrid.Layout.Row = row;
            reducedFieldGrid.Layout.Column = 2;

            gui.UIControls.workingConditions.reducedFieldBorderPanel = uipanel(reducedFieldGrid, ...
                'BorderType', 'none', ...
                'BackgroundColor', [0.94 0.94 0.94]);
            gui.UIControls.workingConditions.reducedFieldBorderPanel.Layout.Row = 1;
            gui.UIControls.workingConditions.reducedFieldBorderPanel.Layout.Column = 1;

            reducedFieldBorderGrid = uigridlayout(gui.UIControls.workingConditions.reducedFieldBorderPanel, [1, 1]);
            reducedFieldBorderGrid.RowHeight = {'1x'};
            reducedFieldBorderGrid.ColumnWidth = {'1x'};
            reducedFieldBorderGrid.Padding = [1 1 1 1];
            reducedFieldBorderGrid.RowSpacing = 0;
            reducedFieldBorderGrid.ColumnSpacing = 0;

            gui.UIControls.workingConditions.reducedField = uieditfield(reducedFieldBorderGrid, 'text', ...
                'Tooltip', 'e.g., 10 or logspace(1,2,10)', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.reducedField', evt.Value));
            gui.UIControls.workingConditions.reducedField.Layout.Row = 1;
            gui.UIControls.workingConditions.reducedField.Layout.Column = 1;

            gui.UIControls.workingConditions.reducedFieldMessage = uilabel(reducedFieldGrid, ...
                'Text', 'Reduced Field is mandatory for EEDF of type Boltzmann', ...
                'FontColor', [0.75 0.1 0.1], ...
                'FontSize', 10, ...
                'Visible', 'off');
            gui.UIControls.workingConditions.reducedFieldMessage.Layout.Row = 2;
            gui.UIControls.workingConditions.reducedFieldMessage.Layout.Column = 1;

            row = row + 1;
            % Electron Temperature [cite: 1]
            electronTempLabel = uilabel(grid, 'Text', 'Electron Temperature (eV):');
            electronTempLabel.Layout.Row = row;
            electronTempLabel.Layout.Column = 1;
            electronTemperatureGrid = uigridlayout(grid, [2, 1]);
            electronTemperatureGrid.RowHeight = {22, 14};
            electronTemperatureGrid.ColumnWidth = {'1x'};
            electronTemperatureGrid.Padding = [0 0 0 0];
            electronTemperatureGrid.RowSpacing = 2;
            electronTemperatureGrid.Layout.Row = row;
            electronTemperatureGrid.Layout.Column = 2;

            gui.UIControls.workingConditions.electronTemperatureBorderPanel = uipanel(electronTemperatureGrid, ...
                'BorderType', 'none', ...
                'BackgroundColor', [0.94 0.94 0.94]);
            gui.UIControls.workingConditions.electronTemperatureBorderPanel.Layout.Row = 1;
            gui.UIControls.workingConditions.electronTemperatureBorderPanel.Layout.Column = 1;

            electronTemperatureBorderGrid = uigridlayout(gui.UIControls.workingConditions.electronTemperatureBorderPanel, [1, 1]);
            electronTemperatureBorderGrid.RowHeight = {'1x'};
            electronTemperatureBorderGrid.ColumnWidth = {'1x'};
            electronTemperatureBorderGrid.Padding = [1 1 1 1];
            electronTemperatureBorderGrid.RowSpacing = 0;
            electronTemperatureBorderGrid.ColumnSpacing = 0;

            gui.UIControls.workingConditions.electronTemperature = uieditfield(electronTemperatureBorderGrid, 'text', ...
                'Tooltip', 'e.g., 1.5 or linspace(0.1, 5, 20)', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.electronTemperature', evt.Value));
            gui.UIControls.workingConditions.electronTemperature.Layout.Row = 1;
            gui.UIControls.workingConditions.electronTemperature.Layout.Column = 1;

            gui.UIControls.workingConditions.electronTemperatureMessage = uilabel(electronTemperatureGrid, ...
                'Text', 'Electron Temperature is mandatory for EEDF of type PrescribedEEDF', ...
                'FontColor', [0.75 0.1 0.1], ...
                'FontSize', 10, ...
                'Visible', 'off');
            gui.UIControls.workingConditions.electronTemperatureMessage.Layout.Row = 2;
            gui.UIControls.workingConditions.electronTemperatureMessage.Layout.Column = 1;

            row = row + 1;
            % Excitation Frequency [cite: 1]
            excitationFreqLabel = uilabel(grid, 'Text', 'Excitation Frequency (Hz):');
            excitationFreqLabel.Layout.Row = row;
            excitationFreqLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.excitationFrequency = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.handleExcitationFrequencyChange(src, evt.Value));
            gui.UIControls.workingConditions.excitationFrequency.Layout.Row = row;
            gui.UIControls.workingConditions.excitationFrequency.Layout.Column = 2;

            row = row + 1;
            % Gas Pressure [cite: 1]
            gasPressureCheckbox = uicheckbox(grid, 'Text', 'Gas Pressure (Pa):', ...
                'ValueChangedFcn', @(src, evt) gui.toggleGasPressureEnable(evt.Value));
            gasPressureCheckbox.Layout.Row = row;
            gasPressureCheckbox.Layout.Column = 1;
            gui.UIControls.workingConditions.gasPressureCheckbox = gasPressureCheckbox;
            
            gui.UIControls.workingConditions.gasPressure = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.gasPressure', evt.Value));
            gui.UIControls.workingConditions.gasPressure.Layout.Row = row;
            gui.UIControls.workingConditions.gasPressure.Layout.Column = 2;

            row = row + 1;
            % Gas Temperature [cite: 1]
            gasTempLabel = uilabel(grid, 'Text', 'Gas Temperature (K):');
            gasTempLabel.Layout.Row = row;
            gasTempLabel.Layout.Column = 1;
            
            gui.UIControls.workingConditions.gasTemperature = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.gasTemperature', evt.Value));
            gui.UIControls.workingConditions.gasTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.gasTemperature.Layout.Column = 2;

            row = row + 1;
            % Wall Temperature [cite: 2]
            wallTempLabel = uilabel(grid, 'Text', 'Wall Temperature (K):');
            wallTempLabel.Layout.Row = row;
            wallTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.wallTemperature = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.wallTemperature', evt.Value));
            gui.UIControls.workingConditions.wallTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.wallTemperature.Layout.Column = 2;

            row = row + 1;
            % External Temperature [cite: 2]
            extTempLabel = uilabel(grid, 'Text', 'External Temperature (K):');
            extTempLabel.Layout.Row = row;
            extTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.extTemperature = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.extTemperature', evt.Value));
            gui.UIControls.workingConditions.extTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.extTemperature.Layout.Column = 2;

            row = row + 1;
            % Surface Site Density [cite: 2]
            surfaceSiteDensityLabel = uilabel(grid, 'Text', 'Surface Site Density (m^-2):');
            surfaceSiteDensityLabel.Layout.Row = row;
            surfaceSiteDensityLabel.Layout.Column = 1;
            
            gui.UIControls.workingConditions.surfaceSiteDensity = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.surfaceSiteDensity', evt.Value));
            gui.UIControls.workingConditions.surfaceSiteDensity.Layout.Row = row;
            gui.UIControls.workingConditions.surfaceSiteDensity.Layout.Column = 2;

            row = row + 1;
            % Electron Density [cite: 3]
            electronDensityCheckbox = uicheckbox(grid, 'Text', 'Electron Density (m^-3):', ...
                'ValueChangedFcn', @(src, evt) gui.toggleElectronDensityEnable(evt.Value));
            electronDensityCheckbox.Layout.Row = row;
            electronDensityCheckbox.Layout.Column = 1;
            gui.UIControls.workingConditions.electronDensityCheckbox = electronDensityCheckbox;
            
            gui.UIControls.workingConditions.electronDensity = uieditfield(grid, 'text', ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.electronDensity', evt.Value));
            gui.UIControls.workingConditions.electronDensity.Layout.Row = row;
            gui.UIControls.workingConditions.electronDensity.Layout.Column = 2;

            row = row + 1;
            % Chamber Length [cite: 3]
            chamberLengthLabel = uilabel(grid, 'Text', 'Chamber Length (m):');
            chamberLengthLabel.Layout.Row = row;
            chamberLengthLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.chamberLength = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.chamberLength', evt.Value));
            gui.UIControls.workingConditions.chamberLength.Layout.Row = row;
            gui.UIControls.workingConditions.chamberLength.Layout.Column = 2;

            row = row + 1;
            % Chamber Radius [cite: 3]
            chamberRadiusLabel = uilabel(grid, 'Text', 'Chamber Radius (m):');
            chamberRadiusLabel.Layout.Row = row;
            chamberRadiusLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.chamberRadius = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.chamberRadius', evt.Value));
            gui.UIControls.workingConditions.chamberRadius.Layout.Row = row;
            gui.UIControls.workingConditions.chamberRadius.Layout.Column = 2;

            row = row + 1;
            % Total SCCM Inflow [cite: 3]
            totalSccmInFlowLabel = uilabel(grid, 'Text', 'Total SCCM Inflow (sccm):');
            totalSccmInFlowLabel.Layout.Row = row;
            totalSccmInFlowLabel.Layout.Column = 1;
            
            gui.UIControls.workingConditions.totalSccmInFlow = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.totalSccmInFlow', evt.Value));
            gui.UIControls.workingConditions.totalSccmInFlow.Layout.Row = row;
            gui.UIControls.workingConditions.totalSccmInFlow.Layout.Column = 2;

            row = row + 1;
            % Total SCCM Outflow [cite: 3]
            totalSccmOutFlowLabel = uilabel(grid, 'Text', 'Total SCCM Outflow (sccm):');
            totalSccmOutFlowLabel.Layout.Row = row;
            totalSccmOutFlowLabel.Layout.Column = 1;
            % One-row overlay behavior:
            % - dropdown always visible
            % - when 'Number' is selected, show a numeric edit field on top,
            %   slightly narrower so the dropdown arrow remains clickable.
            outFlowPanel = uipanel(grid, 'BorderType', 'none');
            outFlowPanel.Layout.Row = row;
            outFlowPanel.Layout.Column = 2;
            gui.UIControls.workingConditions.totalSccmOutFlowPanel = outFlowPanel;
            % Allow SizeChangedFcn to execute (MATLAB suppresses it when AutoResizeChildren is on)
            outFlowPanel.AutoResizeChildren = 'off';
            % Re-layout overlay whenever this cell is resized (more reliable than figure SizeChanged alone)
            outFlowPanel.SizeChangedFcn = @(src, evt) gui.layoutTotalSccmOutFlowOverlay();

            dd = uidropdown(outFlowPanel, ...
                'Items', {'ensureIsobaric', 'totalSccmInFlow', 'Number'}, ...
                'Value', 'ensureIsobaric', ...
                'Enable', 'off');
            gui.UIControls.workingConditions.totalSccmOutFlowType = dd;

            ef = uieditfield(outFlowPanel, 'numeric', ...
                'Limits', [0, Inf], ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.totalSccmOutFlow', evt.Value), ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left');
            ef.Visible = 'off'; % start hidden behind dropdown
            gui.UIControls.workingConditions.totalSccmOutFlow = ef;

            % Initial positioning (will be updated on resize)
            gui.layoutTotalSccmOutFlowOverlay();

            row = row + 1;

            % % Discharge Current [cite: 3] - Optional
            dischargeCurrentLabel = uilabel(grid, 'Text', 'Discharge Current (A):');
            dischargeCurrentLabel.Layout.Row = row;
            dischargeCurrentLabel.Layout.Column = 1;

            gui.UIControls.workingConditions.dischargeCurrent = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.dischargeCurrent', evt.Value));
            gui.UIControls.workingConditions.dischargeCurrent.Layout.Row = row;
            gui.UIControls.workingConditions.dischargeCurrent.Layout.Column = 2;

            row = row + 1;
            % Discharge Power Density [cite: 3] - Optional
            dischargePowerLabel = uilabel(grid, 'Text', 'Discharge Power Density (W/m³):');
            dischargePowerLabel.Layout.Row = row;
            dischargePowerLabel.Layout.Column = 1;

            gui.UIControls.workingConditions.dischargePowerDensity = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'HorizontalAlignment', 'left', ...
                'Enable', 'off', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.dischargePowerDensity', evt.Value));
            gui.UIControls.workingConditions.dischargePowerDensity.Layout.Row = row;
            gui.UIControls.workingConditions.dischargePowerDensity.Layout.Column = 2;

            % --- Information Section (row 16, spans both columns) ---
            row = row + 1;
            
            % Info section: image aligned left, text takes remaining horizontal space
            infoContainer = uigridlayout(grid, [1, 2]);
            infoContainer.Layout.Row = row;
            infoContainer.Layout.Column = [1, 2];
            % Column widths are made responsive in handleFigureSizeChanged()
            infoContainer.ColumnWidth = {150, '1x'}; % {Image, Text}
            infoContainer.RowHeight = {'1x'};
            % Slight top padding so the image doesn't stick to the top edge
            infoContainer.Padding = [0 0 0 20];
            infoContainer.ColumnSpacing = 15;
            gui.UIControls.workingConditions.infoContainer = infoContainer;
            
            % LoKI image
            lokiImagePath = fullfile('figs', 'LoKI.png');
            if isfile(lokiImagePath)
                lokiImage = uiimage(infoContainer, 'ImageSource', lokiImagePath);
                lokiImage.Layout.Row = 1;
                lokiImage.Layout.Column = 1;
                % Expand to fill available space without distorting
                lokiImage.ScaleMethod = 'fit';
                gui.UIControls.workingConditions.infoImage = lokiImage;
            else
                % Placeholder if image not found
                lokiPlaceholder = uilabel(infoContainer, 'Text', '[LoKI]', 'FontSize', 16, 'FontWeight', 'bold');
                lokiPlaceholder.Layout.Row = 1;
                lokiPlaceholder.Layout.Column = 1;
            end

            % Description text
            textContainer = uigridlayout(infoContainer, [2, 1]);
            textContainer.Layout.Row = 1;
            textContainer.Layout.Column = 2;
            textContainer.RowHeight = {'fit', '1x'};
            textContainer.Padding = [0 0 0 0];
            textContainer.RowSpacing = 5;

            titleText = 'Input Graphic User Interface for LoKI-B';
            titleLabel = uilabel(textContainer, 'Text', titleText, ...
                'FontSize', 16, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left');
            titleLabel.Layout.Row = 1;
            titleLabel.Layout.Column = 1;

            descText = sprintf(['\n'...
                                'Developed by group N-PRiME of IPFN/IST (Lisbon, Portugal)\n' ...
                'Build hash: 0x4D522D41525047']);
            descLabel = uilabel(textContainer, 'Text', descText, ...
                'FontSize', 15, ...
                'WordWrap', 'on', ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'left');
            descLabel.Layout.Row = 2;
            descLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.infoDescLabel = descLabel;
            gui.UIControls.workingConditions.infoTitleLabel = titleLabel;
        end

        function createElectronKineticsPanel(gui, parent)
            % Main grid for electron kinetics panel
            grid = uigridlayout(parent, [3, 1]);
            grid.RowHeight = {'fit', '1x', 'fit'}; 
            grid.Padding = [5 5 5 5];
            grid.RowSpacing = 10;

            % --- General Settings Panel ---
            generalPanel = uipanel(grid, 'Title', 'General Kinetics Settings');
            generalPanel.Layout.Row = 1;
            generalGrid = uigridlayout(generalPanel, [6, 2]); % Rows, Columns
            generalGrid.ColumnWidth = {'fit', '1x'};
            generalGrid.RowHeight = repmat({'fit'}, 1, 6);
            generalGrid.Padding = [10 10 10 10];
            generalGrid.RowSpacing = 5;

            row = 1;
            % Is On [cite: 4]
            gui.UIControls.electronKinetics.isOn = uicheckbox(generalGrid, ...
                'Text', 'Enable Electron Kinetics', ...
                'Value', true, ...
                'ValueChangedFcn', @(src, evt) gui.toggleElectronKineticsEnable(evt.Value, true));
            gui.UIControls.electronKinetics.isOn.Layout.Row = row;
            gui.UIControls.electronKinetics.isOn.Layout.Column = [1, 2]; % Span columns

            row = row + 1;
            % EEDF Type [cite: 4]
            eedfTypeLabel = uilabel(generalGrid, 'Text', 'EEDF Type:');
            eedfTypeLabel.Layout.Row = row;
            eedfTypeLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.eedfType = uidropdown(generalGrid, 'Items', {'boltzmann', 'prescribedEedf'}, 'ValueChangedFcn', @(src, evt) gui.handleEedfTypeChange(evt.Value));
            gui.UIControls.electronKinetics.eedfType.Layout.Row = row;
            gui.UIControls.electronKinetics.eedfType.Layout.Column = 2;

            row = row + 1;
            % Shape Parameter [cite: 4] - Only visible for prescribedEedf 
            gui.UIControls.electronKinetics.shapeParameterLabel = uilabel(generalGrid, 'Text', 'Shape Parameter:');
            shapeParamGrid = uigridlayout(generalGrid, [1, 2]);
            shapeParamGrid.ColumnWidth = {'0.6x', '0.4x'}; % Give more space to text field
            shapeParamGrid.Padding = [0 0 0 0];
            shapeParamGrid.ColumnSpacing = 5;
            
            gui.UIControls.electronKinetics.shapeParameter = uislider(shapeParamGrid, ...
                'Limits', [1, 2], 'Value', 1, 'MajorTicks', [1, 2], 'MajorTickLabels', {'1 (Maxwellian)', '2 (Druyvesteyn)'}, 'ValueChangedFcn', ...
                @(src, evt) gui.sliderShapeParameterChanged(src, evt.Value), 'Visible', 'off');
            gui.UIControls.electronKinetics.shapeParameter.Layout.Row = 1;
            gui.UIControls.electronKinetics.shapeParameter.Layout.Column = 1;

            gui.UIControls.electronKinetics.shapeParameterPrecise = uieditfield(shapeParamGrid, 'numeric', 'Limits', ...
                 [1, 2], 'Value', 1, 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.textShapeParameterChanged(src, evt.Value), ...
                 'Visible', 'off');
            gui.UIControls.electronKinetics.shapeParameterPrecise.Layout.Row = 1;
            gui.UIControls.electronKinetics.shapeParameterPrecise.Layout.Column = 2;
            
            shapeParamGrid.Layout.Row = row;
            shapeParamGrid.Layout.Column = 2;

            row = row + 1;
            % Ionization Operator Type [cite: 5]
            ionizationOperatorLabel = uilabel(generalGrid, 'Text', 'Ionization Operator:');
            ionizationOperatorLabel.Layout.Row = row;
            ionizationOperatorLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.ionizationOperatorType = uidropdown(generalGrid, ...
                'Items', {'conservative', 'oneTakesAll', 'equalSharing', 'usingSDCS'}, ...
                'ValueChangedFcn', @(src, evt) gui.handleIonizationOperatorChange(evt.Value));
            gui.UIControls.electronKinetics.ionizationOperatorType.Layout.Row = row;
            gui.UIControls.electronKinetics.ionizationOperatorType.Layout.Column = 2;

            row = row + 1;
            % Growth Model Type [cite: 5]
            growthModelLabel = uilabel(generalGrid, 'Text', 'Growth Model:');
            growthModelLabel.Layout.Row = row;
            growthModelLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.growthModelType = uidropdown(generalGrid, 'Items', {'temporal', 'spatial'}, 'ValueChangedFcn', ...
                        @(src, evt) gui.updateField(src, 'electronKinetics.growthModelType', evt.Value));
            gui.UIControls.electronKinetics.growthModelType.Layout.Row = row;
            gui.UIControls.electronKinetics.growthModelType.Layout.Column = 2;

            row = row + 1;
            % Include e-e Collisions [cite: 5]
            gui.UIControls.electronKinetics.includeEECollisions = uicheckbox(generalGrid, ...
                'Text', 'Include e-e Collisions', ...
                'ValueChangedFcn', @(src, evt) gui.handleEECollisionsChange(evt.Value));
            gui.UIControls.electronKinetics.includeEECollisions.Layout.Row = row;
            gui.UIControls.electronKinetics.includeEECollisions.Layout.Column = [1, 2];

            % --- Files Panel ---
            filesPanel = uipanel(grid, 'Title', 'Cross Sections');
            filesPanel.Layout.Row = 2;
            filesGrid = uigridlayout(filesPanel, [3, 4]); % Use consistent columns
            filesGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            filesGrid.RowHeight = {'1x', 30, 30}; % LXCat Files grows with window, others fixed ~1 item
            filesGrid.Padding = [10 10 10 10];
            filesGrid.RowSpacing = 10;
            filesGrid.ColumnSpacing = 5;

            % LXCat Files [cite: 5]
            lxcatFilesLabel = uilabel(filesGrid, 'Text', 'LXCat Files:');
            lxcatFilesLabel.Layout.Row = 1;
            lxcatFilesLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.LXCatFiles = uilistbox(filesGrid, 'Multiselect', 'on');
            gui.UIControls.electronKinetics.LXCatFiles.Layout.Row = 1;
            gui.UIControls.electronKinetics.LXCatFiles.Layout.Column = 2;
            lxcatFilesAddButton = uibutton(filesGrid, 'Text', 'Add', 'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.LXCatFiles', true));
            lxcatFilesAddButton.Layout.Row = 1;
            lxcatFilesAddButton.Layout.Column = 3;
            lxcatFilesRemoveButton = uibutton(filesGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.LXCatFiles'));
            lxcatFilesRemoveButton.Layout.Row = 1;
            lxcatFilesRemoveButton.Layout.Column = 4;

            % LXCat Extra Files [cite: 5] - Optional
            extraCheckbox = uicheckbox(filesGrid, 'Text', 'LXCat Extra Files:', ...
                'ValueChangedFcn', @(src, evt) gui.toggleLXCatExtraEnable(evt.Value));
            extraCheckbox.Layout.Row = 2;
            extraCheckbox.Layout.Column = 1;
            gui.UIControls.electronKinetics.LXCatExtraCheckbox = extraCheckbox;
            
            gui.UIControls.electronKinetics.LXCatExtraFiles = uilistbox(filesGrid, 'Multiselect', 'on', 'Enable', 'off');
            gui.UIControls.electronKinetics.LXCatExtraFiles.Layout.Row = 2;
            gui.UIControls.electronKinetics.LXCatExtraFiles.Layout.Column = 2;
            extraAddButton = uibutton(filesGrid, 'Text', 'Add', ...
                'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.LXCatExtraFiles', true), ...
                'Enable', 'off');
            extraAddButton.Layout.Row = 2;
            extraAddButton.Layout.Column = 3;
            extraRemoveButton = uibutton(filesGrid, 'Text', 'Remove', ...
                'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.LXCatExtraFiles'), ...
                'Enable', 'off');
            extraRemoveButton.Layout.Row = 2;
            extraRemoveButton.Layout.Column = 4;

            % Effective Cross Section Populations [cite: 5] - Optional
            effPopCheckbox = uicheckbox(filesGrid, 'Text', 'Effective Cross Section Populations:', ...
                'ValueChangedFcn', @(src, evt) gui.toggleEffectivePopEnable(evt.Value));
            effPopCheckbox.Layout.Row = 3;
            effPopCheckbox.Layout.Column = 1;
            gui.UIControls.electronKinetics.effectivePopCheckbox = effPopCheckbox;
            
            gui.UIControls.electronKinetics.effectiveCrossSectionPopulations = uilistbox(filesGrid, 'Multiselect', 'on', 'Enable', 'off');
            gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Layout.Row = 3;
            gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Layout.Column = 2;
            effPopAddButton = uibutton(filesGrid, 'Text', 'Add', ...
                'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.effectiveCrossSectionPopulations', true), ...
                'Enable', 'off');
            effPopAddButton.Layout.Row = 3;
            effPopAddButton.Layout.Column = 3;
            effPopRemoveButton = uibutton(filesGrid, 'Text', 'Remove', ...
                'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.effectiveCrossSectionPopulations'), ...
                'Enable', 'off');
            effPopRemoveButton.Layout.Row = 3;
            effPopRemoveButton.Layout.Column = 4;

            % --- CAR Gases Panel ---
            carPanel = uipanel(grid, 'Title', 'CAR Gases');
            carPanel.Layout.Row = 3;
            carGrid = uigridlayout(carPanel, [1, 4]);
            carGrid.ColumnWidth = {'fit', '1x', 'fit', 'fit'};
            carGrid.RowHeight = {80}; % Fixed height
            carGrid.Padding = [10 10 10 10];
            carGrid.RowSpacing = 5;
            carGrid.ColumnSpacing = 5;

            carCheckbox = uicheckbox(carGrid, 'Text', 'CAR Gases:', ...
                'ValueChangedFcn', @(src, evt) gui.toggleCARGasEnable(evt.Value));
            carCheckbox.Layout.Row = 1;
            carCheckbox.Layout.Column = 1;
            gui.UIControls.electronKinetics.CARcheckbox = carCheckbox;
            
            gui.UIControls.electronKinetics.CARgases = uilistbox(carGrid, ...
                'Multiselect', 'on', ...
                'Enable', 'off', ...
                'DoubleClickedFcn', @(src, evt) gui.editListItem(src, 'electronKinetics.CARgases'));
            gui.UIControls.electronKinetics.CARgases.Layout.Row = 1;
            gui.UIControls.electronKinetics.CARgases.Layout.Column = 2;
            carAddButton = uibutton(carGrid, 'Text', 'Add', ...
                'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.CARgases', false), ...
                'Enable', 'off');
            carAddButton.Layout.Row = 1;
            carAddButton.Layout.Column = 3;
            gui.UIControls.electronKinetics.CARaddButton = carAddButton;
            carRemoveButton = uibutton(carGrid, 'Text', 'Remove', ...
                'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.CARgases'), ...
                'Enable', 'off');
            carRemoveButton.Layout.Row = 1;
            carRemoveButton.Layout.Column = 4;
            gui.UIControls.electronKinetics.CARremoveButton = carRemoveButton;
        end

        function handleFigureSizeChanged(gui)
            % Responsively adjust UI for small windows:
            % - Autoscale the Working Conditions info text font size so it stays visible
            if isempty(gui) || ~isvalid(gui) || isempty(gui.UIControls)
                return;
            end

            % Responsive layout for the info block: when the window gets narrow,
            % shrink the image column so the text gets more horizontal space.
            if isfield(gui.UIControls, 'workingConditions') && isfield(gui.UIControls.workingConditions, 'infoContainer')
                c = gui.UIControls.workingConditions.infoContainer;
                if ~isempty(c) && isvalid(c)
                    w = gui.Fig.Position(3); % figure width
                    if w < 900
                        imgW = 110;
                    elseif w < 1200
                        imgW = 130;
                    else
                        imgW = 150;
                    end
                    c.ColumnWidth = {imgW, '1x'};
                end
            end

            if isfield(gui.UIControls, 'workingConditions') && isfield(gui.UIControls.workingConditions, 'infoDescLabel')
                lbl = gui.UIControls.workingConditions.infoDescLabel;
                if ~isempty(lbl) && isvalid(lbl)
                    w = gui.Fig.Position(3); % figure width
                    % Simple piecewise scaling: clamp to keep readability while fitting small windows.
                    if w < 900
                        fs = 10;
                    elseif w < 1200
                        fs = 12;
                    else
                        fs = 14;
                    end
                    lbl.FontSize = fs;
                end
            end

            % Keep SCCM outflow overlay aligned with panel size
            if isfield(gui.UIControls, 'workingConditions') && isfield(gui.UIControls.workingConditions, 'totalSccmOutFlowPanel')
                gui.layoutTotalSccmOutFlowOverlay();
            end
        end

        function layoutTotalSccmOutFlowOverlay(gui)
            % Layout helper for the "overlay numeric field on top of dropdown" behavior
            if isempty(gui) || ~isvalid(gui) || isempty(gui.UIControls) || ~isfield(gui.UIControls, 'workingConditions')
                return;
            end
            wc = gui.UIControls.workingConditions;
            if ~isfield(wc, 'totalSccmOutFlowPanel') || ~isfield(wc, 'totalSccmOutFlowType') || ~isfield(wc, 'totalSccmOutFlow')
                return;
            end

            p = wc.totalSccmOutFlowPanel;
            dd = wc.totalSccmOutFlowType;
            ef = wc.totalSccmOutFlow;
            if isempty(p) || isempty(dd) || isempty(ef) || ~isvalid(p) || ~isvalid(dd) || ~isvalid(ef)
                return;
            end

            % Panel is grid-managed; use its current pixel size to position children.
            panelPos = p.Position;
            W = max(1, panelPos(3));
            H = max(1, panelPos(4));

            h = min(22, H);
            y = max(0, (H - h) / 2);

            dd.Position = [0, y, W, h];

            % Leave room for dropdown arrow (~32px) so user can still change mode
            arrowW = 34;
            efW = max(60, W - arrowW);
            ef.Position = [0, y, efW, h];
        end

        function createGasPropertiesPanel(gui, parent)
            % Main grid for gas properties panel
            grid = uigridlayout(parent, [2, 1]);
            grid.RowHeight = {'fit', '1x'};
            grid.Padding = [5 5 5 5];
            grid.RowSpacing = 10;

            % --- Gas Files Panel ---
            gasPropsPanel = uipanel(grid, 'Title', 'Gas Properties Files');
            gasPropsPanel.Layout.Row = 1;
            gasPropsGrid = uigridlayout(gasPropsPanel, [6, 3]);
            gasPropsGrid.ColumnWidth = {'fit', '1x', 'fit'};
            gasPropsGrid.RowHeight = repmat({50}, 1, 6);
            gasPropsGrid.Padding = [10 10 10 10];
            gasPropsGrid.RowSpacing = 5;
            gasPropsGrid.ColumnSpacing = 10;

            gasFileSpecs = {
                'mass', 'Mass Files:';
                'harmonicFrequency', 'Harmonic Frequency Files:';
                'anharmonicFrequency', 'Anharmonic Frequency Files:';
                'rotationalConstant', 'Rotational Constant Files:';
                'electricQuadrupoleMoment', 'Electric Quadrupole Moment Files:';
                'OPBParameter', 'OPB Parameter Files:'
            };

            for row = 1:size(gasFileSpecs, 1)
                fieldName = gasFileSpecs{row, 1};
                labelText = gasFileSpecs{row, 2};
                fieldPath = ['electronKinetics.gasProperties.', fieldName];

                fileLabel = uilabel(gasPropsGrid, 'Text', labelText);
                fileLabel.Layout.Row = row;
                fileLabel.Layout.Column = 1;

                gui.UIControls.gasProperties.(fieldName) = uilistbox(gasPropsGrid, ...
                    'Multiselect', 'on', ...
                    'Items', {}, ...
                    'DoubleClickedFcn', @(src, evt) gui.editListItem(src, fieldPath));
                gui.UIControls.gasProperties.(fieldName).Layout.Row = row;
                gui.UIControls.gasProperties.(fieldName).Layout.Column = 2;

                fileBtnGrid = uigridlayout(gasPropsGrid, [1, 2]);
                fileBtnGrid.Layout.Row = row;
                fileBtnGrid.Layout.Column = 3;
                fileBtnGrid.ColumnWidth = {'fit', 'fit'};
                fileBtnGrid.RowHeight = {'fit'};
                fileBtnGrid.Padding = [0 0 0 0];
                fileBtnGrid.ColumnSpacing = 5;

                addButton = uibutton(fileBtnGrid, 'Text', 'Browse', ...
                    'ButtonPushedFcn', @(src, evt) gui.addListItem(fieldPath, true));
                addButton.Layout.Row = 1;
                addButton.Layout.Column = 1;

                removeButton = uibutton(fileBtnGrid, 'Text', 'Remove', ...
                    'ButtonPushedFcn', @(src, evt) gui.removeListItem(fieldPath));
                removeButton.Layout.Row = 1;
                removeButton.Layout.Column = 2;
            end

            % --- Gas Fractions ---
            fractionsPanel = uipanel(grid, 'Title', 'Gas Fractions');
            fractionsPanel.Layout.Row = 2;
            fractionsGrid = uigridlayout(fractionsPanel, [2, 3]);
            fractionsGrid.ColumnWidth = {'fit', '1x', 'fit'};
            fractionsGrid.RowHeight = {'fit', '1x'};
            fractionsGrid.Padding = [10 10 10 10];
            fractionsGrid.RowSpacing = 5;
            fractionsGrid.ColumnSpacing = 5;

            % Fractions
            fractionsLabel = uilabel(fractionsGrid, 'Text', 'Gas Fractions:');
            fractionsLabel.Layout.Row = 1;
            fractionsLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.fraction = uilistbox(fractionsGrid, ...
                'Multiselect', 'on', ...
                'Items', {}, ...
                'DoubleClickedFcn', @(src, evt) gui.editListItem(src, 'electronKinetics.gasProperties.fraction'));
            gui.UIControls.gasProperties.fraction.Layout.Row = 2;
            gui.UIControls.gasProperties.fraction.Layout.Column = 2;
            fractionsBtnGrid = uigridlayout(fractionsGrid, [2, 1]);
            fractionsBtnGrid.Layout.Row = 2;
            fractionsBtnGrid.Layout.Column = 3;
            fractionsBtnGrid.RowHeight = {'fit', 'fit'};
            fractionsBtnGrid.Padding = [0 0 0 0];
            fractionsAddButton = uibutton(fractionsBtnGrid, 'Text', 'Add', 'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.gasProperties.fraction', false));
            fractionsRemoveButton = uibutton(fractionsBtnGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.gasProperties.fraction'));
        end

        function createStatePropertiesPanel(gui, parent)
            % Main grid for state properties panel
            grid = uigridlayout(parent, [3, 3]); % 3 rows for the three lists
            grid.ColumnWidth = {'fit', '1x', 'fit'};
            grid.RowHeight = {200, 200, 200}; % Increased height for each listbox
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 10;

            % Energy listbox
            row = 1;
            energyLabel = uilabel(grid, 'Text', 'Energy:');
            energyLabel.Layout.Row = row;
            energyLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.stateProperties.energy = uilistbox(grid, ...
                'Multiselect', 'on', ...
                'Items', {}, ...
                'DoubleClickedFcn', @(src, evt) gui.editListItem(src, 'electronKinetics.stateProperties.energy'));
            gui.UIControls.electronKinetics.stateProperties.energy.Layout.Row = row;
            gui.UIControls.electronKinetics.stateProperties.energy.Layout.Column = 2;
            energyBtnGrid = uigridlayout(grid, [2, 1]);
            energyBtnGrid.Layout.Row = row;
            energyBtnGrid.Layout.Column = 3;
            energyBtnGrid.RowHeight = {'fit', 'fit'};
            energyBtnGrid.Padding = [0 0 0 0];
            energyAddButton = uibutton(energyBtnGrid, 'Text', 'Add', 'ButtonPushedFcn', ...
                    @(src,evt) gui.addListItem('electronKinetics.stateProperties.energy', false));
            energyRemoveButton = uibutton(energyBtnGrid, 'Text', 'Remove', 'ButtonPushedFcn', ...
                    @(src,evt) gui.removeListItem('electronKinetics.stateProperties.energy'));

            % Statistical Weight listbox
            row = row + 1;
            statisticalWeightLabel = uilabel(grid, 'Text', 'Statistical Weight:');
            statisticalWeightLabel.Layout.Row = row;
            statisticalWeightLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.stateProperties.statisticalWeight = uilistbox(grid, ...
                'Multiselect', 'on', ...
                'Items', {}, ...
                'DoubleClickedFcn', @(src, evt) gui.editListItem(src, 'electronKinetics.stateProperties.statisticalWeight'));
            gui.UIControls.electronKinetics.stateProperties.statisticalWeight.Layout.Row = row;
            gui.UIControls.electronKinetics.stateProperties.statisticalWeight.Layout.Column = 2;
            statWeightBtnGrid = uigridlayout(grid, [2, 1]);
            statWeightBtnGrid.Layout.Row = row;
            statWeightBtnGrid.Layout.Column = 3;
            statWeightBtnGrid.RowHeight = {'fit', 'fit'};
            statWeightBtnGrid.Padding = [0 0 0 0];
            statWeightAddButton = uibutton(statWeightBtnGrid, 'Text', 'Add', 'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.stateProperties.statisticalWeight', false));
            statWeightRemoveButton = uibutton(statWeightBtnGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.stateProperties.statisticalWeight'));

            % Population listbox
            row = row + 1;
            statePopulationsLabel = uilabel(grid, 'Text', 'Population:');
            statePopulationsLabel.Layout.Row = row;
            statePopulationsLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.stateProperties.population = uilistbox(grid, ...
                'Multiselect', 'on', ...
                'Items', {}, ...
                'DoubleClickedFcn', @(src, evt) gui.editListItem(src, 'electronKinetics.stateProperties.population'));
            gui.UIControls.electronKinetics.stateProperties.population.Layout.Row = row;
            gui.UIControls.electronKinetics.stateProperties.population.Layout.Column = 2;
            btnGrid = uigridlayout(grid, [2, 1]);
            btnGrid.Layout.Row = row;
            btnGrid.Layout.Column = 3;
            btnGrid.RowHeight = {'fit', 'fit'};
            btnGrid.Padding = [0 0 0 0];
            statePopulationsAddButton = uibutton(btnGrid, 'Text', 'Add', 'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.stateProperties.population', false));
            statePopulationsAddButton.Layout.Row = 1;
            statePopulationsAddButton.Layout.Column = 1;
            statePopulationsRemoveButton = uibutton(btnGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.stateProperties.population'));
            statePopulationsRemoveButton.Layout.Row = 2;
            statePopulationsRemoveButton.Layout.Column = 1;

        end

        function createNumericsPanel(gui, parent)
            % Main grid for numerics panel
            grid = uigridlayout(parent, [2, 1]); % Energy Grid, Non-Linear Routines
            grid.RowHeight = {'fit', 'fit'};
            grid.Padding = [5 5 5 5];

            % --- Energy Grid Panel ---
            energyGridPanel = uipanel(grid, 'Title', 'Energy Grid Configuration');
            energyGridPanel.Layout.Row = 1;
            energyGridGrid = uigridlayout(energyGridPanel, [8, 2]); % Rows, Columns - increased for variable grid checkbox
            energyGridGrid.ColumnWidth = {'fit', '1x'};
            energyGridGrid.RowHeight = repmat({'fit'}, 1, 8);
            energyGridGrid.Padding = [10 10 10 10];
            energyGridGrid.RowSpacing = 5;

            row = 1;
            % % Variable Energy Grid (Coming Soon) - Mockup checkbox
            % variableGridCheckbox = uicheckbox(energyGridGrid, 'Text', 'Variable Energy Grid (Soon...)', 'Value', false, 'Enable', 'off');
            % variableGridCheckbox.Layout.Row = row;
            % variableGridCheckbox.Layout.Column = [1, 2];

            % row = row + 1;
            % Max Energy [cite: 10]
            maxEnergyLabel = uilabel(energyGridGrid, 'Text', 'Max Energy (eV):');
            maxEnergyLabel.Layout.Row = row;
            maxEnergyLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.maxEnergy', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy.Layout.Column = 2;

            row = row + 1;
            % Cell Number [cite: 10]
            cellNumberLabel = uilabel(energyGridGrid, 'Text', 'Energy Cell Number:');
            cellNumberLabel.Layout.Row = row;
            cellNumberLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.cellNumber = uieditfield(energyGridGrid, 'numeric', 'Limits', [1, Inf], 'ValueDisplayFormat', '%.0f', 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.cellNumber', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.cellNumber.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.cellNumber.Layout.Column = 2;

            % --- Smart Grid Panel ---
            row = row + 1;
            % Smart Grid [cite: 10]
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.isOn = uicheckbox(energyGridGrid, 'Text', 'Smart Grid', 'ValueChangedFcn', @(src, evt) gui.toggleSmartGridEnable(evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.isOn.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.isOn.Layout.Column = [1, 2];

            row = row + 1;
            % Min EEDF Decay
            minEedfDecayLabel = uilabel(energyGridGrid, 'Text', 'Min EEDF Decay:');
            minEedfDecayLabel.Layout.Row = row;
            minEedfDecayLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueDisplayFormat', '%.0f', 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay.Layout.Column = 2;

            row = row + 1;
            % Max EEDF Decay
            maxEedfDecayLabel = uilabel(energyGridGrid, 'Text', 'Max EEDF Decay:');
            maxEedfDecayLabel.Layout.Row = row;
            maxEedfDecayLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueDisplayFormat', '%.0f', 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay.Layout.Column = 2;

            row = row + 1;
            % Update Factor
            updateFactorLabel = uilabel(energyGridGrid, 'Text', 'Update Factor:');
            updateFactorLabel.Layout.Row = row;
            updateFactorLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor = uieditfield(energyGridGrid, ...
                'numeric', ...
                'Limits', [0.001, Inf], ...
                'Value', 0.05, ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.handleUpdateFactorChange(evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Layout.Column = 2;

            row = row + 1;
            % Max Power Balance Rel Error [cite: 12]
            maxPowerBalanceLabel = uilabel(energyGridGrid, 'Text', 'Max Power Bal. Rel. Error:');
            maxPowerBalanceLabel.Layout.Row = row;
            maxPowerBalanceLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.maxPowerBalanceRelError = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.maxPowerBalanceRelError', evt.Value));
            gui.UIControls.electronKinetics.numerics.maxPowerBalanceRelError.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.maxPowerBalanceRelError.Layout.Column = 2;

            % --- Non-Linear Routines Panel ---
            nonLinearPanel = uipanel(grid, 'Title', 'Non-Linear Routines');
            nonLinearPanel.Layout.Row = 2;
            nonLinearGrid = uigridlayout(nonLinearPanel, [6, 2]); % Rows, Columns - increased for advanced panel
            nonLinearGrid.ColumnWidth = {'fit', '1x'};
            nonLinearGrid.RowHeight = repmat({'fit'}, 1, 6);
            nonLinearGrid.Padding = [10 10 10 10];
            nonLinearGrid.RowSpacing = 5;

            row = 1;
            % Algorithm [cite: 12]
            uilabel(nonLinearGrid, 'Text', 'Non-Linear Algorithm:');
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.algorithm = uidropdown(nonLinearGrid, 'Items', {'mixingDirectSolutions', 'temporalIntegration'}, 'ValueChangedFcn', @(src, evt) gui.handleAlgorithmChange(evt.Value));

            row = row + 1;
            % Mixing Parameter [cite: 12] - Simplified slider with precision input
            uilabel(nonLinearGrid, 'Text', 'Mixing Parameter:');
            mixingParamGrid = uigridlayout(nonLinearGrid, [1, 2]);
            mixingParamGrid.ColumnWidth = {'0.6x', '0.4x'}; % Give more space to text field
            mixingParamGrid.Padding = [0 0 0 0];
            mixingParamGrid.ColumnSpacing = 5;
            
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameter = uislider(mixingParamGrid, 'Limits', [0, 1], 'Value', 0.7, 'MajorTicks', [0, 1], 'MajorTickLabels', {'0', '1'}, 'ValueChangedFcn', @(src, evt) gui.sliderMixingParameterChanged(src, evt.Value));
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise = uieditfield(mixingParamGrid, 'numeric', 'Limits', [0, 1], 'Value', 0.7, 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.textMixingParameterChanged(src, evt.Value));
            
            mixingParamGrid.Layout.Row = row;
            mixingParamGrid.Layout.Column = 2;

            row = row + 1;
            % Max EEDF Rel Error [cite: 12]
            uilabel(nonLinearGrid, 'Text', 'Max EEDF Rel. Error:');
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError = uieditfield(nonLinearGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.maxEedfRelError', evt.Value));

            row = row + 1;
            % Advanced Parameters Section (Expandable) - Only visible when algorithm is temporalIntegration
            advancedButton = uibutton(nonLinearGrid, 'Text', 'Advanced (Optional) ▼', 'ButtonPushedFcn', @(src, evt) gui.toggleAdvancedSection(src), 'Visible', 'off');
            advancedButton.Layout.Column = [1, 2];
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton = advancedButton;
            
            % Create expandable advanced panel (initially hidden)
            advancedPanel = uipanel(nonLinearGrid, 'Title', '', 'Visible', 'off');
            advancedPanel.Layout.Row = [row + 1, row + 2]; % Span two rows to accommodate content
            advancedPanel.Layout.Column = [1, 2];
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel = advancedPanel;
            
            % Create the content grid first
            contentGrid = uigridlayout(advancedPanel, [4, 2]); % 4 rows: title + 3 parameters
            contentGrid.ColumnWidth = {'fit', '1x'};
            contentGrid.RowHeight = {'fit', 'fit', 'fit', 'fit'};
            contentGrid.Padding = [10 10 10 10];
            contentGrid.RowSpacing = 5;
            
            % Create title row with checkbox
            titleRow = uigridlayout(contentGrid, [1, 2]);
            titleRow.ColumnWidth = {'1x', 'fit'};
            titleRow.Padding = [0 0 0 0];
            titleRow.RowSpacing = 0;
            titleRow.ColumnSpacing = 10;
            titleRow.Layout.Row = 1;
            titleRow.Layout.Column = [1, 2];
            
            % Title label
            titleLabel = uilabel(titleRow, 'Text', 'Advanced ODE Solver Parameters');
            titleLabel.FontWeight = 'bold';
            
                         % Checkbox to control if ODE parameters can be edited
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn = uicheckbox(titleRow, 'Text', '', 'Value', false, 'ValueChangedFcn', @(src, evt) gui.toggleOdeParametersEnable(evt.Value));
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Layout.Column = 2;
            
            % Absolute Tolerance
            uilabel(contentGrid, 'Text', 'Absolute Tolerance:');
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol', evt.Value));
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Layout.Row = 2;
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Layout.Column = 2;

                         % Relative Tolerance
             uilabel(contentGrid, 'Text', 'Relative Tolerance:');
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol', evt.Value));
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Layout.Row = 3;
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Layout.Column = 2;

             % Max Step
             uilabel(contentGrid, 'Text', 'Max Step:');
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'HorizontalAlignment', 'left', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep', evt.Value));
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep.Layout.Row = 4;
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep.Layout.Column = 2;



        end

        function createOutputPanel(gui, parent)
            % Create main grid with fixed row heights
            grid = uigridlayout(parent, [8, 3]); 
            grid.ColumnWidth = {'fit', '1x', 'fit'}; 
            grid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'}; % All fixed height instead of '1x' for last row
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 8;
            grid.ColumnSpacing = 10;

            % Row 1: GUI Settings Title
            guiTitleLabel = uilabel(grid, 'Text', 'GUI Settings:', 'FontWeight', 'bold');
            guiTitleLabel.Layout.Row = 1;
            guiTitleLabel.Layout.Column = 1;

            % Row 2: Enable GUI checkbox
            gui.UIControls.gui.isOn = uicheckbox(grid, ...
                'Text', 'Enable GUI', ...
                'ValueChangedFcn', @(src, evt) gui.toggleGuiControls(evt.Value));
            gui.UIControls.gui.isOn.Layout.Row = 2;
            gui.UIControls.gui.isOn.Layout.Column = 1;

            % Row 3: Refresh Frequency
            refreshFreqLabel = uilabel(grid, 'Text', 'Refresh Frequency:');
            refreshFreqLabel.Layout.Row = 3;
            refreshFreqLabel.Layout.Column = 1;
            gui.UIControls.gui.refreshFrequency = uieditfield(grid, 'numeric', ...
                'Limits', [1, Inf], ...
                'HorizontalAlignment', 'left', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'gui.refreshFrequency', evt.Value));
            gui.UIControls.gui.refreshFrequency.Layout.Row = 3;
            gui.UIControls.gui.refreshFrequency.Layout.Column = [2, 3];

            % Row 4: Output Settings Title  
            outputTitleLabel = uilabel(grid, 'Text', 'Output Settings:', 'FontWeight', 'bold');
            outputTitleLabel.Layout.Row = 4;
            outputTitleLabel.Layout.Column = 1;

            % Row 5: Enable Output checkbox
            gui.UIControls.output.isOn = uicheckbox(grid, ...
                'Text', 'Enable Output', ...
                'Value', false, ...
                'ValueChangedFcn', @(src, evt) gui.toggleOutputEnable(evt.Value));
            gui.UIControls.output.isOn.Layout.Row = 5;
            gui.UIControls.output.isOn.Layout.Column = 1;

            % Row 6: Data Format
            dataFormatLabel = uilabel(grid, 'Text', 'Data Format:');
            dataFormatLabel.Layout.Row = 6;
            dataFormatLabel.Layout.Column = 1;
            
            gui.UIControls.output.dataFormat = uidropdown(grid, ...
                'Items', {'txt', 'hdf5', 'hdf5+txt'}, ...
                'Enable', 'on', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'output.dataFormat', evt.Value));
            gui.UIControls.output.dataFormat.Layout.Row = 6;
            gui.UIControls.output.dataFormat.Layout.Column = [2, 3];

            % Row 7: Output Folder
            folderLabelGrid = uigridlayout(grid, [1, 2]);
            folderLabelGrid.ColumnWidth = {'fit', 'fit'};
            folderLabelGrid.RowHeight = {'fit'};
            folderLabelGrid.Padding = [0 0 0 0];
            folderLabelGrid.ColumnSpacing = 4;
            folderLabelGrid.Layout.Row = 7;
            folderLabelGrid.Layout.Column = 1;

            gui.UIControls.output.useCustomFolder = uicheckbox(folderLabelGrid, ...
                'Text', '', ...
                'Value', true, ...
                'ValueChangedFcn', @(src, evt) gui.toggleOutputFolderMode(evt.Value));
            gui.UIControls.output.useCustomFolder.Layout.Row = 1;
            gui.UIControls.output.useCustomFolder.Layout.Column = 1;

            folderLabel = uilabel(folderLabelGrid, 'Text', 'Output Folder:');
            folderLabel.Layout.Row = 1;
            folderLabel.Layout.Column = 2;
            gui.UIControls.output.folder = uieditfield(grid, 'text', ...
                'Enable', 'on', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'output.folder', evt.Value));
            gui.UIControls.output.folder.Layout.Row = 7;
            gui.UIControls.output.folder.Layout.Column = 2;
            
            gui.UIControls.output.browseButton = uibutton(grid, 'Text', 'Browse...', ...
                'ButtonPushedFcn', @gui.browseFolder);
            gui.UIControls.output.browseButton.Layout.Row = 7;
            gui.UIControls.output.browseButton.Layout.Column = 3;

            % Row 8: Data Sets - Use a more compact layout
            dataSetLabel = uilabel(grid, 'Text', 'Data Sets to Save:');
            dataSetLabel.Layout.Row = 8;
            dataSetLabel.Layout.Column = 1;
            dataSetLabel.VerticalAlignment = 'top'; % Align to top instead of center
            
            % Create checkbox panel with minimal padding
            checkboxPanel = uipanel(grid, 'Title', '', 'BorderType', 'none');
            checkboxPanel.Layout.Row = 8;
            checkboxPanel.Layout.Column = [2, 3];
            
            % Compact checkbox grid - adjust row height to be minimal
            checkboxGrid = uigridlayout(checkboxPanel, [2, 4]);
            checkboxGrid.ColumnWidth = repmat({'1x'}, 1, 4);
            checkboxGrid.RowHeight = {'fit', 'fit'}; % Both rows fit content
            checkboxGrid.Padding = [0 5 5 5]; % Minimal top padding
            checkboxGrid.RowSpacing = 3; % Reduced spacing
            checkboxGrid.ColumnSpacing = 8; % Reduced spacing
            
            % Create checkboxes
            dataSetNames = {'inputs', 'log', 'eedf', 'swarmParameters', ...
                        'rateCoefficients', 'powerBalance', 'lookUpTables'};
            
            for i = 1:length(dataSetNames)
                row_idx = ceil(i/4);
                col_idx = mod(i-1, 4) + 1;
                
                checkbox = uicheckbox(checkboxGrid, ...
                    'Text', dataSetNames{i}, ...
                    'Value', true, ...
                    'ValueChangedFcn', @(src, evt) gui.updateDataSetsSelection(src, dataSetNames{i}, evt.Value));
                checkbox.Layout.Row = row_idx;
                checkbox.Layout.Column = col_idx;
                
                gui.UIControls.output.dataSets.(dataSetNames{i}) = checkbox;
            end
        end
        
        function populateGUIFromSetup(gui)
            % Populate all UI controls with values from gui.Setup
            % --- Top-level fields (e.g., workingConditions, electronKinetics, output) ---
            topLevelFields = fieldnames(gui.UIControls);
            for i = 1:length(topLevelFields)
                sectionName = topLevelFields{i}; % e.g., 'workingConditions'

                % --- Second-level fields (controls directly under the section) ---
                controlsInSection = fieldnames(gui.UIControls.(sectionName));
                for j = 1:length(controlsInSection)
                    controlName = controlsInSection{j}; % e.g., 'reducedField' or 'numerics'

                    control = gui.UIControls.(sectionName).(controlName);

                    % Check if this is a nested structure of controls (like 'numerics')
                    if isstruct(control)
                        % --- Third-level fields (e.g., under 'numerics') ---
                        subSectionName = controlName; % e.g., 'numerics'
                        controlsInSubSection = fieldnames(gui.UIControls.(sectionName).(subSectionName));
                        for k = 1:length(controlsInSubSection)
                            subControlName = controlsInSubSection{k}; % e.g., 'maxPowerBalanceRelError' or 'energyGrid'
                            subControl = gui.UIControls.(sectionName).(subSectionName).(subControlName);

                            if isstruct(subControl)
                                % --- Fourth-level fields (e.g., under 'energyGrid' or 'nonLinearRoutines') ---
                                subSubSectionName = subControlName; % e.g., 'energyGrid'
                                controlsInSubSubSection = fieldnames(gui.UIControls.(sectionName).(subSectionName).(subSubSectionName));
                                for l = 1:length(controlsInSubSubSection)
                                    subSubControlName = controlsInSubSubSection{l}; % e.g., 'variableGrid'
                                    subSubControl = gui.UIControls.(sectionName).(subSectionName).(subSubSectionName).(subSubControlName);
                                    
                                    % Skip UI containers and controls that don't have corresponding Setup data
                                    if isa(subSubControl, 'matlab.ui.control.Button') || ...
                                       isa(subSubControl, 'matlab.ui.control.Panel') || ...
                                       strcmp(subSubControlName, 'mixingParameterPrecise') || ...
                                       strcmp(subSubControlName, 'advancedButton') || ...
                                       strcmp(subSubControlName, 'advancedPanel')
                                        continue;
                                    end
                                    
                                    % Skip output.dataSets fields during initial population to avoid warnings
                                    if strcmp(sectionName, 'output') && strcmp(subSectionName, 'dataSets')
                                        continue;
                                    end
                                    
                                    % Skip listboxes - they're handled specially below
                                    if isa(subSubControl, 'matlab.ui.control.ListBox')
                                        continue;
                                    end
                                    
                                    if isstruct(subSubControl)
                                        % --- Fifth-level fields (e.g., under 'smartGrid' or 'odeSetParameters') ---
                                        fifthLevelFields = fieldnames(subSubControl);
                                        for m = 1:length(fifthLevelFields)
                                            fifthLevelName = fifthLevelFields{m};
                                            fifthLevelControl = subSubControl.(fifthLevelName);
                                            
                                            % Skip UI containers
                                            if isa(fifthLevelControl, 'matlab.ui.control.Button') || ...
                                               isa(fifthLevelControl, 'matlab.ui.control.Panel')
                                                continue;
                                            end
                                            
                                            dataPath = sprintf('%s.%s.%s.%s.%s', sectionName, subSectionName, subSubSectionName, subSubControlName, fifthLevelName);
                                            gui.setControlValue(fifthLevelControl, dataPath);
                                        end
                                    else
                                        dataPath = sprintf('%s.%s.%s.%s', sectionName, subSectionName, subSubSectionName, subSubControlName);
                                        gui.setControlValue(subSubControl, dataPath);
                                    end
                                end
                            else % Control is at third level (e.g., 'maxPowerBalanceRelError')
                                % Skip listboxes - they're handled specially below
                                if isa(subControl, 'matlab.ui.control.ListBox')
                                    continue;
                                end
                                dataPath = sprintf('%s.%s.%s', sectionName, subSectionName, subControlName);
                                gui.setControlValue(subControl, dataPath);
                            end
                        end
                    else % Control is at second level (e.g., 'reducedField')
                        % Skip listboxes - they're handled specially below
                        if isa(control, 'matlab.ui.control.ListBox')
                            continue;
                        end
                        % Skip totalSccmOutFlow - it's handled specially below
                        if strcmp(sectionName, 'workingConditions') && strcmp(controlName, 'totalSccmOutFlow')
                            continue;
                        end
                        % Skip output controls that don't map directly to input data
                        if strcmp(sectionName, 'output') && ...
                                (strcmp(controlName, 'browseButton') || strcmp(controlName, 'useCustomFolder'))
                            continue;
                        end
                        dataPath = sprintf('%s.%s', sectionName, controlName);
                        gui.setControlValue(control, dataPath);
                    end
                end
            end

            % Special handling for list boxes (setting Items source and Value)
            % These paths need to match your gui.Setup structure exactly
            try
                gui.UIControls.electronKinetics.LXCatFiles.Items = gui.Setup.electronKinetics.LXCatFiles;
                if ~isempty(gui.UIControls.electronKinetics.LXCatFiles.Items)
                    gui.UIControls.electronKinetics.LXCatFiles.Value = {}; % Clear selection initially
                end
            catch ME
                warning('Error setting LXCatFiles list: %s', ME.message);
            end
            try
                gui.UIControls.electronKinetics.LXCatExtraFiles.Items = gui.Setup.electronKinetics.LXCatExtraFiles;
                if ~isempty(gui.UIControls.electronKinetics.LXCatExtraFiles.Items)
                    gui.UIControls.electronKinetics.LXCatExtraFiles.Value = {}; % Clear selection initially
                end
            catch ME
                warning('Error setting LXCatExtraFiles list: %s', ME.message);
            end
            try
                gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Items = gui.Setup.electronKinetics.effectiveCrossSectionPopulations;
                if ~isempty(gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Items)
                    gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Value = {}; % Clear selection initially
                end
            catch ME
                warning('Error setting effectiveCrossSectionPopulations list: %s', ME.message);
            end
            try
                gasFileFields = {'mass', 'harmonicFrequency', 'anharmonicFrequency', ...
                    'rotationalConstant', 'electricQuadrupoleMoment', 'OPBParameter'};
                for idx = 1:length(gasFileFields)
                    fieldName = gasFileFields{idx};
                    if isfield(gui.UIControls.gasProperties, fieldName) && ...
                            isfield(gui.Setup.electronKinetics.gasProperties, fieldName)
                        items = gui.asCellArray(gui.Setup.electronKinetics.gasProperties.(fieldName));
                        gui.UIControls.gasProperties.(fieldName).Items = items;
                        gui.UIControls.gasProperties.(fieldName).Value = {};
                    end
                end
            catch ME
                warning('Error setting gas property file lists: %s', ME.message);
            end
            try
                gui.UIControls.gasProperties.fraction.Items = gui.Setup.electronKinetics.gasProperties.fraction;
                if ~isempty(gui.UIControls.gasProperties.fraction.Items)
                    gui.UIControls.gasProperties.fraction.Value = {}; % Clear selection initially
                end
            catch ME
                warning('Error setting gas fractions list: %s', ME.message);
            end
            try
                carItems = gui.Setup.electronKinetics.CARgases;
                if isempty(carItems)
                    gui.UIControls.electronKinetics.CARgases.Items = {};
                else
                    gui.UIControls.electronKinetics.CARgases.Items = carItems;
                    gui.UIControls.electronKinetics.CARgases.Value = {}; % Clear selection initially
                    % Enable the checkbox and controls since we have default items
                    if isfield(gui.UIControls.electronKinetics, 'CARcheckbox')
                        gui.UIControls.electronKinetics.CARcheckbox.Value = true;
                        gui.toggleCARGasEnable(true);
                    end
                end
            catch ME
                warning('Error setting CARgases list: %s', ME.message);
            end
            try
                gui.UIControls.electronKinetics.stateProperties.population.Items = gui.Setup.electronKinetics.stateProperties.population;
                if ~isempty(gui.UIControls.electronKinetics.stateProperties.population.Items)
                    gui.UIControls.electronKinetics.stateProperties.population.Value = {}; % Clear selection initially
                end
            catch ME
                warning('Error setting state populations list: %s', ME.message);
            end
            try
                items = gui.Setup.electronKinetics.stateProperties.energy;
                if ~isempty(items)
                    gui.UIControls.electronKinetics.stateProperties.energy.Items = items;
                    gui.UIControls.electronKinetics.stateProperties.energy.Value = {};
                end
                % If empty, leave it as is - Items will remain empty
            catch ME
                warning('Error setting state energy list: %s', ME.message);
            end
            try
                items = gui.Setup.electronKinetics.stateProperties.statisticalWeight;
                if ~isempty(items)
                    gui.UIControls.electronKinetics.stateProperties.statisticalWeight.Items = items;
                    gui.UIControls.electronKinetics.stateProperties.statisticalWeight.Value = {};
                end
                % If empty, leave it as is - Items will remain empty
            catch ME
                warning('Error setting state statistical weight list: %s', ME.message);
            end
            try
                % Configure data sets checkboxes - skip during initial population to avoid warnings
                if isfield(gui.UIControls, 'output') && isfield(gui.UIControls.output, 'dataSets')
                    try
                        if isfield(gui.Setup, 'output') && isfield(gui.Setup.output, 'dataSets')
                            selectedDataSets = gui.Setup.output.dataSets;
                            dataSetNames = {'inputs', 'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTables'};
                            
                            for i = 1:length(dataSetNames)
                                dataSetName = dataSetNames{i};
                                if isfield(gui.UIControls.output.dataSets, dataSetName)
                                    checkbox = gui.UIControls.output.dataSets.(dataSetName);
                                    checkbox.Value = ismember(dataSetName, selectedDataSets);
                                end
                            end
                        end
                    catch ME
                        % Skip if dataSets structure doesn't exist yet
                    end
                end
            catch ME
                warning('Error configuring data sets checkboxes: %s', ME.message);
            end

            % Configure initial state of smart grid controls
            try
                gui.toggleSmartGridEnable(gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.isOn);
            catch ME
                warning('Error configuring smart grid controls: %s', ME.message);
            end

            % Configure initial state of advanced button visibility
            % Button is only visible when algorithm is temporalIntegration
            try
                advancedButton = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton;
                advancedPanel = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel;
                
                % Check current algorithm
                algorithm = gui.Setup.electronKinetics.numerics.nonLinearRoutines.algorithm;
                
                if strcmp(algorithm, 'temporalIntegration')
                    % Show button when algorithm is temporalIntegration
                    advancedButton.Visible = 'on';
                    
                    % Check if odeSetParameters exists to determine if panel should be expanded
                    hasOdeParams = isfield(gui.Setup.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters') && ...
                                   isstruct(gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters) && ...
                                   (isfield(gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'AbsTol') || ...
                                    isfield(gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'RelTol') || ...
                                    isfield(gui.Setup.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'MaxStep'));
                    
                    if hasOdeParams
                        % Expand panel if odeSetParameters exists
                        advancedPanel.Visible = 'on';
                        advancedButton.Text = 'Advanced (Optional) ▲';
                    else
                        % Collapse panel if odeSetParameters doesn't exist
                        advancedPanel.Visible = 'off';
                        advancedButton.Text = 'Advanced (Optional) ▼';
                    end
                else
                    % Hide button when algorithm is not temporalIntegration
                    advancedButton.Visible = 'off';
                    advancedPanel.Visible = 'off';
                    advancedButton.Text = 'Advanced (Optional) ▼';
                end
            catch ME
                warning('Error configuring advanced button visibility: %s', ME.message);
            end

            % Configure initial state of ODE parameters (initially disabled)
            try
                gui.toggleOdeParametersEnable(false);
            catch ME
                warning('Error configuring ODE parameters: %s', ME.message);
            end

            % Configure initial state of GUI controls
            try
                gui.toggleGuiControls(gui.Setup.gui.isOn);
            catch ME
                warning('Error configuring GUI controls: %s', ME.message);
            end

            % Configure initial state of Working Conditions optional fields (inactive by default)
            try
                gui.toggleGasPressureEnable(false);
                gui.toggleElectronDensityEnable(false);
                gui.toggleDischargeCurrentEnable(false);
                gui.toggleDischargePowerEnable(false);
                gui.activateGasPressureForExcitationFrequency(gui.Setup.workingConditions.excitationFrequency);
            catch ME
                warning('Error configuring optional working conditions: %s', ME.message);
            end

            % Configure initial state of Electron Kinetics (enabled but inactive by default)
            try
                gui.toggleElectronKineticsEnable(true);
                % Ensure growth model respects initial ionization operator choice
                gui.handleIonizationOperatorChange(gui.Setup.electronKinetics.ionizationOperatorType);
                % Ensure e-e collisions state is applied (unlocks Working Conditions fields if active)
                gui.handleEECollisionsChange(gui.Setup.electronKinetics.includeEECollisions);
            catch ME
                warning('Error configuring electron kinetics controls: %s', ME.message);
            end

            % Configure initial state of Output controls (respect Enable Output flag)
            try
                gui.toggleOutputEnable(gui.Setup.output.isOn);
            catch ME
                warning('Error configuring output controls: %s', ME.message);
            end

            
            % Configure data sets checkboxes after all other controls are populated
            try
                if isfield(gui.UIControls, 'output') && isfield(gui.UIControls.output, 'dataSets')
                    if isfield(gui.Setup, 'output') && isfield(gui.Setup.output, 'dataSets')
                        selectedDataSets = gui.Setup.output.dataSets;
                        dataSetNames = {'inputs', 'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTables'};
                        
                        for i = 1:length(dataSetNames)
                            dataSetName = dataSetNames{i};
                            if isfield(gui.UIControls.output.dataSets, dataSetName)
                                checkbox = gui.UIControls.output.dataSets.(dataSetName);
                                checkbox.Value = ismember(dataSetName, selectedDataSets);
                            end
                        end
                    end
                end
            catch ME
                % Skip if dataSets structure doesn't exist yet
            end

            % Configure initial state of EEDF type visibility
            try
                eedfType = gui.Setup.electronKinetics.eedfType;
                gui.handleEedfTypeChangeForInit(eedfType);
            catch ME
                warning('Error configuring EEDF type visibility: %s', ME.message);
            end
            
            % Sync paired controls that are skipped by generic population.
            gui.syncMixingParameterControls();

            % Update run button state after initial population
            gui.updateRunButtonState();
        end

        function setControlValue(gui, control, dataPath)
            % Helper function to get value from Setup and set control value
            
            % Skip output.dataSets fields completely to avoid warnings
            if contains(dataPath, 'output.dataSets.')
                return;
            end
            
            % Skip UI-only controls that don't exist in Setup
            if contains(dataPath, 'tabs.tabGroup') || ...
               contains(dataPath, 'tabs.kineticsTab') || ...
               startsWith(dataPath, 'runButton') || ...
               startsWith(dataPath, 'saveInputButton') || ...
               contains(dataPath, 'totalSccmOutFlowType') || ...
               contains(dataPath, 'workingConditions.infoTitleLabel') || ...
               contains(dataPath, 'workingConditions.infoDescLabel') || ...
               contains(dataPath, 'workingConditions.reducedFieldMessage') || ...
               contains(dataPath, 'workingConditions.reducedFieldBorderPanel') || ...
               contains(dataPath, 'workingConditions.electronTemperatureMessage') || ...
               contains(dataPath, 'workingConditions.electronTemperatureBorderPanel') || ...
               contains(dataPath, 'workingConditions.infoContainer') || ...
               contains(dataPath, 'workingConditions.infoImage') || ...
               contains(dataPath, 'workingConditions.totalSccmOutFlowPanel') || ...
               contains(dataPath, 'gasPressureCheckbox') || ...
               contains(dataPath, 'electronDensityCheckbox') || ...
               contains(dataPath, 'LXCatExtraCheckbox') || ...
               contains(dataPath, 'effectivePopCheckbox') || ...
               contains(dataPath, 'CARcheckbox') || ...
               contains(dataPath, 'CARaddButton') || ...
               contains(dataPath, 'CARremoveButton') || ...
               contains(dataPath, 'odeSetParameters.isOn') || ...
               contains(dataPath, 'smartGrid.isOn') || ...
               contains(dataPath, 'electronKinetics.shapeParameterPrecise') || ...
               contains(dataPath, 'electronKinetics.shapeParameterLabel')
                return;
            end
            
            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(dataPath, 'gasProperties.')
                dataPath = strrep(dataPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end
            
            try
                value = gui.getNestedField(gui.Setup, dataPath);
                if isempty(value) && ~ischar(value) && ~iscell(value) % Allow empty char and empty cell arrays
                     % Only show debug message for non-output fields to reduce noise
                     if ~contains(dataPath, 'output.dataSets.') && ...
                        ~contains(dataPath, 'dischargeCurrent') && ...
                        ~contains(dataPath, 'dischargePowerDensity')
                          fprintf('Debug: No value found for %s in Setup struct.\n', dataPath);
                     end
                     return;
                end

                if isa(control, 'matlab.ui.control.ListBox')
                    % For list boxes, we need to set Items first
                    if strcmp(dataPath, 'electronKinetics.LXCatFiles') || ...
                       strcmp(dataPath, 'electronKinetics.LXCatExtraFiles') || ...
                       strcmp(dataPath, 'electronKinetics.effectiveCrossSectionPopulations') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.fraction') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.mass') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.harmonicFrequency') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.anharmonicFrequency') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.rotationalConstant') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.electricQuadrupoleMoment') || ...
                       strcmp(dataPath, 'electronKinetics.gasProperties.OPBParameter') || ...
                       strcmp(dataPath, 'electronKinetics.stateProperties.population') || ...
                       strcmp(dataPath, 'electronKinetics.stateProperties.energy') || ...
                       strcmp(dataPath, 'electronKinetics.stateProperties.statisticalWeight') || ...
                       strcmp(dataPath, 'electronKinetics.CARgases')
                        control.Items = value;
                        control.Value = {}; % Clear selection initially
                    else
                        % For other list boxes, just set the value
                        control.Value = value;
                    end
                elseif isa(control, 'matlab.ui.control.DropDown')
                    % Items set separately, here we set the selected value
                    control.Value = value;
                elseif isa(control, 'matlab.ui.control.CheckBox')
                    control.Value = logical(value); % Ensure its logical
                elseif isa(control, 'matlab.ui.control.NumericEditField')
                    if isnumeric(value)
                        control.Value = value;
                    else % Handle non-numeric stored value if necessary
                        numVal = str2double(value);
                        if isnan(numVal)
                            control.Value = 0;
                        else
                            control.Value = numVal;
                        end
                        warning('Non-numeric value found for numeric field %s. Attempted conversion.', dataPath);
                    end
                elseif isa(control, 'matlab.ui.control.EditField') % Text edit field
                    if isnumeric(value)
                        control.Value = num2str(value);
                    elseif ischar(value) || isstring(value)
                        control.Value = value;
                    elseif islogical(value)
                        if value
                            control.Value = 'true';
                        else
                            control.Value = 'false';
                        end
                    else
                        % Try to convert other types (like arrays) to string
                        try
                            control.Value = mat2str(value);
                        catch
                            control.Value = ''; % Fallback
                            warning('Could not convert value for %s to string.', dataPath);
                        end
                    end
                elseif isa(control, 'matlab.ui.control.Slider')
                    % For sliders, we need to handle them specially
                    if strcmp(dataPath, 'electronKinetics.numerics.nonLinearRoutines.mixingParameter')
                        % This is handled by the mixing parameter functions
                        return;
                    end
                    % For other sliders, try to set the value
                    if isnumeric(value)
                        control.Value = value;
                    end
                elseif isa(control, 'matlab.ui.control.Button') || isa(control, 'matlab.ui.control.Panel')
                    % These are UI containers, not data controls
                    return;
                elseif isstruct(control)
                    % Handle nested struct controls (like smartGrid.isOn, odeSetParameters.isOn)
                    % Skip these during initial population to avoid warnings
                    return;
                else
                    fprintf('Debug: Control type for %s not explicitly handled in setControlValue.\n', dataPath);
                end
            catch ME
                warning('Error setting control value for path %s: %s', dataPath, ME.message);
            end
        end


        % --- Callback Functions ---

        function updateField(gui, control, fieldPath, value)
            % Update a field in the Setup struct using dot notation path

            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(fieldPath, 'gasProperties.')
                fieldPath = strrep(fieldPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end

            % Check if the control supports BackgroundColor for visual feedback
            supportsBackgroundColor = isa(control, 'matlab.ui.control.EditField') || ...
                                      isa(control, 'matlab.ui.control.NumericEditField'); % Add other types if they support it

            if supportsBackgroundColor
                originalColor = control.BackgroundColor; % Store original color
            end

            try
                % Handle specific types if necessary (e.g., numeric conversion)
                currentValue = gui.getNestedField(gui.Setup, fieldPath);
                % --- Type Conversion Logic ---
                % Check if the target field in Setup is numeric
                if isnumeric(currentValue) && ~isa(value, 'logical') % Don't convert logicals to numeric
                    if ischar(value) || isstring(value) % Handle text input for numeric fields
                        numericValue = str2double(value);
                        if isnan(numericValue)
                            error('Invalid numeric input: "%s"', value);
                        end
                        value = numericValue; % Use the converted value
                    elseif ~isnumeric(value) % If it's neither string nor numeric (e.g., unexpected type)
                        error('Expected numeric or string input for numeric field, got %s', class(value));
                    end
                % Add elseif blocks here if other specific type conversions are needed
                end
                % --- End Type Conversion ---

                % Set the nested field in the Setup struct
                gui.setNestedField(fieldPath, value);

                % Provide visual feedback if supported
                if supportsBackgroundColor
                    control.BackgroundColor = [0.9, 1.0, 0.9]; % Greenish tint on success
                    drawnow; % Ensure color update is visible briefly
                    pause(0.1);
                    control.BackgroundColor = originalColor; % Restore original color
                end

            catch ME
                warning('Error updating field "%s": %s', fieldPath, ME.message);
                % Indicate error on the control if supported
                if supportsBackgroundColor
                    control.BackgroundColor = [1.0, 0.8, 0.8]; % Reddish tint on error
                else
                    % Alternative feedback for unsupported controls (e.g., brief message)
                    origTooltip = control.Tooltip;
                    control.Tooltip = sprintf('Error: %s', ME.message);
                    pause(1.5); % Show tooltip briefly
                    control.Tooltip = origTooltip;
                end
                % Optional: Restore the previous valid value from gui.Setup to the control
                % try
                %    previousValue = gui.getNestedField(gui.Setup, fieldPath);
                %    control.Value = previousValue; % Revert UI (careful with control types)
                % catch % Ignore errors during revert
                % end
                return; % Stop further processing if update failed
            end

            % Optional: Re-enable/disable controls based on the change
            if strcmp(fieldPath, 'electronKinetics.numerics.energyGrid.smartGrid.isOn')
                gui.toggleSmartGridEnable(value);
            end

            if strcmp(fieldPath, 'workingConditions.reducedField') || ...
               strcmp(fieldPath, 'workingConditions.electronTemperature')
                gui.updateRunButtonState();
            end
        end

        function handleExcitationFrequencyChange(gui, control, value)
            % Update excitation frequency and require gas pressure for RF cases.
            gui.updateField(control, 'workingConditions.excitationFrequency', value);
            gui.activateGasPressureForExcitationFrequency(value);
        end

        function activateGasPressureForExcitationFrequency(gui, excitationFrequency)
            if ~isempty(excitationFrequency) && isnumeric(excitationFrequency) && any(excitationFrequency(:) ~= 0) && ...
               isfield(gui.UIControls, 'workingConditions') && ...
               isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                gui.UIControls.workingConditions.gasPressureCheckbox.Value = true;
                gui.toggleGasPressureEnable(true);
                gui.setTemporalGrowthModelForExcitationFrequency();
                gui.UIControls.electronKinetics.growthModelType.Enable = 'off';
            elseif isfield(gui.UIControls, 'workingConditions') && ...
                   isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                gui.UIControls.workingConditions.gasPressureCheckbox.Value = false;
                gui.toggleGasPressureEnable(false);
                gui.UIControls.electronKinetics.growthModelType.Enable = 'on';
            end
        end

        function setTemporalGrowthModelForExcitationFrequency(gui)
            if isfield(gui.UIControls, 'electronKinetics') && ...
               isfield(gui.UIControls.electronKinetics, 'growthModelType')
                gui.UIControls.electronKinetics.growthModelType.Value = 'temporal';
            end

            if isfield(gui.Setup, 'electronKinetics') && ...
               isfield(gui.Setup.electronKinetics, 'growthModelType') && ...
               ~strcmp(gui.Setup.electronKinetics.growthModelType, 'temporal')
                gui.setNestedField('electronKinetics.growthModelType', 'temporal');
            end
        end
        
        function handleAlgorithmChange(gui, algorithmValue)
            % Handle algorithm change: show/hide Advanced button based on algorithm
            % Also update the Setup struct
            
            % Update the Setup struct
            gui.updateField([], 'electronKinetics.numerics.nonLinearRoutines.algorithm', algorithmValue);
            
            % Show/hide Advanced button based on algorithm
            try
                advancedButton = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton;
                advancedPanel = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel;
                
                if strcmp(algorithmValue, 'temporalIntegration')
                    % Show button when algorithm is temporalIntegration
                    advancedButton.Visible = 'on';
                else
                    % Hide button and collapse panel when algorithm is not temporalIntegration
                    advancedButton.Visible = 'off';
                    advancedPanel.Visible = 'off';
                    advancedButton.Text = 'Advanced (Optional) ▼';
                    % Also uncheck the checkbox and disable controls
                    if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn')
                        gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Value = false;
                        gui.toggleOdeParametersEnable(false);
                    end
                end
            catch ME
                warning('Error handling algorithm change: %s', ME.message);
            end
        end

        function setNestedField(gui, fieldPath, value)
            % Sets a value in a nested struct using a dot-separated path string
            parts = strsplit(fieldPath, '.');
            s = gui.Setup; % Start with the main struct

            % Traverse the structure except for the last part
            for i = 1:length(parts)-1
                if isfield(s, parts{i})
                    s = s.(parts{i});
                else
                    % If a field doesnt exist, maybe create it (careful!)
                    error('Field "%s" not found in path "%s"', parts{i}, fieldPath);
                    % Alternatively, create nested structs if needed:
                    % s.(parts{i}) = struct();
                    % s = s.(parts{i});
                end
            end

            % Set the value of the final field
            % Use eval to set the nested field dynamically (use with caution)
            assignCmd = sprintf('gui.Setup.%s = value;', fieldPath);
            eval(assignCmd);
            % Alternative (safer if possible):
            % s.(parts{end}) = value; % This only works if 's' points to the correct sub-struct
            % To make the alternative work, you need to pass sub-structs by reference,
            % which MATLAB doesnt do directly for structs. Handle classes would work.
            % Or reconstruct the assignment.
        end

        function value = getNestedField(gui, startStruct, fieldPath)
            % Gets a value from a nested struct using a dot-separated path string
            parts = strsplit(fieldPath, '.');
            currentValue = startStruct;
            try
                for i = 1:length(parts)
                    currentValue = currentValue.(parts{i});
                end
                value = currentValue;
            catch ME
                warning('Could not retrieve field: %s', fieldPath);
                value = []; % Return empty or default if not found
            end
        end



        function toggleSmartGridEnable(gui, shouldEnable)
            % Function to enable/disable smart grid related controls
            controls = {
                gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay,
                gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay,
                gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor
            };
            
            labels = {
                findobj(controls{1}.Parent, 'Text', 'Min EEDF Decay:'),
                findobj(controls{2}.Parent, 'Text', 'Max EEDF Decay:'),
                findobj(controls{3}.Parent, 'Text', 'Update Factor:')
            };
            
            for i = 1:length(controls)
                if shouldEnable
                    controls{i}.Enable = 'on';
                    if ~isempty(labels{i})
                        labels{i}.Enable = 'on';
                    end
                else
                    controls{i}.Enable = 'off';
                    if ~isempty(labels{i})
                        labels{i}.Enable = 'off';
                    end
                end
            end
            
            % Also update the Setup struct to reflect the checkbox state
            % We'll add a temporary isOn field just for the checkbox state
            if shouldEnable
                gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.isOn = true;
            else
                gui.Setup.electronKinetics.numerics.energyGrid.smartGrid.isOn = false;
            end
        end

        function toggleAdvancedSection(gui, button)
            % Function to expand/collapse the advanced parameters section
            advancedPanel = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel;
               
            if strcmp(advancedPanel.Visible, 'off')
                % Expand the panel
                advancedPanel.Visible = 'on';
                button.Text = 'Advanced (Optional) ▲';
            else
                % Collapse the panel
                advancedPanel.Visible = 'off';
                button.Text = 'Advanced (Optional) ▼';
            end
        end

        function toggleOdeParametersEnable(gui, shouldEnable)
            % Function to enable/disable ODE parameter editing
            controls = {
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol,
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol,
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep
            };
            
            labels = {
                findobj(controls{1}.Parent, 'Text', 'Absolute Tolerance:'),
                findobj(controls{2}.Parent, 'Text', 'Relative Tolerance:'),
                findobj(controls{3}.Parent, 'Text', 'Max Step:')
            };
            
            for i = 1:length(controls)
                if shouldEnable
                    controls{i}.Enable = 'on';
                    if ~isempty(labels{i})
                        labels{i}.Enable = 'on';
                    end
                else
                    controls{i}.Enable = 'off';
                    if ~isempty(labels{i})
                        labels{i}.Enable = 'off';
                    end
                end
            end
        end



        function syncMixingParameterControls(gui)
            try
                value = gui.Setup.electronKinetics.numerics.nonLinearRoutines.mixingParameter;
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameter.Value = value;
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise.Value = value;
            catch ME
                warning('Error syncing mixing parameter controls: %s', ME.message);
            end
        end

        function sliderMixingParameterChanged(gui, ~, value)
            % Callback for slider value change
            try
                % Update the text field
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise.Value = value;
                % Update the setup struct
                gui.updateField(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameter, 'electronKinetics.numerics.nonLinearRoutines.mixingParameter', value);
            catch ME
                warning('Error updating mixing parameter from slider: %s', ME.message);
            end
        end

        function textMixingParameterChanged(gui, ~, value)
            % Callback for text field value change
            try
                % Validate the value
                if value < 0 || value > 1
                    warning('Mixing parameter must be between 0 and 1. Value reset to 0.7');
                    value = 0.7;
                    gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise.Value = value;
                end
                % Update the slider
                gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameter.Value = value;
                % Update the setup struct
                gui.updateField(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise, 'electronKinetics.numerics.nonLinearRoutines.mixingParameter', value);
            catch ME
                warning('Error updating mixing parameter from text field: %s', ME.message);
            end
        end

        function sliderShapeParameterChanged(gui, ~, value)
            % Callback for slider value change
            try
                % Update the text field
                gui.UIControls.electronKinetics.shapeParameterPrecise.Value = value;
                % Update the setup struct
                gui.updateField(gui.UIControls.electronKinetics.shapeParameter, 'electronKinetics.shapeParameter', value);
                gui.updateField(gui.UIControls.electronKinetics.shapeParameterPrecise, 'electronKinetics.shapeParameter', value);
            catch ME
                warning('Error updating shape parameter from slider: %s', ME.message);
            end
        end

        function textShapeParameterChanged(gui, ~, value)
            % Callback for text field value change
            try
                % Validate the value
                if value < 1 || value > 2
                    warning('Shape parameter must be between 1 and 2. Value reset to 1');
                    value = 1;
                    gui.UIControls.electronKinetics.shapeParameterPrecise.Value = value;
                end
                % Update the slider
                gui.UIControls.electronKinetics.shapeParameter.Value = value;
                % Update the setup struct
                gui.updateField(gui.UIControls.electronKinetics.shapeParameter, 'electronKinetics.shapeParameter', value);
                gui.updateField(gui.UIControls.electronKinetics.shapeParameterPrecise, 'electronKinetics.shapeParameter', value);
            catch ME
                warning('Error updating shape parameter from text field: %s', ME.message);
            end
        end



        function [newItem, wasCancelled] = promptStatePropertyItem(gui, fieldPath, itemToEdit, dlgtitle)
            newItem = '';
            wasCancelled = true;

            initialState = '';
            initialValue = '';
            initialFile = '';
            if ~isempty(itemToEdit)
                equalSignIndex = strfind(itemToEdit, '=');
                if ~isempty(equalSignIndex)
                    splitIndex = equalSignIndex(end);
                    initialState = strtrim(itemToEdit(1:splitIndex-1));
                    initialValue = strtrim(itemToEdit(splitIndex+1:end));
                else
                    initialFile = itemToEdit;
                end
            end

            if startsWith(dlgtitle, 'Edit')
                promptText = ['Edit item for ', fieldPath, ':'];
            else
                promptText = ['Enter new item for ', fieldPath, ':'];
            end

            oldState = gui.Fig.WindowState;
            screenSize = get(0, 'ScreenSize');
            dlgWidth = 620;
            dlgHeight = 215;
            dlgX = (screenSize(3) - dlgWidth) / 2;
            dlgY = (screenSize(4) - dlgHeight) / 2;

            dlg = uifigure('Name', dlgtitle, ...
                'WindowStyle', 'modal', ...
                'Position', [dlgX, dlgY, dlgWidth, dlgHeight], ...
                'Resize', 'off', ...
                'KeyPressFcn', @(src, evt) handleDialogKeyPress(evt), ...
                'CloseRequestFcn', @(src, evt) cancelDialog());

            grid = uigridlayout(dlg, [6, 3]);
            grid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
            grid.ColumnWidth = {'1x', 'fit', '1x'};
            grid.Padding = [18 16 18 16];
            grid.RowSpacing = 10;
            grid.ColumnSpacing = 10;

            promptLabel = uilabel(grid, 'Text', promptText, 'FontSize', 12);
            promptLabel.Layout.Row = 1;
            promptLabel.Layout.Column = [1, 3];

            stateLabel = uilabel(grid, 'Text', 'State', 'FontSize', 12, 'FontWeight', 'bold');
            stateLabel.Layout.Row = 2;
            stateLabel.Layout.Column = 1;


            valueLabel = uilabel(grid, 'Text', gui.getStatePropertyValueLabel(fieldPath), 'FontSize', 12, 'FontWeight', 'bold');
            valueLabel.Layout.Row = 2;
            valueLabel.Layout.Column = 3;

            stateField = uieditfield(grid, 'text', 'Value', initialState, 'FontSize', 12);
            stateField.Layout.Row = 3;
            stateField.Layout.Column = 1;

            equalsValueLabel = uilabel(grid, 'Text', '=', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            equalsValueLabel.Layout.Row = 3;
            equalsValueLabel.Layout.Column = 2;

            valueField = uieditfield(grid, 'text', 'Value', initialValue, 'FontSize', 12);
            valueField.Layout.Row = 3;
            valueField.Layout.Column = 3;

            fileLabel = uilabel(grid, 'Text', 'or call a text file:', 'FontSize', 12);
            fileLabel.Layout.Row = 4;
            fileLabel.Layout.Column = [1, 3];

            fileField = uieditfield(grid, 'text', 'Value', initialFile, 'FontSize', 12);
            fileField.Layout.Row = 5;
            fileField.Layout.Column = [1, 3];

            buttonGrid = uigridlayout(grid, [1, 3]);
            buttonGrid.Layout.Row = 6;
            buttonGrid.Layout.Column = [1, 3];
            buttonGrid.ColumnWidth = {'1x', 'fit', 'fit'};
            buttonGrid.Padding = [0 0 0 0];
            buttonGrid.ColumnSpacing = 8;

            okButton = uibutton(buttonGrid, 'Text', 'OK', 'ButtonPushedFcn', @(src, evt) acceptDialog());
            okButton.Layout.Row = 1;
            okButton.Layout.Column = 2;

            cancelButton = uibutton(buttonGrid, 'Text', 'Cancel', 'ButtonPushedFcn', @(src, evt) cancelDialog());
            cancelButton.Layout.Row = 1;
            cancelButton.Layout.Column = 3;

            drawnow;
            figure(dlg);
            uiwait(dlg);

            try
                gui.Fig.WindowState = oldState;
                figure(gui.Fig);
                drawnow;
            catch
            end

            function acceptDialog()
                fileValue = strtrim(fileField.Value);
                if ~isempty(fileValue)
                    newItem = fileValue;
                    wasCancelled = false;
                    delete(dlg);
                    return;
                end

                stateValue = strtrim(stateField.Value);
                propertyValue = strtrim(valueField.Value);
                if isempty(stateValue) || isempty(propertyValue)
                    uialert(dlg, ['State and ', gui.getStatePropertyValueLabel(fieldPath), ' cannot be empty.'], 'Invalid Input');
                    return;
                end
                newItem = [stateValue, ' = ', propertyValue];
                wasCancelled = false;
                delete(dlg);
            end

            function cancelDialog()
                wasCancelled = true;
                if isvalid(dlg)
                    delete(dlg);
                end
            end

            function handleDialogKeyPress(evt)
                if strcmp(evt.Key, 'return') || strcmp(evt.Key, 'enter')
                    acceptDialog();
                end
            end
        end

        function isStructured = isStructuredStatePropertyField(~, fieldPath)
            structuredFields = {
                'electronKinetics.stateProperties.energy', ...
                'electronKinetics.stateProperties.statisticalWeight', ...
                'electronKinetics.stateProperties.population' ...
            };
            isStructured = any(strcmp(fieldPath, structuredFields));
        end

        function label = getStatePropertyValueLabel(~, fieldPath)
            if strcmp(fieldPath, 'electronKinetics.stateProperties.energy')
                label = 'Energy';
            elseif strcmp(fieldPath, 'electronKinetics.stateProperties.statisticalWeight')
                label = 'Statistical Weight';
            elseif strcmp(fieldPath, 'electronKinetics.stateProperties.population')
                label = 'Population';
            else
                label = 'Value';
            end
        end

        function addListItem(gui, fieldPath, browseFile)
            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(fieldPath, 'gasProperties.')
                fieldPath = strrep(fieldPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end
            
            % For UIControls access, remove electronKinetics prefix if present
            uiPath = fieldPath;
            if startsWith(uiPath, 'electronKinetics.gasProperties.')
                uiPath = strrep(uiPath, 'electronKinetics.gasProperties.', 'gasProperties.');
            end
            
            listBox = gui.getNestedField(gui.UIControls, uiPath);
            
            % Check if listBox is empty (path not found)
            if isempty(listBox) || ~isobject(listBox)
                uialert(gui.Fig, 'Could not find the list control.', 'Error');
                return;
            end
            
            currentItems = listBox.Items;

            newItem = '';
            if browseFile
                % Start in Input folder
                inputFolder = fullfile(pwd, 'Input');
                if ~isfolder(inputFolder)
                    inputFolder = pwd;
                end
                
                % Save window state and bring to front
                drawnow;
                oldState = gui.Fig.WindowState;
                drawnow;
                pause(0.05);
                
                [file, path] = uigetfile('*.*', ['Select file for ', fieldPath], inputFolder);
                if isequal(file, 0) || isequal(path, 0)
                    % Restore window state after dialog cancellation
                    drawnow;
                    pause(0.1);
                    gui.Fig.WindowState = 'maximized';
                    figure(gui.Fig);
                    drawnow;
                    return; % User cancelled
                end
                
                % Restore window state after file dialog
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = oldState; % Restore previous (typically maximized)
                figure(gui.Fig); % Bring to front
                drawnow;
                fullPath = fullfile(path, file);
                
                % Convert to relative path from Input folder
                inputFolderFull = strrep(fullfile(pwd, 'Input'), '\', '/');
                fullPathNormalized = strrep(fullPath, '\', '/');
                inputFolderNormalized = strrep(inputFolderFull, '\', '/');
                
                if contains(fullPathNormalized, inputFolderNormalized)
                    % Extract the relative path
                    relativePath = erase(fullPathNormalized, inputFolderNormalized);
                    % Remove leading slash
                    if startsWith(relativePath, '/')
                        newItem = relativePath(2:end);
                    elseif startsWith(relativePath, '\')
                        newItem = relativePath(2:end);
                    else
                        newItem = relativePath;
                    end
                else
                    % If not in Input folder, just use the filename
                    newItem = file;
                end
            else
                % State Properties entries use a custom dialog to collect State = Value parts separately.
                if gui.isStructuredStatePropertyField(fieldPath)
                    [newItem, wasCancelled] = gui.promptStatePropertyItem(fieldPath, '', 'Add List Item');
                    if wasCancelled
                        return;
                    end
                else
                    prompt = {['Enter new item for ', fieldPath, ':']};
                    dlgtitle = 'Add List Item';
                    dims = [1 60];
                    answer = inputdlg(prompt, dlgtitle, dims);
                    if isempty(answer)
                        % Restore window state after dialog cancellation
                        drawnow;
                        pause(0.1);
                        gui.Fig.WindowState = 'maximized';
                        figure(gui.Fig);
                        drawnow;
                        return; % User cancelled
                    end
                    newItem = answer{1};
                end
            end

            if ~isempty(newItem) && ~ismember(newItem, currentItems)
                % For gas fractions, validate the format and check if sum equals 1
                if endsWith(fieldPath, 'gasProperties.fraction') || contains(fieldPath, 'gasProperties.fraction')
                    % Validate format: should be "species = number"
                    pattern = '^\s*(\w+)\s*=\s*([-+]?[0-9]+\.?[0-9]*)\s*$';
                    match = regexp(newItem, pattern, 'tokens');
                    if isempty(match)
                        uialert(gui.Fig, ['Invalid format for gas fraction. Expected format: "species = number"\n' ...
                            'Example: "N2 = 1" or "H2 = 0.5"'], 'Invalid Format');
                        return;
                    end
                    
                    fractionValue = str2double(match{1}{2});
                    if isnan(fractionValue) || fractionValue < 0 || fractionValue > 1
                        uialert(gui.Fig, 'Gas fraction must be between 0 and 1.', 'Invalid Gas Fraction');
                        return;
                    end
                    % Extract all numbers from current items and new item
                    allItems = [currentItems(:); {newItem}];
                    totalSum = 0;
                    for i = 1:length(allItems)
                        matchTokens = regexp(allItems{i}, pattern, 'tokens');
                        if ~isempty(matchTokens)
                            totalSum = totalSum + str2double(matchTokens{1}{2});
                        end
                    end
                    
                end
                
                % Ensure currentItems is a column cell array
                if isempty(currentItems)
                    newItems = {newItem};
                else
                    % Convert to column cell array to ensure consistent concatenation
                    currentItemsCol = currentItems(:);
                    newItemCol = {newItem};
                    newItems = [currentItemsCol; newItemCol];
                end
                listBox.Items = newItems;
                % Update the Setup struct as well
                gui.setNestedField(fieldPath, newItems);
            end
            
            % Update run button state if LXCat files were modified
            if strcmp(fieldPath, 'electronKinetics.LXCatFiles')
                gui.updateRunButtonState();
            end
         end

        function removeListItem(gui, fieldPath)
            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(fieldPath, 'gasProperties.')
                fieldPath = strrep(fieldPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end
            
            % For UIControls access, remove electronKinetics prefix if present
            uiPath = fieldPath;
            if startsWith(uiPath, 'electronKinetics.gasProperties.')
                uiPath = strrep(uiPath, 'electronKinetics.gasProperties.', 'gasProperties.');
            end
            
            listBox = gui.getNestedField(gui.UIControls, uiPath);
            
            % Check if listBox is empty (path not found)
            if isempty(listBox) || ~isobject(listBox)
                uialert(gui.Fig, 'Could not find the list control.', 'Error');
                return;
            end
            
            selectedValues = listBox.Value; % Selected items (cell array when Multiselect is 'on')

            if isempty(selectedValues)
                uialert(gui.Fig, 'No item selected to remove.', 'Selection Error');
                return;
            end

            currentItems = listBox.Items;
            
            % Convert selected values (cell array) to indices
            if iscell(selectedValues)
                indicesToRemove = [];
                for i = 1:length(selectedValues)
                    idx = find(strcmp(currentItems, selectedValues{i}));
                    if ~isempty(idx)
                        indicesToRemove = [indicesToRemove; idx];
                    end
                end
            else
                indicesToRemove = selectedValues;
            end
            
            % Remove selected items in reverse order to avoid index shifting
            currentItems(indicesToRemove) = [];
            listBox.Items = currentItems;
            listBox.Value = {}; % Clear selection

            % For gas fractions, check if sum equals 1 after removal
            if endsWith(fieldPath, 'gasProperties.fraction') || contains(fieldPath, 'gasProperties.fraction')
                pattern = '^\s*(\w+)\s*=\s*([0-9]+\.?[0-9]*)\s*$';
                totalSum = 0;
                for i = 1:length(currentItems)
                    matchTokens = regexp(currentItems{i}, pattern, 'tokens');
                    if ~isempty(matchTokens)
                        totalSum = totalSum + str2double(matchTokens{1}{2});
                    end
                end
                
            end
            
            % Update the Setup struct as well
            gui.setNestedField(fieldPath, currentItems);
            
            % Update run button state if LXCat files were modified
            if strcmp(fieldPath, 'electronKinetics.LXCatFiles')
                gui.updateRunButtonState();
            end
        end

        function editListItem(gui, listBox, fieldPath)
            % Edit an existing item in a listbox (called on double-click)
            selectedValues = listBox.Value;
            
            if isempty(selectedValues)
                uialert(gui.Fig, 'Please select an item to edit.', 'Selection Error');
                return;
            end
            
            % For multiselect, only edit the first selected item
            if iscell(selectedValues)
                itemToEdit = selectedValues{1};
            else
                itemToEdit = selectedValues;
            end
            
            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(fieldPath, 'gasProperties.')
                fieldPath = strrep(fieldPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end

            if gui.isStructuredStatePropertyField(fieldPath)
                [newItem, wasCancelled] = gui.promptStatePropertyItem(fieldPath, itemToEdit, 'Edit List Item');
                if wasCancelled
                    return;
                end
            else
                % Show loading spinner briefly - use classic figure to match inputdlg style
                spinnerFig = [];
            try
                screenSize = get(0, 'ScreenSize');
                dlgWidth = 250;
                dlgHeight = 80;
                dlgX = (screenSize(3) - dlgWidth) / 2;
                dlgY = (screenSize(4) - dlgHeight) / 2 + 100; % Pull down more for better vertical centering
                
                spinnerFig = figure('Name', 'Edit List Item', ...
                    'Position', [dlgX, dlgY, dlgWidth, dlgHeight], ...
                    'MenuBar', 'none', 'ToolBar', 'none', ...
                    'NumberTitle', 'off', 'Resize', 'off', ...
                    'Color', [0.94 0.94 0.94], 'Visible', 'on');
                
                uicontrol(spinnerFig, 'Style', 'text', 'String', 'Loading...', ...
                    'Position', [10 40 230 25], 'FontSize', 10, ...
                    'BackgroundColor', [0.94 0.94 0.94], 'HorizontalAlignment', 'center');
                
                spinnerText = uicontrol(spinnerFig, 'Style', 'text', 'String', '●', ...
                    'Position', [10 10 230 25], 'FontSize', 18, ...
                    'BackgroundColor', [0.94 0.94 0.94], 'ForegroundColor', [0.18 0.70 0.25], ...
                    'HorizontalAlignment', 'center');
                
                drawnow;
                spinnerChars = {'●', '◐', '◓', '◑', '◒'};
                for i = 1:5
                    spinnerText.String = spinnerChars{mod(i-1, 5) + 1};
                    drawnow; pause(0.04);
                end
            catch
                spinnerFig = [];
            end
            
            if ~isempty(spinnerFig) && ishandle(spinnerFig)
                close(spinnerFig);
                drawnow;
            end
            
            % Use original inputdlg for safety and correct functionality
            prompt = {['Edit item for ', fieldPath, ':']};
            dlgtitle = 'Edit List Item';
            dims = [1, 50];
            answer = inputdlg(prompt, dlgtitle, dims, {itemToEdit});
            
            if isempty(answer)
                % Restore window state after dialog cancellation
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = 'maximized';
                figure(gui.Fig);
                drawnow;
                return; % User cancelled
            end
            
                newItem = answer{1};
            end
            if isempty(newItem)
                uialert(gui.Fig, 'Item cannot be empty.', 'Invalid Input');
                return;
            end
            
            currentItems = listBox.Items;
            
            % Find the index of the item to replace
            if iscell(currentItems)
                itemIndex = find(strcmp(currentItems, itemToEdit), 1);
            else
                itemIndex = find(strcmp(currentItems, itemToEdit), 1);
            end
            
            if isempty(itemIndex)
                uialert(gui.Fig, 'Could not find the selected item.', 'Error');
                return;
            end
            
            % Validate gas fractions format if applicable
            if endsWith(fieldPath, 'gasProperties.fraction') || contains(fieldPath, 'gasProperties.fraction')
                pattern = '^\s*(\w+)\s*=\s*([-+]?[0-9]+\.?[0-9]*)\s*$';
                match = regexp(newItem, pattern, 'tokens');
                if isempty(match)
                    uialert(gui.Fig, ['Invalid format for gas fraction. Expected format: "species = number"\n' ...
                        'Example: "N2 = 1" or "H2 = 0.5"'], 'Invalid Format');
                    return;
                end

                fractionValue = str2double(match{1}{2});
                if isnan(fractionValue) || fractionValue < 0 || fractionValue > 1
                    uialert(gui.Fig, 'Gas fraction must be between 0 and 1.', 'Invalid Gas Fraction');
                    return;
                end
                
            end
            
            % Update the list and selection
            currentItems{itemIndex} = newItem;
            listBox.Items = currentItems;
            listBox.Value = newItem;
            
            % Update the Setup struct as well
            gui.setNestedField(fieldPath, currentItems);
            
            % Update run button state if LXCat files were modified
            if strcmp(fieldPath, 'electronKinetics.LXCatFiles')
                gui.updateRunButtonState();
            end
        end

        function browseFolder(gui, ~, ~)
            folderPath = uigetdir(pwd, 'Select Output Folder'); % Start in current directory
            if ~isequal(folderPath, 0)
                % Update the edit field and the setup struct
                control = gui.UIControls.output.folder;
                control.Value = folderPath;
                gui.setNestedField('output.folder', folderPath); % Update setup directly
            end
        end

        function browseFile(gui, fieldPath)
            % Map old gasProperties paths to new electronKinetics.gasProperties paths
            if startsWith(fieldPath, 'gasProperties.')
                fieldPath = strrep(fieldPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end
            
            % Start in Input folder for gas property files
            inputFolder = fullfile(pwd, 'Input');
            if ~isfolder(inputFolder)
                inputFolder = pwd;
            end
            
            % Save window state and bring to front
            drawnow;
            oldState = gui.Fig.WindowState;
            drawnow;
            pause(0.05);
            
            [file, path] = uigetfile('*.*', 'Select File', inputFolder); % Allow any file type
            if isequal(file, 0) || isequal(path, 0)
                % Restore window state after dialog cancellation
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = 'maximized';
                figure(gui.Fig);
                drawnow;
                return;
            end
            
            % Restore window state after file dialog
            drawnow;
            pause(0.1);
            gui.Fig.WindowState = oldState; % Restore previous (typically maximized)
            figure(gui.Fig); % Bring to front
            drawnow;
            fullPath = fullfile(path, file);
            
            % Convert to relative path from Input folder for gas properties
            if contains(fieldPath, 'gasProperties')
                inputFolderFull = strrep(fullfile(pwd, 'Input'), '\', '/');
                fullPathNormalized = strrep(fullPath, '\', '/');
                inputFolderNormalized = strrep(inputFolderFull, '\', '/');
                
                if contains(fullPathNormalized, inputFolderNormalized)
                    % Extract the relative path
                    relativePath = erase(fullPathNormalized, inputFolderNormalized);
                    % Remove leading slash
                    if startsWith(relativePath, '/')
                        filePath = relativePath(2:end);
                    elseif startsWith(relativePath, '\')
                        filePath = relativePath(2:end);
                    else
                        filePath = relativePath;
                    end
                else
                    % If not in Input folder, just use the filename
                    filePath = file;
                end
            else
                % For other files, store the full path
                filePath = fullPath;
            end
            
            uiPath = fieldPath;
            if startsWith(uiPath, 'electronKinetics.gasProperties.')
                uiPath = strrep(uiPath, 'electronKinetics.gasProperties.', 'gasProperties.');
            end

            control = gui.getNestedField(gui.UIControls, uiPath);
            control.Value = filePath;
            gui.setNestedField(fieldPath, filePath); % Update setup directly
        end


        function handleIonizationOperatorChange(gui, value)
            % Update Setup and enable/disable growth model depending on operator
            try
                if isfield(gui.UIControls.electronKinetics, 'eedfType') && ...
                        strcmp(gui.UIControls.electronKinetics.eedfType.Value, 'prescribedEedf')
                    value = 'conservative';
                end
            catch
            end
            % Ensure the UI dropdown reflects the chosen value when called programmatically
            try
                if isfield(gui.UIControls.electronKinetics, 'ionizationOperatorType') && isprop(gui.UIControls.electronKinetics.ionizationOperatorType, 'Value')
                    gui.UIControls.electronKinetics.ionizationOperatorType.Value = value;
                end
            catch
            end
            gui.updateField(gui.UIControls.electronKinetics.ionizationOperatorType, 'electronKinetics.ionizationOperatorType', value);

            excitationFrequency = gui.Setup.workingConditions.excitationFrequency;
            hasExcitationFrequency = ~isempty(excitationFrequency) && ...
                isnumeric(excitationFrequency) && any(excitationFrequency(:) ~= 0);

            if hasExcitationFrequency
                gui.setTemporalGrowthModelForExcitationFrequency();
                gui.UIControls.electronKinetics.growthModelType.Enable = 'off';
            elseif strcmp(value, 'conservative')
                gui.UIControls.electronKinetics.growthModelType.Enable = 'off';
            else
                gui.UIControls.electronKinetics.growthModelType.Enable = 'on';
            end
        end

        function handleEECollisionsChange(gui, value)
            % Update Setup
            gui.setNestedField('electronKinetics.includeEECollisions', value);
            
            % When e-e collisions are enabled, unlock Gas pressure and Electron density
            if value
                % Enable the fields (but keep checkboxes as they are - user can still control them)
                gui.UIControls.workingConditions.gasPressure.Enable = 'on';
                gui.UIControls.workingConditions.electronDensity.Enable = 'on';
                % Also check the checkboxes to indicate they're active
                gui.UIControls.workingConditions.gasPressureCheckbox.Value = true;
                gui.UIControls.workingConditions.electronDensityCheckbox.Value = true;
            else
                % When disabled, block the fields and uncheck the checkboxes
                gui.UIControls.workingConditions.gasPressureCheckbox.Value = false;
                gui.UIControls.workingConditions.electronDensityCheckbox.Value = false;
                gui.toggleGasPressureEnable(false);
                gui.toggleElectronDensityEnable(false);
            end
        end

        function updateRunButtonState(gui)
            % Update the state of the "Generate & Run" button based on LXCat files
            % Button is disabled if Electron Kinetics is disabled, or if it is
            % enabled but no LXCat files are present
            try
                reducedFieldIsMissing = false;
                try
                    reducedFieldValue = gui.Setup.workingConditions.reducedField;
                    if strcmpi(gui.Setup.electronKinetics.eedfType, 'boltzmann') && ...
                       (isempty(reducedFieldValue) || ...
                        ((ischar(reducedFieldValue) || isstring(reducedFieldValue)) && isempty(strtrim(char(reducedFieldValue)))))
                        reducedFieldIsMissing = true;
                    end
                catch
                end

                electronTemperatureIsMissing = false;
                try
                    electronTemperatureValue = gui.Setup.workingConditions.electronTemperature;
                    if strcmpi(gui.Setup.electronKinetics.eedfType, 'prescribedEedf') && ...
                       (isempty(electronTemperatureValue) || ...
                        ((ischar(electronTemperatureValue) || isstring(electronTemperatureValue)) && isempty(strtrim(char(electronTemperatureValue)))))
                        electronTemperatureIsMissing = true;
                    end
                catch
                end

                gui.updateReducedFieldMandatoryMessage(reducedFieldIsMissing);
                gui.updateElectronTemperatureMandatoryMessage(electronTemperatureIsMissing);

                if isfield(gui.UIControls, 'saveInputButton')
                    if reducedFieldIsMissing || electronTemperatureIsMissing
                        gui.UIControls.saveInputButton.Enable = 'off';
                    else
                        gui.UIControls.saveInputButton.Enable = 'on';
                    end
                end

                if isfield(gui.UIControls, 'runButton')
                    % If Electron Kinetics is missing or explicitly disabled -> disable run
                    if ~isfield(gui.Setup, 'electronKinetics') || ~gui.Setup.electronKinetics.isOn
                        gui.UIControls.runButton.Enable = 'off';
                    else
                        % When Electron Kinetics is enabled, require at least one LXCat file
                        items = {};
                        try
                            items = gui.UIControls.electronKinetics.LXCatFiles.Items;
                        catch
                        end
                        if isempty(items) || reducedFieldIsMissing || electronTemperatureIsMissing
                            gui.UIControls.runButton.Enable = 'off';
                        else
                            gui.UIControls.runButton.Enable = 'on';
                        end
                    end
                end
            catch
                % Fail safe: enable button if something goes wrong
                try
                    if isfield(gui.UIControls, 'runButton')
                        gui.UIControls.runButton.Enable = 'on';
                    end
                    if isfield(gui.UIControls, 'saveInputButton')
                        gui.UIControls.saveInputButton.Enable = 'on';
                    end
                catch
                end
            end
        end

        function updateReducedFieldMandatoryMessage(gui, showMessage)
            try
                if isfield(gui.UIControls, 'workingConditions') && ...
                   isfield(gui.UIControls.workingConditions, 'reducedFieldMessage')
                    if showMessage
                        gui.UIControls.workingConditions.reducedFieldMessage.Visible = 'on';
                    else
                        gui.UIControls.workingConditions.reducedFieldMessage.Visible = 'off';
                    end
                end

                if isfield(gui.UIControls, 'workingConditions') && ...
                   isfield(gui.UIControls.workingConditions, 'reducedFieldBorderPanel')
                    if showMessage
                        gui.UIControls.workingConditions.reducedFieldBorderPanel.BackgroundColor = [0.75 0.1 0.1];
                    else
                        gui.UIControls.workingConditions.reducedFieldBorderPanel.BackgroundColor = [0.94 0.94 0.94];
                    end
                end
            catch
            end
        end

        function updateElectronTemperatureMandatoryMessage(gui, showMessage)
            try
                if isfield(gui.UIControls, 'workingConditions') && ...
                   isfield(gui.UIControls.workingConditions, 'electronTemperatureMessage')
                    if showMessage
                        gui.UIControls.workingConditions.electronTemperatureMessage.Visible = 'on';
                    else
                        gui.UIControls.workingConditions.electronTemperatureMessage.Visible = 'off';
                    end
                end

                if isfield(gui.UIControls, 'workingConditions') && ...
                   isfield(gui.UIControls.workingConditions, 'electronTemperatureBorderPanel')
                    if showMessage
                        gui.UIControls.workingConditions.electronTemperatureBorderPanel.BackgroundColor = [0.75 0.1 0.1];
                    else
                        gui.UIControls.workingConditions.electronTemperatureBorderPanel.BackgroundColor = [0.94 0.94 0.94];
                    end
                end
            catch
            end
        end

        function handleTabChange(gui, evt)
            % Prevent leaving Electron Kinetics tab if LXCat files list is empty
            try
                oldTab = evt.OldValue;
                newTab = evt.NewValue;

                % If Electron Kinetics is disabled, prevent entering several tabs
                try
                    if isfield(gui.Setup, 'electronKinetics') && ~gui.Setup.electronKinetics.isOn
                        blockedTitles = {'Gas Properties', 'State Properties', 'Numerics', 'Output'};
                        try
                            if isprop(newTab, 'Title') && any(strcmp(newTab.Title, blockedTitles))
                                evt.Source.SelectedTab = gui.UIControls.tabs.kineticsTab;
                                uialert(gui.Fig, sprintf('%s is unavailable when Electron Kinetics is disabled.', newTab.Title), 'Tab Disabled');
                                return;
                            end
                        catch
                        end
                    end
                catch
                end

                % Only enforce when leaving the Electron Kinetics tab
                if isequal(oldTab, gui.UIControls.tabs.kineticsTab) && ~isequal(newTab, gui.UIControls.tabs.kineticsTab)
                    % Only enforce if Electron Kinetics is enabled
                    if gui.Setup.electronKinetics.isOn
                        items = gui.UIControls.electronKinetics.LXCatFiles.Items;
                        if isempty(items)
                            % Revert selection and warn user
                            evt.Source.SelectedTab = gui.UIControls.tabs.kineticsTab;
                            uialert(gui.Fig, 'At least one LXCat file is required. Please add a file before leaving the Electron Kinetics tab.', ...
                                'Missing LXCat Files');
                        end
                    end
                end
                
                % Update run button state after tab change
                gui.updateRunButtonState();
            catch
                % Fail safe: do not block tab change if something goes wrong
            end
        end

        function toggleGuiControls(gui, isEnabled)
            % Enable/disable GUI controls based on checkbox state
            if isEnabled
                gui.UIControls.gui.refreshFrequency.Enable = 'on';
            else
                gui.UIControls.gui.refreshFrequency.Enable = 'off';
            end
            
            % Update the setup struct
            gui.setNestedField('gui.isOn', isEnabled);
        end

        function toggleOutputEnable(gui, isEnabled)
            % Enable/disable output-related controls based on checkbox state
            
            % Update Setup flag
            gui.setNestedField('output.isOn', isEnabled);

            if isEnabled
                enableState = 'on';
            else
                enableState = 'off';
            end

            if isfield(gui.UIControls, 'output')
                if isfield(gui.UIControls.output, 'dataFormat') && isvalid(gui.UIControls.output.dataFormat)
                    gui.UIControls.output.dataFormat.Enable = enableState;
                end
                if isfield(gui.UIControls.output, 'useCustomFolder') && isvalid(gui.UIControls.output.useCustomFolder)
                    gui.UIControls.output.useCustomFolder.Enable = enableState;
                end
                gui.toggleOutputFolderMode(isEnabled && ...
                    (~isfield(gui.UIControls.output, 'useCustomFolder') || ...
                    gui.UIControls.output.useCustomFolder.Value));
            end

            if isfield(gui.UIControls, 'output') && ...
                    isfield(gui.UIControls.output, 'dataSets')
                dataSetNames = fieldnames(gui.UIControls.output.dataSets);
                for i = 1:length(dataSetNames)
                    checkbox = gui.UIControls.output.dataSets.(dataSetNames{i});
                    if isvalid(checkbox)
                        checkbox.Enable = enableState;
                        if ~isEnabled
                            checkbox.Value = false;
                        end
                    end
                end
                if ~isEnabled
                    gui.Setup.output.dataSets = {};
                end
            end
        end

        function toggleOutputFolderMode(gui, useCustomFolder)
            % Enable/disable manual output folder controls.
            if useCustomFolder
                enableState = 'on';
            else
                enableState = 'off';
            end

            if isfield(gui.UIControls, 'output')
                if isfield(gui.UIControls.output, 'folder') && isvalid(gui.UIControls.output.folder)
                    gui.UIControls.output.folder.Enable = enableState;
                end
                if isfield(gui.UIControls.output, 'browseButton') && isvalid(gui.UIControls.output.browseButton)
                    gui.UIControls.output.browseButton.Enable = enableState;
                end
            end
        end

        function toggleCARGasEnable(gui, isEnabled)
            % Enable/disable CAR gas controls based on checkbox state
            if isEnabled
                gui.UIControls.electronKinetics.CARgases.Enable = 'on';
            else
                gui.UIControls.electronKinetics.CARgases.Enable = 'off';
            end
            % Enable/disable add/remove buttons using stored references
            gui.UIControls.electronKinetics.CARaddButton.Enable = isEnabled;
            gui.UIControls.electronKinetics.CARremoveButton.Enable = isEnabled;
        end

        function toggleDischargeCurrentEnable(gui, isEnabled)
            % Enable/disable discharge current field based on checkbox state
            if isEnabled
                gui.UIControls.workingConditions.dischargeCurrent.Enable = 'on';
            else
                gui.UIControls.workingConditions.dischargeCurrent.Enable = 'off';
            end
        end

        function toggleDischargePowerEnable(gui, isEnabled)
            % Enable/disable discharge power density field based on checkbox state
            if isEnabled
                gui.UIControls.workingConditions.dischargePowerDensity.Enable = 'on';
            else
                gui.UIControls.workingConditions.dischargePowerDensity.Enable = 'off';
            end
        end

        function toggleGasPressureEnable(gui, isEnabled)
            % Enable/disable gas pressure based on checkbox state
            if isEnabled
                gui.UIControls.workingConditions.gasPressure.Enable = 'on';
            else
                gui.UIControls.workingConditions.gasPressure.Enable = 'off';
            end
        end

        function toggleElectronDensityEnable(gui, isEnabled)
            % Enable/disable electron density based on checkbox state
            if isEnabled
                gui.UIControls.workingConditions.electronDensity.Enable = 'on';
            else
                gui.UIControls.workingConditions.electronDensity.Enable = 'off';
            end
        end


        function handleUpdateFactorChange(gui, value)
            % Validate Update Factor: cannot be 0
            if value == 0
                uialert(gui.Fig, 'Update Factor cannot be 0. Using default value 0.05.', ...
                    'Invalid Value', 'Icon', 'error');
                gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Value = 0.05;
                gui.setNestedField('electronKinetics.numerics.energyGrid.smartGrid.updateFactor', 0.05);
            else
                gui.setNestedField('electronKinetics.numerics.energyGrid.smartGrid.updateFactor', value);
            end
        end

        function toggleElectronKineticsEnable(gui, isEnabled, resetToDefaults)
            % Enable/disable electron kinetics controls based on master checkbox
            if nargin < 3
                resetToDefaults = false;
            end

            if isEnabled && resetToDefaults && ~gui.IsApplyingDefaultSetup
                gui.IsApplyingDefaultSetup = true;
                try
                    gui.initializeDefaultSetup();
                    gui.populateGUIFromSetup();
                    gui.IsApplyingDefaultSetup = false;
                    return;
                catch ME
                    gui.IsApplyingDefaultSetup = false;
                    rethrow(ME);
                end
            end

            % Update Setup flag
            gui.Setup.electronKinetics.isOn = isEnabled;

            % Core controls that always follow the master flag
            coreControls = {
                'eedfType', ...
                'shapeParameter', ...
                'shapeParameterPrecise', ...
                'ionizationOperatorType', ...
                'growthModelType', ...
                'includeEECollisions', ...
                'LXCatFiles' ...
            };

            if isEnabled
                state = 'on';
            else
                state = 'off';
            end

            % Apply state to core controls
            for i = 1:numel(coreControls)
                try
                    ctrl = gui.UIControls.electronKinetics.(coreControls{i});
                    ctrl.Enable = state;
                catch
                    % Ignore missing or non-enable-able controls
                end
            end
            
            % Update run button state when Electron Kinetics is toggled
            gui.updateRunButtonState();

            % If electronKinetics is disabled, LXCat cross sections should not be added or removed
            try
                % LXCat Files
                if isfield(gui.UIControls.electronKinetics, 'LXCatFiles')
                    gui.UIControls.electronKinetics.LXCatFiles.Enable = state;
                    % Enable/disable add/remove buttons (they're at columns 3 and 4, row 1)
                    parent = gui.UIControls.electronKinetics.LXCatFiles.Parent;
                    for i = 1:length(parent.Children)
                        child = parent.Children(i);
                        if strcmp(child.Type, 'uibutton') && child.Layout.Row == 1 && (child.Layout.Column == 3 || child.Layout.Column == 4)
                            child.Enable = state;
                        end
                    end
                end

            catch
            end

            % Optional groups (extra LXCat, effective pops, CAR gases) depend on both
            % the master flag AND their individual checkboxes
            try
                % LXCat Extra Files
                if isfield(gui.UIControls.electronKinetics, 'LXCatExtraCheckbox')
                    gui.UIControls.electronKinetics.LXCatExtraCheckbox.Enable = state;
                    if isEnabled
                        gui.toggleLXCatExtraEnable(gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value);
                    else
                        gui.toggleLXCatExtraEnable(false);
                    end
                end
            catch
            end

            try
                % Effective Cross Section Populations
                if isfield(gui.UIControls.electronKinetics, 'effectivePopCheckbox')
                    gui.UIControls.electronKinetics.effectivePopCheckbox.Enable = state;
                    if isEnabled
                        gui.toggleEffectivePopEnable(gui.UIControls.electronKinetics.effectivePopCheckbox.Value);
                    else
                        gui.toggleEffectivePopEnable(false);
                    end
                end
            catch
            end

            try
                % CAR Gases
                if isfield(gui.UIControls.electronKinetics, 'CARcheckbox')
                    gui.UIControls.electronKinetics.CARcheckbox.Enable = state;
                    if isEnabled
                        gui.toggleCARGasEnable(gui.UIControls.electronKinetics.CARcheckbox.Value);
                    else
                        gui.toggleCARGasEnable(false);
                    end
                end
            catch
            end

            % If Enable Electron Kinetics is turned off and then on again, reset all of the "General Kinetics Settings" 
            % fields to their default values
            % if isEnabled
                % gui.Setup.electronKinetics.shapeParameter = 1; % Only used for prescribedEedf, but reset it anyway
                % gui.handleEedfTypeChange('boltzmann'); % This will also update the setup field for eedfType

                % gui.Setup.electronKinetics.ionizationOperatorType = 'usingSDCS'; % Default to using SDCS
                % gui.handleIonizationOperatorChange(gui.Setup.electronKinetics.ionizationOperatorType); % This will also enable/disable growth model
                % gui.UIControls.electronKinetics.growthModelType.Value = 'temporal';
                % gui.Setup.electronKinetics.growthModelType = 'temporal';

                % gui.UIControls.electronKinetics.includeEECollisions.Value = false; % Default to no e-e collisions
                % gui.handleEECollisionsChange(false); % This will also enable/disable related fields
            % end
        end

        function toggleLXCatExtraEnable(gui, isEnabled)
            % Enable/disable LXCat extra files controls based on checkbox state
            gui.UIControls.electronKinetics.LXCatExtraFiles.Enable = isEnabled;
            % Enable/disable add/remove buttons (they're at columns 3 and 4, row 2)
            parent = gui.UIControls.electronKinetics.LXCatExtraFiles.Parent;
            for i = 1:length(parent.Children)
                child = parent.Children(i);
                if strcmp(child.Type, 'uibutton') && child.Layout.Row == 2 && (child.Layout.Column == 3 || child.Layout.Column == 4)
                    child.Enable = isEnabled;
                end
            end
        end

        function toggleEffectivePopEnable(gui, isEnabled)
            % Enable/disable effective cross section populations based on checkbox state
            gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Enable = isEnabled;
            % Enable/disable add/remove buttons (they're at columns 3 and 4, row 3)
            parent = gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Parent;
            for i = 1:length(parent.Children)
                child = parent.Children(i);
                if strcmp(child.Type, 'uibutton') && child.Layout.Row == 3 && (child.Layout.Column == 3 || child.Layout.Column == 4)
                    child.Enable = isEnabled;
                end
            end
        end

        function handleEedfTypeChange(gui, eedfType)
            % Show/hide shape parameter based on EEDF type
            % Ensure the UI dropdown reflects the chosen type when called programmatically
            try
                if isfield(gui.UIControls.electronKinetics, 'eedfType') && isprop(gui.UIControls.electronKinetics.eedfType, 'Value')
                    gui.UIControls.electronKinetics.eedfType.Value = eedfType;
                end
            catch
            end
            if strcmp(eedfType, 'prescribedEedf')
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'on';
                gui.UIControls.electronKinetics.shapeParameterPrecise.Visible = 'on';
                gui.UIControls.electronKinetics.shapeParameterLabel.Visible = 'on';

                gui.UIControls.electronKinetics.includeEECollisions.Value = false;
                gui.UIControls.electronKinetics.includeEECollisions.Enable = 'off';
                gui.handleIonizationOperatorChange('conservative');

                % Update the label visibility as well
                parent = gui.UIControls.electronKinetics.shapeParameter.Parent;
                children = parent.Children;
                for i = 1:length(children)
                    if strcmp(children(i).Type, 'uilabel') && children(i).Layout.Row == 3 && children(i).Layout.Column == 1
                        children(i).Visible = 'on';
                    end
                end
            else
                gui.Setup.electronKinetics.shapeParameter = 1; % Reset to default when not using prescribed EEDF
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'off';
                gui.UIControls.electronKinetics.shapeParameterPrecise.Visible = 'off';
                gui.UIControls.electronKinetics.shapeParameterLabel.Visible = 'off';

                gui.UIControls.electronKinetics.includeEECollisions.Enable = 'on';

                % Directly update the UI controls so the change is visible
                try
                    if isfield(gui.UIControls.electronKinetics, 'shapeParameter') && isprop(gui.UIControls.electronKinetics.shapeParameter, 'Value')
                        gui.UIControls.electronKinetics.shapeParameter.Value = 1;
                    end
                    if isfield(gui.UIControls.electronKinetics, 'shapeParameterPrecise') && isprop(gui.UIControls.electronKinetics.shapeParameterPrecise, 'Value')
                        gui.UIControls.electronKinetics.shapeParameterPrecise.Value = 1;
                    end
                catch
                end
                % Also update the underlying Setup field
                gui.setNestedField('electronKinetics.shapeParameter', 1);
                drawnow;
                
                % Update the label visibility as well
                parent = gui.UIControls.electronKinetics.shapeParameter.Parent;
                children = parent.Children;
                for i = 1:length(children)
                    if strcmp(children(i).Type, 'uilabel') && children(i).Layout.Row == 3 && children(i).Layout.Column == 1
                        children(i).Visible = 'off';
                    end
                end
            end
            % Also update the setup
            gui.setNestedField('electronKinetics.eedfType', eedfType);
            gui.updateRunButtonState();
        end

        function handleEedfTypeChangeForInit(gui, eedfType)
            % Show/hide shape parameter based on EEDF type (for initialization)
            if strcmp(eedfType, 'prescribedEedf')
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'on';
                gui.UIControls.electronKinetics.shapeParameterPrecise.Visible = 'on';
                gui.UIControls.electronKinetics.shapeParameterLabel.Visible = 'on';
                gui.handleIonizationOperatorChange('conservative');

                % Update the label visibility as well
                parent = gui.UIControls.electronKinetics.shapeParameter.Parent;
                children = parent.Children;
                for i = 1:length(children)
                    if strcmp(children(i).Type, 'uilabel') && children(i).Layout.Row == 3 && children(i).Layout.Column == 1
                        children(i).Visible = 'on';
                    end
                end
            else
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'off';
                gui.UIControls.electronKinetics.shapeParameterPrecise.Visible = 'off';
                gui.UIControls.electronKinetics.shapeParameterLabel.Visible = 'off';

                % Update the label visibility as well
                parent = gui.UIControls.electronKinetics.shapeParameter.Parent;
                children = parent.Children;
                for i = 1:length(children)
                    if strcmp(children(i).Type, 'uilabel') && children(i).Layout.Row == 3 && children(i).Layout.Column == 1
                        children(i).Visible = 'off';
                    end
                end
            end
        end

        function clearAllFields(gui, ~, ~)
            % Clear user-entered setup fields while keeping dropdowns at defaults.
            currentSetup = gui.Setup;
            gui.initializeDefaultSetup();
            defaultSetup = gui.Setup;
            gui.Setup = currentSetup;

            gui.clearControlTree(gui.UIControls, '', defaultSetup);
            gui.applyClearedControlStates();
            gui.returnToElectronKineticsTab();
            gui.updateRunButtonState();
        end

        function returnToElectronKineticsTab(gui)
            try
                if isfield(gui.UIControls, 'tabs') && ...
                        isfield(gui.UIControls.tabs, 'tabGroup') && ...
                        isfield(gui.UIControls.tabs, 'kineticsTab') && ...
                        isvalid(gui.UIControls.tabs.tabGroup) && ...
                        isvalid(gui.UIControls.tabs.kineticsTab)
                    gui.UIControls.tabs.tabGroup.SelectedTab = gui.UIControls.tabs.kineticsTab;
                end
            catch ME
                warning('Could not return to Electron Kinetics tab: %s', ME.message);
            end
        end

        function applyClearedControlStates(gui)
            try gui.toggleGuiControls(false); catch, end
            try gui.toggleGasPressureEnable(false); catch, end
            try gui.toggleElectronDensityEnable(false); catch, end
            try gui.toggleDischargeCurrentEnable(false); catch, end
            try gui.toggleDischargePowerEnable(false); catch, end
            try gui.toggleOutputEnable(false); catch, end
            try gui.toggleSmartGridEnable(false); catch, end
            try gui.toggleOdeParametersEnable(false); catch, end
            try gui.toggleLXCatExtraEnable(false); catch, end
            try gui.toggleEffectivePopEnable(false); catch, end
            try gui.toggleCARGasEnable(false); catch, end
            try gui.handleEedfTypeChange('boltzmann'); catch, end
            try gui.toggleElectronKineticsEnable(false); catch, end
        end

        function clearControlTree(gui, controls, prefix, defaultSetup)
            if ~isstruct(controls)
                return;
            end

            controlNames = fieldnames(controls);
            for i = 1:numel(controlNames)
                controlName = controlNames{i};
                control = controls.(controlName);

                if isempty(prefix)
                    uiPath = controlName;
                else
                    uiPath = [prefix, '.', controlName];
                end

                if isstruct(control)
                    gui.clearControlTree(control, uiPath, defaultSetup);
                else
                    gui.clearSingleControl(control, uiPath, defaultSetup);
                end
            end
        end

        function clearSingleControl(gui, control, uiPath, defaultSetup)
            try
                if isempty(control) || ~isvalid(control)
                    return;
                end
            catch
                return;
            end

            setupPath = gui.uiPathToSetupPath(uiPath);

            try
                if isa(control, 'matlab.ui.control.NumericEditField')
                    if strcmp(setupPath, 'electronKinetics.shapeParameter') || ...
                            strcmp(uiPath, 'electronKinetics.shapeParameterPrecise')
                        control.Value = 1;
                        gui.clearSetupField('electronKinetics.shapeParameter', 1);
                    else
                        try
                            control.Value = [];
                            gui.clearSetupField(setupPath, []);
                        catch
                            fallbackValue = gui.getNumericClearFallback(control);
                            control.Value = fallbackValue;
                            gui.clearSetupField(setupPath, fallbackValue);
                        end
                    end
                elseif isa(control, 'matlab.ui.control.EditField')
                    control.Value = '';
                    gui.clearSetupField(setupPath, '');
                elseif isa(control, 'matlab.ui.control.ListBox')
                    control.Items = {};
                    control.Value = {};
                    gui.clearSetupField(setupPath, {});
                elseif isa(control, 'matlab.ui.control.CheckBox')
                    control.Value = false;
                    gui.clearSetupField(setupPath, false);
                elseif isa(control, 'matlab.ui.control.DropDown')
                    [hasDefault, defaultValue] = gui.tryGetNestedField(defaultSetup, setupPath);
                    if hasDefault && any(strcmp(control.Items, defaultValue))
                        control.Value = defaultValue;
                        gui.clearSetupField(setupPath, defaultValue);
                    elseif ~isempty(control.Items)
                        control.Value = control.Items{1};
                        gui.clearSetupField(setupPath, control.Value);
                    end
                end
            catch ME
                warning('Could not clear control %s: %s', uiPath, ME.message);
            end
        end

        function setupPath = uiPathToSetupPath(~, uiPath)
            setupPath = uiPath;
            if startsWith(setupPath, 'gasProperties.')
                setupPath = strrep(setupPath, 'gasProperties.', 'electronKinetics.gasProperties.');
            end
        end

        function fallbackValue = getNumericClearFallback(~, control)
            fallbackValue = 0;
            try
                limits = control.Limits;
                if numel(limits) == 2 && (fallbackValue < limits(1) || fallbackValue > limits(2))
                    if isfinite(limits(1))
                        fallbackValue = limits(1);
                    elseif isfinite(limits(2))
                        fallbackValue = limits(2);
                    end
                end
            catch
            end
        end

        function clearSetupField(gui, setupPath, value)
            if isempty(setupPath)
                return;
            end

            [fieldExists, ~] = gui.tryGetNestedField(gui.Setup, setupPath);
            if fieldExists
                gui.setNestedField(setupPath, value);
            end
        end

        function [fieldExists, value] = tryGetNestedField(~, startStruct, fieldPath)
            fieldExists = false;
            value = [];
            if isempty(fieldPath) || ~isstruct(startStruct)
                return;
            end

            parts = strsplit(fieldPath, '.');
            currentValue = startStruct;
            for i = 1:numel(parts)
                if ~isstruct(currentValue) || ~isfield(currentValue, parts{i})
                    return;
                end
                currentValue = currentValue.(parts{i});
            end

            fieldExists = true;
            value = currentValue;
        end

        function loadSettings(gui, ~, ~)
            % Load settings from a .in or .json input file
            inputFolder = fullfile(pwd, 'Input');
            if ~isfolder(inputFolder)
                inputFolder = pwd;
            end

            % Save window state before opening dialog
            drawnow;
            oldState = gui.Fig.WindowState;
            drawnow;
            pause(0.1);
            
            % Allow both .in and .json files
            [file, path] = uigetfile({'*.in;*.json', 'LoKI Input Files (*.in, *.json)'; '*.in', 'LoKI Input Files (*.in)'; '*.json', 'LoKI JSON Files (*.json)'}, 'Load LoKI Input File', inputFolder);
            if isequal(file, 0) || isequal(path, 0)
                % Restore window state after dialog cancellation
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = oldState;
                figure(gui.Fig);
                drawnow;
                return; % User cancelled
            end
            
            % Restore window state immediately after file selection
            drawnow;
            pause(0.1);
            gui.Fig.WindowState = oldState;
            figure(gui.Fig);
            drawnow;
            
            inputFile = fullfile(path, file);
            [~, ~, ext] = fileparts(file);
            
            try
                % Reset to defaults first to clear any unused fields
                gui.initializeDefaultSetup();
                
                % Parse based on file extension
                if strcmpi(ext, '.json')
                    % Load from JSON file
                    [gui.Setup, foundFields] = gui.parseJSONFile(inputFile);
                else
                    % Load from .in file (default)
                    [gui.Setup, foundFields] = gui.parseInputFile(inputFile);
                end
                
                gui.populateGUIFromSetup(); % Update UI
                % Disable optional fields that were not found in the file
                gui.disableUnusedOptionalFields(foundFields);
                % Activate and expand fields that were found in the file
                gui.activateFoundOptionalFields(foundFields);
                % Ensure window stays maximized
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = oldState;
                figure(gui.Fig);
                drawnow;
                gui.handleFigureSizeChanged();
                uialert(gui.Fig, ['Settings loaded from ', file], 'Load Successful', 'Icon', 'success');
            catch ME
                % Ensure window stays maximized even on error
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = oldState;
                figure(gui.Fig);
                drawnow;
                uialert(gui.Fig, ['Error loading file: ', ME.message], 'Load Error');
            end
        end

        function [setup, foundFields] = parseInputFile(gui, inputFile)
            % Parse a LoKI-B .in file written by generateInputFile (YAML-like subset).
            %
            % Strategy:
            % - Start from the current gui.Setup (assumed defaults already initialized)
            % - Override fields found in the file
            % - Track which optional fields were found
            txt = fileread(inputFile);
            lines = regexp(txt, '\r\n|\n|\r', 'split');

            setup = gui.Setup;
            foundFields = containers.Map('KeyType', 'char', 'ValueType', 'logical'); % Track found fields

            keyStack = {};
            indentStack = [-1];

            i = 1;
            while i <= numel(lines)
                raw = lines{i};
                if isempty(raw)
                    i = i + 1;
                    continue;
                end
                % Strip full-line comments (but keep leading indentation for structure)
                if startsWith(strtrim(raw), '%')
                    i = i + 1;
                    continue;
                end
                if startsWith(strtrim(raw), '#')
                    i = i + 1;
                    continue;
                end

                % Strip inline comments (LoKI-B input files commonly use "%" at end of line)
                % Important: do this AFTER skipping full-line comments.
                pctIdx = strfind(raw, '%');
                if ~isempty(pctIdx)
                    raw = raw(1:pctIdx(1)-1);
                end
                hashIdx = strfind(raw, '#');
                if ~isempty(hashIdx)
                    raw = raw(1:hashIdx(1)-1);
                end

                indent = numel(regexp(raw, '^\s*', 'match', 'once'));
                line = strtrim(raw);
                if isempty(line)
                    i = i + 1;
                    continue;
                end

                % Pop stack if indentation decreased
                while ~isempty(indentStack) && indent <= indentStack(end) && numel(indentStack) > 1
                    indentStack(end) = [];
                    keyStack(end) = [];
                end

                % List item: "- something"
                if startsWith(line, '-')
                    % Append to last list key
                    item = strtrim(line(2:end));
                    if isempty(keyStack)
                        i = i + 1;
                        continue;
                    end
                    pathParts = keyStack;
                    % Get current list (should be empty or already started from "key:" line)
                    currentList = getByPath(setup, pathParts);
                    if isempty(currentList) || ~iscell(currentList)
                        currentList = {};
                    end
                    % Append item to list (completely replacing, not adding to defaults)
                    currentList{end+1} = item;
                    setup = setByPath(setup, pathParts, currentList);
                    % Mark this list field as found
                    fieldPath = strjoin(pathParts, '.');
                    foundFields(fieldPath) = true;
                    i = i + 1;
                    continue;
                end

                % Key/value line: "key: value" or "key:"
                tok = regexp(line, '^([^:]+):\s*(.*)$', 'tokens', 'once');
                if isempty(tok)
                    i = i + 1;
                    continue;
                end
                key = strtrim(tok{1});
                valStr = strtrim(tok{2}); % Trim value string to remove trailing spaces from inline comments

                if isempty(valStr)
                    % Lookahead to decide struct vs list
                    j = i + 1;
                    while j <= numel(lines) && isempty(strtrim(lines{j}))
                        j = j + 1;
                    end
                    isList = false;
                    if j <= numel(lines)
                        nxtRaw = lines{j};
                        nxtIndent = numel(regexp(nxtRaw, '^\s*', 'match', 'once'));
                        nxtLine = strtrim(nxtRaw);
                        if nxtIndent > indent && startsWith(nxtLine, '-')
                            isList = true;
                        end
                    end

                    if isList
                        % Initialize list as empty (will be populated by "- " lines)
                        pathParts = [keyStack, {key}];
                        setup = setByPath(setup, pathParts, {}); % Start with empty list
                        % Mark this list field as found
                        fieldPath = strjoin(pathParts, '.');
                        foundFields(fieldPath) = true;
                        % Push key so "- " lines append to it
                        keyStack = pathParts;
                        indentStack(end+1) = indent; %#ok<AGROW>
                    else
                        % Struct header
                        keyStack{end+1} = key; %#ok<AGROW>
                        indentStack(end+1) = indent; %#ok<AGROW>
                        % Ensure struct exists
                        pathParts = keyStack;
                        cur = getByPath(setup, pathParts);
                        if isempty(cur)
                            setup = setByPath(setup, pathParts, struct());
                        end
                        % Mark this struct field as found (even if empty, it was present in file)
                        fieldPath = strjoin(pathParts, '.');
                        foundFields(fieldPath) = true;
                    end
                else
                    pathParts = [keyStack, {key}];
                    setup = setByPath(setup, pathParts, parseScalar(valStr));
                    % Mark this field as found
                    fieldPath = strjoin(pathParts, '.');
                    foundFields(fieldPath) = true;
                end

                i = i + 1;
            end

            % Local helpers
            function v = parseScalar(s)
                s = strtrim(s);
                % Booleans
                if strcmpi(s, 'true')
                    v = true; return;
                elseif strcmpi(s, 'false')
                    v = false; return;
                end
                % Numeric
                n = str2double(s);
                if ~isnan(n) && ~contains(s, ' ') && ~contains(lower(s), 'logspace') && ~contains(lower(s), 'linspace')
                    v = n; return;
                end
                % Keep strings as-is (expressions like logspace(...) are expected)
                v = s;
            end

            function cur = getByPath(s, parts)
                cur = s;
                for k = 1:numel(parts)
                    p = parts{k};
                    if ~isstruct(cur) || ~isfield(cur, p)
                        cur = [];
                        return;
                    end
                    cur = cur.(p);
                end
            end

            function s = setByPath(s, parts, value)
                if isempty(parts)
                    s = value;
                    return;
                end
                p = parts{1};
                if numel(parts) == 1
                    s.(p) = value;
                    return;
                end
                if ~isfield(s, p) || ~isstruct(s.(p))
                    s.(p) = struct();
                end
                s.(p) = setByPath(s.(p), parts(2:end), value);
            end

            function out = appendList(cur, item)
                if isempty(cur)
                    out = {item};
                elseif iscell(cur)
                    out = cur;
                    out{end+1} = item;
                else
                    out = {item};
                end
            end
        end

        function [setup, foundFields] = parseJSONFile(gui, jsonFile)
            % Parse a LoKI-B .json file into a MATLAB struct.
            % Validates the JSON format and converts it to the same structure as parseInputFile.
            
            % Read and parse JSON
            try
                jsonText = fileread(jsonFile);
                jsonData = jsondecode(jsonText);
            catch ME
                error('Invalid JSON file format: %s', ME.message);
            end
            
            % Validate required top-level sections
            requiredSections = {'workingConditions', 'electronKinetics'};
            for i = 1:length(requiredSections)
                if ~isfield(jsonData, requiredSections{i})
                    error('Missing required section: %s', requiredSections{i});
                end
            end
            
            % Initialize with defaults
            setup = gui.Setup;
            foundFields = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            
            % Convert JSON structure to Setup structure, tracking found fields
            setup = gui.convertJSONToSetup(jsonData, setup, foundFields, '');
        end
        
        function setup = convertJSONToSetup(gui, jsonData, setup, foundFields, prefix)
            % Recursively convert JSON structure to Setup structure
            % prefix: current path prefix (e.g., 'workingConditions.')
            
            if ~isstruct(jsonData)
                return;
            end
            
            fields = fieldnames(jsonData);
            
            for i = 1:length(fields)
                fieldName = fields{i};
                value = jsonData.(fieldName);
                
                % Build current path
                if isempty(prefix)
                    currentPath = fieldName;
                else
                    currentPath = [prefix, '.', fieldName];
                end
                
                % Mark field as found
                foundFields(currentPath) = true;
                
                if isstruct(value)
                    % Nested structure - recurse
                    if ~isfield(setup, fieldName) || ~isstruct(setup.(fieldName))
                        setup.(fieldName) = struct();
                    end
                    setup.(fieldName) = gui.convertJSONToSetup(value, setup.(fieldName), foundFields, currentPath);
                elseif iscell(value)
                    % Cell array (list) - convert to cell array of strings
                    setup.(fieldName) = cellfun(@(x) gui.jsonValueToString(x), value, 'UniformOutput', false);
                elseif isnumeric(value) || islogical(value)
                    % Numeric or boolean - keep as is
                    setup.(fieldName) = value;
                else
                    % String - keep as string
                    if ischar(value) || isstring(value)
                        setup.(fieldName) = char(value);
                    else
                        setup.(fieldName) = value;
                    end
                end
            end
        end
        
        function str = jsonValueToString(gui, value)
            % Convert JSON value to string representation
            if isnumeric(value)
                str = num2str(value);
            elseif islogical(value)
                if value
                    str = 'true';
                else
                    str = 'false';
                end
            else
                str = char(value);
            end
        end
        
        function generateJSONFile(gui, filename)
            % Generates a JSON input file from the gui.Setup struct
            
            
            % Convert Setup struct to JSON structure
            jsonStruct = gui.convertSetupToJSON(gui.getSetupForExport());
            
            % Write JSON file with pretty formatting
            try
                jsonText = jsonencode(jsonStruct, 'PrettyPrint', true);
                fid = fopen(filename, 'w');
                if fid == -1
                    error('Cannot open file "%s" for writing.', filename);
                end
                fprintf(fid, '%s', jsonText);
                fclose(fid);
                fprintf('JSON file generated: %s\n', filename);
            catch ME
                error('Error writing JSON file: %s', ME.message);
            end
        end

        function setup = getSetupForExport(gui)
            % Return setup data with export-only UI choices applied.
            setup = gui.Setup;
            if isfield(setup, 'workingConditions')
                fieldsToRemoveIfBlank = {'reducedField', 'electronTemperature'};
                for idx = 1:length(fieldsToRemoveIfBlank)
                    fieldName = fieldsToRemoveIfBlank{idx};
                    if isfield(setup.workingConditions, fieldName)
                        value = setup.workingConditions.(fieldName);
                        if isempty(value) || ...
                           ((ischar(value) || isstring(value)) && isempty(strtrim(char(value))))
                            setup.workingConditions = rmfield(setup.workingConditions, fieldName);
                        end
                    end
                end
            end
            if isfield(setup, 'output') && isfield(setup.output, 'folder') && ...
                    isfield(gui.UIControls, 'output') && ...
                    isfield(gui.UIControls.output, 'useCustomFolder') && ...
                    isvalid(gui.UIControls.output.useCustomFolder) && ...
                    ~gui.UIControls.output.useCustomFolder.Value
                setup.output.folder = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
            end
            if isfield(gui.UIControls, 'gasProperties')
                gasFileFields = {'mass', 'harmonicFrequency', 'anharmonicFrequency', ...
                    'rotationalConstant', 'electricQuadrupoleMoment', 'OPBParameter'};
                for idx = 1:length(gasFileFields)
                    fieldName = gasFileFields{idx};
                    if isfield(gui.UIControls.gasProperties, fieldName) && ...
                            isvalid(gui.UIControls.gasProperties.(fieldName))
                        setup.electronKinetics.gasProperties.(fieldName) = gui.UIControls.gasProperties.(fieldName).Items;
                    end
                end
            end
        end
        
        function jsonStruct = convertSetupToJSON(gui, setup)
            % Convert Setup struct to JSON-compatible structure
            % In JSON format, all fields that exist in Setup are included
            % Only remove UI-only fields like isOn within smartGrid and odeSetParameters
            % IMPORTANT: Only include optional fields if their checkboxes are checked (blocked fields should not appear)
            
            jsonStruct = struct();
            fields = fieldnames(setup);
            
            for i = 1:length(fields)
                fieldName = fields{i};
                value = setup.(fieldName);
                
                % Skip UI-only isOn fields in nested structures (handled separately)
                if strcmp(fieldName, 'isOn') && isstruct(value)
                    continue;
                end
                
                if isstruct(value)
                    % Recursively convert nested structures
                    nestedStruct = gui.convertSetupToJSON(value);
                    % Only add if not empty
                    if ~isempty(fieldnames(nestedStruct))
                        jsonStruct.(fieldName) = nestedStruct;
                    end
                elseif iscell(value)
                    % Cell array - keep as cell array (will be JSON array)
                    % Include all arrays (even empty ones, as per JSON format)
                    jsonStruct.(fieldName) = value;
                elseif isnumeric(value) || islogical(value)
                    % Numeric or boolean - keep as is
                    jsonStruct.(fieldName) = value;
                else
                    % String - convert to char if needed
                    jsonStruct.(fieldName) = char(value);
                end
            end
            
            % Special handling: Remove UI-only isOn fields from smartGrid and odeSetParameters
            % But keep the structures themselves if they have other data
            if isfield(jsonStruct, 'electronKinetics') && isfield(jsonStruct.electronKinetics, 'numerics')
                if isfield(jsonStruct.electronKinetics.numerics, 'energyGrid')
                    if isfield(jsonStruct.electronKinetics.numerics.energyGrid, 'smartGrid')
                        % Remove isOn field from smartGrid if it exists (UI-only)
                        if isfield(jsonStruct.electronKinetics.numerics.energyGrid.smartGrid, 'isOn')
                            jsonStruct.electronKinetics.numerics.energyGrid.smartGrid = rmfield(jsonStruct.electronKinetics.numerics.energyGrid.smartGrid, 'isOn');
                        end
                    end
                end
                if isfield(jsonStruct.electronKinetics.numerics, 'nonLinearRoutines')
                    if isfield(jsonStruct.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters')
                        % Remove isOn field from odeSetParameters if it exists (UI-only)
                        if isfield(jsonStruct.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn')
                            jsonStruct.electronKinetics.numerics.nonLinearRoutines.odeSetParameters = rmfield(jsonStruct.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn');
                        end
                    end
                end
            end
            
            % Remove optional fields that are blocked (checkboxes unchecked)
            % Working Conditions optional fields
            if isfield(jsonStruct, 'workingConditions')
                % dischargeCurrent - only if checkbox is checked
                if isfield(gui.UIControls.workingConditions, 'dischargeCurrentCheckbox')
                    if ~gui.UIControls.workingConditions.dischargeCurrentCheckbox.Value
                        if isfield(jsonStruct.workingConditions, 'dischargeCurrent')
                            jsonStruct.workingConditions = rmfield(jsonStruct.workingConditions, 'dischargeCurrent');
                        end
                    end
                end             
                
                
                % gasPressure - only if checkbox is checked
                if isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                    if ~gui.UIControls.workingConditions.gasPressureCheckbox.Value
                        if isfield(jsonStruct.workingConditions, 'gasPressure')
                            jsonStruct.workingConditions = rmfield(jsonStruct.workingConditions, 'gasPressure');
                        end
                    end
                end
                                                
                % electronDensity - only if checkbox is checked
                if isfield(gui.UIControls.workingConditions, 'electronDensityCheckbox')
                    if ~gui.UIControls.workingConditions.electronDensityCheckbox.Value
                        if isfield(jsonStruct.workingConditions, 'electronDensity')
                            jsonStruct.workingConditions = rmfield(jsonStruct.workingConditions, 'electronDensity');
                        end
                    end
                end
            end
            
            % Electron Kinetics optional fields
            if isfield(jsonStruct, 'electronKinetics')
                % shapeParameter - only if eedfType is 'prescribedEedf'
                if isfield(gui.UIControls.electronKinetics, 'eedfType')
                    if ~strcmp(gui.UIControls.electronKinetics.eedfType.Value, 'prescribedEedf')
                        if isfield(jsonStruct.electronKinetics, 'shapeParameter')
                            jsonStruct.electronKinetics = rmfield(jsonStruct.electronKinetics, 'shapeParameter');
                        end
                    end
                end
                
                % LXCatExtraFiles - only if checkbox is checked
                if isfield(gui.UIControls.electronKinetics, 'LXCatExtraCheckbox')
                    if ~gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value
                        if isfield(jsonStruct.electronKinetics, 'LXCatExtraFiles')
                            jsonStruct.electronKinetics = rmfield(jsonStruct.electronKinetics, 'LXCatExtraFiles');
                        end
                    end
                end
                
                % effectiveCrossSectionPopulations - only if checkbox is checked
                if isfield(gui.UIControls.electronKinetics, 'effectivePopCheckbox')
                    if ~gui.UIControls.electronKinetics.effectivePopCheckbox.Value
                        if isfield(jsonStruct.electronKinetics, 'effectiveCrossSectionPopulations')
                            jsonStruct.electronKinetics = rmfield(jsonStruct.electronKinetics, 'effectiveCrossSectionPopulations');
                        end
                    end
                end
                
                % CARgases - only if checkbox is checked
                if isfield(gui.UIControls.electronKinetics, 'CARcheckbox')
                    if ~gui.UIControls.electronKinetics.CARcheckbox.Value
                        if isfield(jsonStruct.electronKinetics, 'CARgases')
                            jsonStruct.electronKinetics = rmfield(jsonStruct.electronKinetics, 'CARgases');
                        end
                    end
                end
                
                % odeSetParameters - only if algorithm is 'temporalIntegration' AND checkbox is checked
                if isfield(jsonStruct.electronKinetics, 'numerics') && isfield(jsonStruct.electronKinetics.numerics, 'nonLinearRoutines')
                    if isfield(jsonStruct.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters')
                        % Check if algorithm is temporalIntegration
                        algorithmIsTemporal = false;
                        if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines, 'algorithm')
                            algorithmIsTemporal = strcmp(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.algorithm.Value, 'temporalIntegration');
                        end
                        
                        % Check if checkbox is checked
                        checkboxChecked = false;
                        if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn')
                            checkboxChecked = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Value;
                        end
                        
                        % Remove if not activated (algorithm not temporalIntegration OR checkbox not checked)
                        if ~algorithmIsTemporal || ~checkboxChecked
                            jsonStruct.electronKinetics.numerics.nonLinearRoutines = rmfield(jsonStruct.electronKinetics.numerics.nonLinearRoutines, 'odeSetParameters');
                        end
                    end
                end
            end
        end

        function disableUnusedOptionalFields(gui, foundFields)
            % Disable optional fields that were not found in the loaded file
            % This ensures that fields not present in the file appear blocked/unselected
            % Also activate fields that WERE found in the file
            
            % Helper function to check if a field was found
            function wasFound = isFieldFound(fieldPath)
                wasFound = foundFields.isKey(fieldPath) && foundFields(fieldPath);
            end
            
            % Working Conditions optional fields - ACTIVATE if found
            
            
            if isFieldFound('workingConditions.gasPressure')
                try
                    gui.toggleGasPressureEnable(true);
                    if isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                        gui.UIControls.workingConditions.gasPressureCheckbox.Value = true;
                    end
                catch
                end
            elseif ~isFieldFound('workingConditions.gasPressure')
                try
                    gui.toggleGasPressureEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                        gui.UIControls.workingConditions.gasPressureCheckbox.Value = false;
                    end
                catch
                end
            end
            
                        
            if isFieldFound('workingConditions.electronDensity')
                try
                    gui.toggleElectronDensityEnable(true);
                    if isfield(gui.UIControls.workingConditions, 'electronDensityCheckbox')
                        gui.UIControls.workingConditions.electronDensityCheckbox.Value = true;
                    end
                catch
                end
            elseif ~isFieldFound('workingConditions.electronDensity')
                try
                    gui.toggleElectronDensityEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'electronDensityCheckbox')
                        gui.UIControls.workingConditions.electronDensityCheckbox.Value = false;
                    end
                catch
                end
            end
            
            % Working Conditions optional fields - DISABLE if not found
            if ~isFieldFound('workingConditions.dischargeCurrent')
                try
                    gui.toggleDischargeCurrentEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'dischargeCurrentCheckbox')
                        gui.UIControls.workingConditions.dischargeCurrentCheckbox.Value = false;
                    end
                catch
                end
            end
            
            if ~isFieldFound('workingConditions.dischargePowerDensity')
                try
                    gui.toggleDischargePowerEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'dischargePowerDensityCheckbox')
                        gui.UIControls.workingConditions.dischargePowerDensityCheckbox.Value = false;
                    end
                catch
                end
            end
                                    
            
            if ~isFieldFound('workingConditions.gasPressure')
                try
                    gui.toggleGasPressureEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'gasPressureCheckbox')
                        gui.UIControls.workingConditions.gasPressureCheckbox.Value = false;
                    end
                catch
                end
            end
            
            if ~isFieldFound('workingConditions.electronDensity')
                try
                    gui.toggleElectronDensityEnable(false);
                    if isfield(gui.UIControls.workingConditions, 'electronDensityCheckbox')
                        gui.UIControls.workingConditions.electronDensityCheckbox.Value = false;
                    end
                catch
                end
            end
            
            % Electron Kinetics optional fields
            % Check both possible field names (LXCatExtraFiles and LXCatFilesExtra)
            if ~isFieldFound('electronKinetics.LXCatExtraFiles') && ~isFieldFound('electronKinetics.LXCatFilesExtra')
                try
                    gui.toggleLXCatExtraEnable(false);
                    if isfield(gui.UIControls.electronKinetics, 'LXCatExtraCheckbox')
                        gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value = false;
                    end
                catch
                end
            end
            
            if ~isFieldFound('electronKinetics.effectiveCrossSectionPopulations')
                try
                    gui.toggleEffectivePopEnable(false);
                    if isfield(gui.UIControls.electronKinetics, 'effectivePopCheckbox')
                        gui.UIControls.electronKinetics.effectivePopCheckbox.Value = false;
                    end
                catch
                end
            end
            
            if ~isFieldFound('electronKinetics.CARgases')
                try
                    gui.toggleCARGasEnable(false);
                    if isfield(gui.UIControls.electronKinetics, 'CARcheckbox')
                        gui.UIControls.electronKinetics.CARcheckbox.Value = false;
                    end
                catch
                end
            end
            
            % Output optional field
            if ~isFieldFound('output.isOn')
                try
                    gui.toggleOutputEnable(false);
                    if isfield(gui.UIControls.output, 'isOn')
                        gui.UIControls.output.isOn.Value = false;
                    end
                catch
                end
            end
            
            % Clear lists that were not found in the file
            % LXCatFiles
            if ~isFieldFound('electronKinetics.LXCatFiles')
                gui.Setup.electronKinetics.LXCatFiles = {};
                try
                    if isfield(gui.UIControls.electronKinetics, 'LXCatFiles')
                        gui.UIControls.electronKinetics.LXCatFiles.Items = {};
                    end
                catch
                end
            end
            
            % Gas Properties lists
            if ~isFieldFound('electronKinetics.gasProperties.fraction')
                gui.Setup.electronKinetics.gasProperties.fraction = {};
                try
                    if isfield(gui.UIControls, 'gasProperties') && isfield(gui.UIControls.gasProperties, 'fraction')
                        gui.UIControls.gasProperties.fraction.Items = {};
                    end
                catch
                end
            end
            
            % State Properties lists
            if ~isFieldFound('electronKinetics.stateProperties.energy')
                gui.Setup.electronKinetics.stateProperties.energy = {};
                try
                    if isfield(gui.UIControls.electronKinetics, 'stateProperties') && isfield(gui.UIControls.electronKinetics.stateProperties, 'energy')
                        gui.UIControls.electronKinetics.stateProperties.energy.Items = {};
                    end
                catch
                end
            end
            
            if ~isFieldFound('electronKinetics.stateProperties.statisticalWeight')
                gui.Setup.electronKinetics.stateProperties.statisticalWeight = {};
                try
                    if isfield(gui.UIControls.electronKinetics, 'stateProperties') && isfield(gui.UIControls.electronKinetics.stateProperties, 'statisticalWeight')
                        gui.UIControls.electronKinetics.stateProperties.statisticalWeight.Items = {};
                    end
                catch
                end
            end
            
            if ~isFieldFound('electronKinetics.stateProperties.population')
                gui.Setup.electronKinetics.stateProperties.population = {};
                try
                    if isfield(gui.UIControls.electronKinetics, 'stateProperties') && isfield(gui.UIControls.electronKinetics.stateProperties, 'population')
                        gui.UIControls.electronKinetics.stateProperties.population.Items = {};
                    end
                catch
                end
            end
            
            % Output dataSets
            if ~isFieldFound('output.dataSets')
                gui.Setup.output.dataSets = {};
            end
        end

        function activateFoundOptionalFields(gui, foundFields)
            % Activate and expand optional fields that were found in the loaded file
            
            % Helper function to check if a field was found (or any of its subfields)
            function wasFound = isFieldFound(fieldPath)
                wasFound = foundFields.isKey(fieldPath) && foundFields(fieldPath);
                % Also check if any subfield was found (e.g., smartGrid.minEedfDecay means smartGrid was found)
                if ~wasFound
                    keys = foundFields.keys();
                    for k = 1:length(keys)
                        if startsWith(keys{k}, [fieldPath, '.'])
                            wasFound = true;
                            break;
                        end
                    end
                end
            end
            
            % Smart Grid - if found (or any subfield), activate checkbox and enable controls
            if isFieldFound('electronKinetics.numerics.energyGrid.smartGrid') || ...
               isFieldFound('electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay') || ...
               isFieldFound('electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay') || ...
               isFieldFound('electronKinetics.numerics.energyGrid.smartGrid.updateFactor')
                try
                    if isfield(gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid, 'isOn')
                        gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.isOn.Value = true;
                        gui.toggleSmartGridEnable(true);
                    end
                catch ME
                    warning('Error activating smartGrid: %s', ME.message);
                end
            end
            
            % ODE Set Parameters - if found (or any subfield), show button and expand panel
            % Button is only visible when algorithm is temporalIntegration
            try
                advancedButton = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton;
                advancedPanel = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel;
                
                % Check current algorithm
                algorithm = gui.Setup.electronKinetics.numerics.nonLinearRoutines.algorithm;
                
                if strcmp(algorithm, 'temporalIntegration')
                    % Show button when algorithm is temporalIntegration
                    advancedButton.Visible = 'on';
                    
                    if isFieldFound('electronKinetics.numerics.nonLinearRoutines.odeSetParameters') || ...
                       isFieldFound('electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol') || ...
                       isFieldFound('electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol') || ...
                       isFieldFound('electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep')
                        % Expand panel if odeSetParameters is found
                        advancedPanel.Visible = 'on';
                        advancedButton.Text = 'Advanced (Optional) ▲';
                        
                        % Activate the checkbox
                        if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn')
                            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Value = true;
                            gui.toggleOdeParametersEnable(true);
                        end
                    else
                        % Collapse panel if odeSetParameters is not found
                        advancedPanel.Visible = 'off';
                        advancedButton.Text = 'Advanced (Optional) ▼';
                    end
                else
                    % Hide button when algorithm is not temporalIntegration
                    advancedButton.Visible = 'off';
                    advancedPanel.Visible = 'off';
                    advancedButton.Text = 'Advanced (Optional) ▼';
                end
            catch ME
                warning('Error configuring advanced options: %s', ME.message);
            end
            
            % Effective Cross Section Populations - if found, activate checkbox
            if isFieldFound('electronKinetics.effectiveCrossSectionPopulations')
                try
                    gui.toggleEffectivePopEnable(true);
                    if isfield(gui.UIControls.electronKinetics, 'effectivePopCheckbox')
                        gui.UIControls.electronKinetics.effectivePopCheckbox.Value = true;
                    end
                catch ME
                    warning('Error activating effectiveCrossSectionPopulations: %s', ME.message);
                end
            end
        end

        function saveInputFile(gui, ~, ~)
            % Save the current setup to a file (.in or .json)
            defaultName = ['LoKI_Setup_', datestr(now,'yyyymmdd_HHMMSS'), '.in'];
            % Point to Code/Input folder
            inputFolder = fullfile(pwd, 'Input');
            if ~isfolder(inputFolder)
                inputFolder = pwd; % Fallback to current directory
            end
            % Allow both .in and .json files
            [file, path] = uiputfile({'*.in;*.json', 'LoKI Input Files (*.in, *.json)'; '*.in', 'LoKI Input Files (*.in)'; '*.json', 'LoKI JSON Files (*.json)'}, 'Save LoKI Input File As', fullfile(inputFolder, defaultName));

            if isequal(file, 0) || isequal(path, 0)
                uialert(gui.Fig, 'File save cancelled.', 'Cancelled');
                return; % User cancelled
            end

            outputFile = fullfile(path, file);
            [~, ~, ext] = fileparts(file);
            
            try
                % Generate based on file extension
                if strcmpi(ext, '.json')
                    gui.generateJSONFile(outputFile);
                else
                    gui.generateInputFile(outputFile);
                end
                uialert(gui.Fig, ['Input file saved successfully: ', outputFile], 'Save Successful', 'Icon', 'success');
            catch ME
                uialert(gui.Fig, ['Error saving file: ', ME.message], 'Save Error');
            end
        end


        function generateInputFile(gui, filename)
            % Generates the text input file from the gui.Setup struct
            fid = fopen(filename, 'w');
            if fid == -1
                error('Cannot open file "%s" for writing.', filename);
            end
            fprintf(fid, '%% LoKI Input File generated by InputGUILokib on %s %%\n\n', datestr(now));

            % CARgases will be populated from the listbox items automatically when writing
            setup = gui.getSetupForExport();

            % Use recursion or explicit handling for nested structs
            sections = fieldnames(setup);
            for i = 1:length(sections)
                sectionName = sections{i};
                if strcmp(sectionName, 'electronKinetics')
                    % Special handling for electronKinetics to reorder fields
                    fprintf(fid, '%s:\n', sectionName); % Section header
                    gui.writeElectronKineticsContent(setup.(sectionName), '  ', fid);
                else
                    fprintf(fid, '%s:\n', sectionName); % Section header
                    gui.writeStructContent(setup.(sectionName), '  ', fid); % Indent content
                end
                fprintf(fid, '\n'); % Blank line between sections
            end

            fclose(fid);
            fprintf('Input file generated: %s\n', filename);
        end

        function writeElectronKineticsContent(gui, ek, indent, fid)
            % Write electronKinetics section with proper field order
            % First write basic fields
            if isfield(ek, 'isOn')
                fprintf(fid, '%sisOn: %s\n', indent, mat2str(ek.isOn));
            end
            if isfield(ek, 'eedfType')
                fprintf(fid, '%seedfType: %s\n', indent, ek.eedfType);
            end
            if isfield(ek, 'shapeParameter')
                if strcmp(gui.UIControls.electronKinetics.eedfType.Value, 'prescribedEedf')
                    fprintf(fid, '%sshapeParameter: %f\n', indent, ek.shapeParameter);
                end
                % Skip if not prescribedEedf
       
            end
            if isfield(ek, 'ionizationOperatorType')
                fprintf(fid, '%sionizationOperatorType: %s\n', indent, ek.ionizationOperatorType);
            end
            if isfield(ek, 'growthModelType')
                fprintf(fid, '%sgrowthModelType: %s\n', indent, ek.growthModelType);
            end
            if isfield(ek, 'includeEECollisions')
                fprintf(fid, '%sincludeEECollisions: %s\n', indent, mat2str(ek.includeEECollisions));
            end
            
            % LXCat files (only write if not empty)
            if isfield(ek, 'LXCatFiles') && ~isempty(ek.LXCatFiles)
                fprintf(fid, '%sLXCatFiles:\n', indent);
                for j = 1:length(ek.LXCatFiles)
                    fprintf(fid, '%s  - %s\n', indent, ek.LXCatFiles{j});
                end
            end
            
            % LXCat Extra Files (only if checked)
            if isfield(ek, 'LXCatExtraFiles')
                try
                    if gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value
                        fprintf(fid, '%sLXCatExtraFiles:\n', indent);
                        extraFiles = gui.UIControls.electronKinetics.LXCatExtraFiles.Items;
                        if isempty(extraFiles)
                            fprintf(fid, '%s   -\n', indent);
                        else
                            for j = 1:length(extraFiles)
                                fprintf(fid, '%s   - %s\n', indent, extraFiles{j});
                            end
                        end
                    end
                catch
                    % Skip
                end
            end
            
            % Effective Cross Section Populations (only if checked)
            if isfield(ek, 'effectiveCrossSectionPopulations')
                try
                    if gui.UIControls.electronKinetics.effectivePopCheckbox.Value
                        fprintf(fid, '%seffectiveCrossSectionPopulations:\n', indent);
                        effectivePops = gui.UIControls.electronKinetics.effectiveCrossSectionPopulations.Items;
                        if isempty(effectivePops)
                            fprintf(fid, '%s   -\n', indent);
                        else
                            for j = 1:length(effectivePops)
                                fprintf(fid, '%s   - %s\n', indent, effectivePops{j});
                            end
                        end
                    end
                catch
                    % Skip
                end
            end
            
            % CAR gases (only if checked)
            if isfield(ek, 'CARgases')
                try
                    if gui.UIControls.electronKinetics.CARcheckbox.Value
                        fprintf(fid, '%sCARgases:\n', indent);
                        carGases = gui.UIControls.electronKinetics.CARgases.Items;
                        if isempty(carGases)
                            fprintf(fid, '%s   -\n', indent);
                        else
                            for j = 1:length(carGases)
                                fprintf(fid, '%s   - %s\n', indent, carGases{j});
                            end
                        end
                    end
                catch
                    % Skip
                end
            end
            
            % Gas Properties (write before stateProperties, only if fraction is not empty)
            if isfield(ek, 'gasProperties')
                % Check if fraction list is not empty
                hasFraction = false;
                if isfield(gui.UIControls, 'gasProperties') && isfield(gui.UIControls.gasProperties, 'fraction')
                    fractionItems = gui.UIControls.gasProperties.fraction.Items;
                    hasFraction = ~isempty(fractionItems);
                elseif isfield(ek.gasProperties, 'fraction')
                    hasFraction = ~isempty(ek.gasProperties.fraction);
                end
                
                if hasFraction
                    fprintf(fid, '%sgasProperties:\n', indent);
                    % Write gas properties manually to handle fraction from UI
                    gasFileFields = {'mass', 'harmonicFrequency', 'anharmonicFrequency', ...
                        'rotationalConstant', 'electricQuadrupoleMoment', 'OPBParameter'};
                    for fileIdx = 1:length(gasFileFields)
                        gasFileField = gasFileFields{fileIdx};
                        if isfield(ek.gasProperties, gasFileField)
                            gui.writeNamedList(fid, indent, gasFileField, ek.gasProperties.(gasFileField));
                        end
                    end
                    fprintf(fid, '%s  fraction:\n', indent);
                    if isfield(gui.UIControls, 'gasProperties') && isfield(gui.UIControls.gasProperties, 'fraction')
                        fractionItems = gui.UIControls.gasProperties.fraction.Items;
                        for j = 1:length(fractionItems)
                            fprintf(fid, '%s    - %s\n', indent, fractionItems{j});
                        end
                    elseif isfield(ek.gasProperties, 'fraction')
                        for j = 1:length(ek.gasProperties.fraction)
                            fprintf(fid, '%s    - %s\n', indent, ek.gasProperties.fraction{j});
                        end
                    end
                end
            end
            
            % State Properties (write normally with writeStructContent, but check if any lists are non-empty)
            if isfield(ek, 'stateProperties')
                % Check if any state property list is non-empty
                hasStateProps = false;
                if isfield(ek.stateProperties, 'energy') && ~isempty(ek.stateProperties.energy)
                    hasStateProps = true;
                elseif isfield(ek.stateProperties, 'statisticalWeight') && ~isempty(ek.stateProperties.statisticalWeight)
                    hasStateProps = true;
                elseif isfield(ek.stateProperties, 'population') && ~isempty(ek.stateProperties.population)
                    hasStateProps = true;
                end
                
                if hasStateProps
                    fprintf(fid, '%sstateProperties:\n', indent);
                    gui.writeStructContent(ek.stateProperties, [indent, '  '], fid);
                end
            end
            
            % Numerics (write normally with writeStructContent)
            if isfield(ek, 'numerics')
                fprintf(fid, '%snumerics:\n', indent);
                gui.writeStructContent(ek.numerics, [indent, '  '], fid);
            end
        end

        function items = asCellArray(~, value)
            if isempty(value)
                items = {};
            elseif iscell(value)
                items = value;
            elseif isstring(value)
                items = cellstr(value(:));
            elseif ischar(value)
                items = {value};
            else
                try
                    items = {char(value)};
                catch
                    items = {num2str(value)};
                end
            end
        end

        function writeNamedList(gui, fid, indent, fieldName, value)
            items = gui.asCellArray(value);
            if isempty(items)
                return;
            end
            fprintf(fid, '%s  %s:\n', indent, fieldName);
            for itemIdx = 1:length(items)
                fprintf(fid, '%s    - %s\n', indent, gui.formatValue(items{itemIdx}));
            end
        end

        function writeStructContent(gui, s, indent, fid)
            % Helper function to recursively write struct content to file
            fields = fieldnames(s);
            for i = 1:length(fields)
                fieldName = fields{i};
                value = s.(fieldName);

                if isstruct(value)
                    % Special handling for smartGrid and odeSetParameters
                    if strcmp(fieldName, 'smartGrid')
                        % Check if smartGrid should be written by checking the checkbox state
                        if gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.isOn.Value
                            % If checkbox is checked, write the section with current values
                            fprintf(fid, '%s%s:\n', indent, fieldName);
                            % Write the actual values from the GUI controls
                            minDecay = gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay.Value;
                            maxDecay = gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay.Value;
                            updateFac = gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Value;
                            fprintf(fid, '%s  minEedfDecay: %s\n', indent, num2str(minDecay));
                            fprintf(fid, '%s  maxEedfDecay: %s\n', indent, num2str(maxDecay));
                            fprintf(fid, '%s  updateFactor: %s\n', indent, num2str(updateFac));
                        end
                        % If checkbox not checked, don't write anything
                    elseif strcmp(fieldName, 'odeSetParameters')
                        % Check if ODE parameters should be written:
                        % 1. Algorithm must be 'temporalIntegration'
                        % 2. Checkbox must be checked
                        algorithmIsTemporal = false;
                        if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines, 'algorithm')
                            algorithmIsTemporal = strcmp(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.algorithm.Value, 'temporalIntegration');
                        end
                        
                        checkboxChecked = false;
                        if isfield(gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters, 'isOn')
                            checkboxChecked = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Value;
                        end
                        
                        if algorithmIsTemporal && checkboxChecked
                            % If algorithm is temporalIntegration AND checkbox is checked, write the section
                            fprintf(fid, '%s%s:\n', indent, fieldName);
                            % Write the actual values from the GUI controls with extra indentation
                            absTol = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Value;
                            relTol = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Value;
                            maxStep = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep.Value;
                            fprintf(fid, '%s   AbsTol: %s\n', indent, num2str(absTol));
                            fprintf(fid, '%s   RelTol: %s\n', indent, num2str(relTol));
                            fprintf(fid, '%s   MaxStep: %s\n', indent, num2str(maxStep));
                        end
                        % If not enabled (algorithm not temporalIntegration OR checkbox not checked), don't write anything
                    else
                        fprintf(fid, '%s%s:\n', indent, fieldName);
                        gui.writeStructContent(value, [indent, '  '], fid); % Recurse with more indent
                    end
                elseif iscell(value) % Handle lists (cell arrays)
                    % Skip empty optional fields (checkboxes control these)
                    if isempty(value)
                        if strcmp(fieldName, 'LXCatExtraFiles')
                            % Check if checkbox is checked
                            try
                                if gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value
                                    fprintf(fid, '%s%s:\n', indent, fieldName);
                                    fprintf(fid, '%s  -\n', indent);
                                end
                            catch
                                % Skip if not checked
                            end
                        elseif strcmp(fieldName, 'effectiveCrossSectionPopulations')
                            try
                                if gui.UIControls.electronKinetics.effectivePopCheckbox.Value
                                    fprintf(fid, '%s%s:\n', indent, fieldName);
                                    fprintf(fid, '%s  -\n', indent);
                                end
                            catch
                                % Skip if not checked
                            end
                        elseif strcmp(fieldName, 'CARgases')
                            try
                                if gui.UIControls.electronKinetics.CARcheckbox.Value
                                    fprintf(fid, '%s%s:\n', indent, fieldName);
                                    fprintf(fid, '%s  -\n', indent);
                                end
                            catch
                                % Skip if not checked
                            end
                        end
                        % If none of the above, don't write empty fields (e.g., energy, statisticalWeight, population if empty)
                    else
                        % Only write non-empty lists
                        fprintf(fid, '%s%s:\n', indent, fieldName);
                        for j = 1:length(value)
                            fprintf(fid, '%s  - %s\n', indent, gui.formatValue(value{j}));
                        end
                    end
                elseif strcmp(fieldName, 'CARgas')
                    % Only write CARgas if checkbox is checked and value is not empty
                    try
                        if gui.UIControls.electronKinetics.CARcheckbox.Value && ~isempty(value)
                            fprintf(fid, '%sCARgases:\n', indent);
                            fprintf(fid, '%s  - %s\n', indent, value);
                        end
                    catch
                        % Skip if CARcheckbox doesn't exist
                    end
                elseif islogical(value)
                    if value
                        fprintf(fid, '%s%s: true\n', indent, fieldName);
                    else
                        fprintf(fid, '%s%s: false\n', indent, fieldName);
                    end
                elseif isnumeric(value)
                    % Handle optional discharge fields (must check checkbox to include)
                    if strcmp(fieldName, 'dischargeCurrent')
                        try
                            if gui.UIControls.workingConditions.dischargeCurrentCheckbox.Value
                                fprintf(fid, '%s%s: %s\n', indent, fieldName, num2str(value));
                            end
                        catch
                            % Skip if checkbox doesn't exist
                        end
                    elseif isempty(value)
                        % For other empty numeric fields, skip writing
                    elseif isscalar(value)
                        fprintf(fid, '%s%s: %s\n', indent, fieldName, num2str(value));
                    else
                        fprintf(fid, '%s%s: %s\n', indent, fieldName, mat2str(value));
                    end
                elseif ischar(value) || isstring(value)
                    fprintf(fid, '%s%s: %s\n', indent, fieldName, value);
                else
                    warning('Unhandled value type for field %s: %s', fieldName, class(value));
                end
            end
        end

        function strValue = formatValue(~, value)
            % Format different value types for the output file
            if ischar(value) || isstring(value)
                strValue = char(value); % Ensure it's char, handle strings/expressions
            elseif isnumeric(value)
                if isscalar(value)
                    strValue = num2str(value, '%.8g'); % Format numbers nicely
                else
                    strValue = mat2str(value); % Use mat2str for arrays/vectors
                end
            elseif islogical(value) % Should be handled before calling formatValue
                if value
                    strValue = 'true';
                else
                    strValue = 'false';
                end
            else
                strValue = ''; % Fallback for unknown types
            end
        end

        function updateDataSetsSelection(gui, checkbox, dataSetName, isSelected)
            % Update the data sets selection in the Setup struct
            currentDataSets = gui.Setup.output.dataSets;
            
            if isSelected
                % Add to selection if not already there
                if ~ismember(dataSetName, currentDataSets)
                    gui.Setup.output.dataSets = [currentDataSets, {dataSetName}];
                end
            else
                % Remove from selection
                gui.Setup.output.dataSets = currentDataSets(~strcmp(currentDataSets, dataSetName));
            end
        end

        function openSetupFileHelp(~, ~, ~)
            web("../Documentation/html/The_setup_file.html");
        end

        function runSimulation(gui, ~, ~)
            % Show loading indicator while launching LoKI Simulation Tool
            d = uiprogressdlg(gui.Fig, 'Title', 'LoKI-B', ...
                'Message', 'Launching LoKI Simulation Tool...', ...
                'Indeterminate', 'on');
            cleanupDlg = onCleanup(@() gui.safeCloseDlg(d));
            if isfield(gui.UIControls, 'runButton') && ~isempty(gui.UIControls.runButton) && isvalid(gui.UIControls.runButton)
                gui.UIControls.runButton.Enable = 'off';
                cleanupBtn = onCleanup(@() gui.safeEnableBtn(gui.UIControls.runButton));
            end

            % 1. Generate the input file (e.g., to a temporary location or specific output)
            
            tempInputFile = fullfile(['Input' filesep 'loki_run_' datestr(now,'yyyymmdd_HHMMSSFFF') '.txt']);
            try
                gui.generateInputFile(tempInputFile);
            catch ME
                uialert(gui.Fig, ['Error generating temporary input file: ', ME.message], 'Run Error');
                return;
            end

            % 2. thanks gpt! Call the actual LoKI Boltzmann solver
            %    This is where you link to your solver code.
            %    You might need to pass the path to the generated input file.
            fprintf('--- Running Simulation (Conceptual) ---\n');
            fprintf('Using input file: %s\n', tempInputFile);
            fprintf('Setup details:\n');
            disp(gui.Setup); % Display current setup in console for debugging

            % try
                % --- Replace this with the actual call to your solver ---
                % Example: status = run_loki_solver(tempInputFile);
                % Launch simulation tool
                lokibcl(tempInputFile(7:end));
                % Run the simulation
                
            % catch ME
            %     uialert(gui.Fig, ['Error occurred during simulation: ', ME.message], 'Simulation Runtime Error');
            % end

            % Optional: Clean up temporary file
            % delete(tempInputFile);
            fprintf('--- Simulation Finished ---\n');
         end

        function safeCloseDlg(~, d)
            try
                if ~isempty(d) && isvalid(d)
                    close(d);
                end
            catch
            end
        end

        function safeEnableBtn(~, btn)
            try
                if ~isempty(btn) && isvalid(btn)
                    btn.Enable = 'on';
                end
            catch
            end
        end
    end % methods
end % classdef
