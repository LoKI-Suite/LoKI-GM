% --- LoKI_GUI.m ---
classdef LoKI_GUI < handle
    %LoKI_GUI Class that provides a user-friendly interface for LoKI setup

    properties
        Setup;      % Struct holding the simulation setup data
        Fig;        % Main figure handle
        UIControls; % Struct to hold handles to UI controls for easy access/update
    end

    methods
        function gui = LoKI_GUI(inputFile)
            % Constructor
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
            gui.createGUI();
            gui.populateGUIFromSetup(); % Populate GUI fields with Setup data
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
            gui.Setup.electronKinetics.shapeParameter = 'Maxwellian'; % Only used for prescribedEedf
            gui.Setup.electronKinetics.ionizationOperatorType = 'usingSDCS'; % 'conservative', 'oneTakesAll', 'usingSDCS'
            gui.Setup.electronKinetics.growthModelType = 'temporal'; % 'temporal' or 'spatial'
            gui.Setup.electronKinetics.includeEECollisions = false;
            gui.Setup.electronKinetics.LXCatFiles = {'Nitrogen/N2_LXCat.txt', 'Nitrogen/N2_rot_LXCat.txt'}; % Cell array of strings
            gui.Setup.electronKinetics.LXCatExtraFiles = {'Nitrogen/extra_LXCat.txt'}; % Optional
            gui.Setup.electronKinetics.effectiveCrossSectionPopulations = {'Nitrogen/N2_effectivePop.txt'}; % Optional
            gui.Setup.electronKinetics.CARgases = {}; % Optional
            
            % Gas Properties [cite: 5] - Inside electronKinetics
            gui.Setup.electronKinetics.gasProperties.mass = 'Databases/masses.txt';
            gui.Setup.electronKinetics.gasProperties.fraction = {'N2 = 1'}; % Cell array of strings
            gui.Setup.electronKinetics.gasProperties.harmonicFrequency = 'Databases/harmonicFrequencies.txt';
            gui.Setup.electronKinetics.gasProperties.anharmonicFrequency = 'Databases/anharmonicFrequencies.txt';
            gui.Setup.electronKinetics.gasProperties.rotationalConstant = 'Databases/rotationalConstants.txt';
            gui.Setup.electronKinetics.gasProperties.electricQuadrupoleMoment = 'Databases/quadrupoleMoment.txt';
            gui.Setup.electronKinetics.gasProperties.OPBParameter = 'Databases/OPBParameter.txt';
            
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
            gui.Setup.output.dataSets = {'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTable'}; % Cell array
        end

        function createGUI(gui)
            % Create the main figure
            iconPath = 'icon.png'; % Relative or absolute path to your icon
            if ~isfile(iconPath)
                warning('Icon file not found: %s. Using default icon.', iconPath);
                iconPath = ''; % Use default if not found
            end

            gui.Fig = uifigure('Name', 'LoKI-B', ...
                'Position', [100, 100, 850, 650], ...
                'NumberTitle', 'off', ...
                'Resize', 'on', ... % Allow resizing
                'Icon', iconPath, ...
                'WindowStyle', 'normal'); % Keep window visible

            % Main grid layout
            mainGrid = uigridlayout(gui.Fig, [2, 1]);
            mainGrid.RowHeight = {'1x', 'fit'}; % Tabs take most space, buttons at bottom

            % Create tabs
            tabGroup = uitabgroup(mainGrid);
            tabGroup.Layout.Row = 1;
            tabGroup.Layout.Column = 1;

            % --- Working Conditions Tab ---
            workingTab = uitab(tabGroup, 'Title', 'Working Conditions', 'Scrollable', 'on');
            gui.createWorkingConditionsPanel(workingTab);

            % --- Electron Kinetics Tab ---
            kineticsTab = uitab(tabGroup, 'Title', 'Electron Kinetics', 'Scrollable', 'on');
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

            % --- Button Panel ---
            buttonPanel = uipanel(mainGrid, 'BorderType', 'none');
            buttonPanel.Layout.Row = 2;
            buttonPanel.Layout.Column = 1;
            buttonGrid = uigridlayout(buttonPanel, [1, 4]);
            % Add some padding/spacing if needed
            buttonGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'}; % Push buttons to right
            buttonGrid.Padding = [10 10 10 10];
            buttonGrid.ColumnSpacing = 10;

            % Load Button (Placeholder)
            uibutton(buttonGrid, 'Text', 'Load Settings (Soon...)', ...
                'ButtonPushedFcn', @gui.loadSettings, 'Enable', 'off'); % Disabled for now

            % Save Button
            uibutton(buttonGrid, 'Text', 'Save Input File', ...
                'ButtonPushedFcn', @gui.saveInputFile);

            % Run Button
            uibutton(buttonGrid, 'Text', 'Generate & Run', ...
                'FontWeight', 'bold', ...
                'ButtonPushedFcn', @gui.runSimulation);
        end

        function createWorkingConditionsPanel(gui, parent)
            % Use grid layout for better alignment and resizing
            grid = uigridlayout(parent, [15, 3]); % Adjust rows as needed
            grid.ColumnWidth = {'fit', '1x', 'fit'}; % Label, Edit, Browse/Unit
            grid.RowHeight = repmat({'fit'}, 1, 15);
            grid.Padding = [10 10 10 10];
            grid.RowSpacing = 5;
            grid.ColumnSpacing = 10;

            row = 1;
            % Reduced Field [cite: 1]
            reducedFieldLabel = uilabel(grid, 'Text', 'Reduced Field (Td):');
            reducedFieldLabel.Layout.Row = row;
            reducedFieldLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.reducedField = uieditfield(grid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.reducedField', evt.Value));
            gui.UIControls.workingConditions.reducedField.Layout.Row = row;
            gui.UIControls.workingConditions.reducedField.Layout.Column = 2;
            reducedFieldHint = uilabel(grid, 'Text', '(e.g., 10 or logspace(1,2,10))'); % Hint
            reducedFieldHint.Layout.Row = row;
            reducedFieldHint.Layout.Column = 3;

            row = row + 1;
            % Electron Temperature [cite: 1]
            electronTempLabel = uilabel(grid, 'Text', 'Electron Temperature (eV):');
            electronTempLabel.Layout.Row = row;
            electronTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.electronTemperature = uieditfield(grid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.electronTemperature', evt.Value));
            gui.UIControls.workingConditions.electronTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.electronTemperature.Layout.Column = 2;
            electronTempHint = uilabel(grid, 'Text', '(e.g., 1.5 or linspace(0.1, 5, 20))'); % Hint
            electronTempHint.Layout.Row = row;
            electronTempHint.Layout.Column = 3;

            row = row + 1;
            % Excitation Frequency [cite: 1]
            excitationFreqLabel = uilabel(grid, 'Text', 'Excitation Frequency (Hz):');
            excitationFreqLabel.Layout.Row = row;
            excitationFreqLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.excitationFrequency = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.excitationFrequency', evt.Value));
            gui.UIControls.workingConditions.excitationFrequency.Layout.Row = row;
            gui.UIControls.workingConditions.excitationFrequency.Layout.Column = 2;
            excitationFreqHint = uilabel(grid, 'Text', ''); % Hint
            excitationFreqHint.Layout.Row = row;
            excitationFreqHint.Layout.Column = 3;

            row = row + 1;
            % Gas Pressure [cite: 1]
            gasPressureLabel = uilabel(grid, 'Text', 'Gas Pressure (Pa):');
            gasPressureLabel.Layout.Row = row;
            gasPressureLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.gasPressure = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.gasPressure', evt.Value));
            gui.UIControls.workingConditions.gasPressure.Layout.Row = row;
            gui.UIControls.workingConditions.gasPressure.Layout.Column = 2;
            gasPressureHint = uilabel(grid, 'Text', ''); % Hint
            gasPressureHint.Layout.Row = row;
            gasPressureHint.Layout.Column = 3;

            row = row + 1;
            % Gas Temperature [cite: 1]
            gasTempLabel = uilabel(grid, 'Text', 'Gas Temperature (K):');
            gasTempLabel.Layout.Row = row;
            gasTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.gasTemperature = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.gasTemperature', evt.Value));
            gui.UIControls.workingConditions.gasTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.gasTemperature.Layout.Column = 2;
            gasTempHint = uilabel(grid, 'Text', ''); % Hint
            gasTempHint.Layout.Row = row;
            gasTempHint.Layout.Column = 3;

            row = row + 1;
            % Wall Temperature [cite: 2]
            wallTempLabel = uilabel(grid, 'Text', 'Wall Temperature (K):');
            wallTempLabel.Layout.Row = row;
            wallTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.wallTemperature = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.wallTemperature', evt.Value));
            gui.UIControls.workingConditions.wallTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.wallTemperature.Layout.Column = 2;
            wallTempHint = uilabel(grid, 'Text', ''); % Hint
            wallTempHint.Layout.Row = row;
            wallTempHint.Layout.Column = 3;

            row = row + 1;
            % External Temperature [cite: 2]
            extTempLabel = uilabel(grid, 'Text', 'External Temperature (K):');
            extTempLabel.Layout.Row = row;
            extTempLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.extTemperature = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.extTemperature', evt.Value));
            gui.UIControls.workingConditions.extTemperature.Layout.Row = row;
            gui.UIControls.workingConditions.extTemperature.Layout.Column = 2;
            extTempHint = uilabel(grid, 'Text', ''); % Hint
            extTempHint.Layout.Row = row;
            extTempHint.Layout.Column = 3;

            row = row + 1;
            % Surface Site Density [cite: 2]
            surfaceSiteDensityLabel = uilabel(grid, 'Text', 'Surface Site Density (m^-2):');
            surfaceSiteDensityLabel.Layout.Row = row;
            surfaceSiteDensityLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.surfaceSiteDensity = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.surfaceSiteDensity', evt.Value));
            gui.UIControls.workingConditions.surfaceSiteDensity.Layout.Row = row;
            gui.UIControls.workingConditions.surfaceSiteDensity.Layout.Column = 2;
            surfaceSiteDensityHint = uilabel(grid, 'Text', ''); % Hint
            surfaceSiteDensityHint.Layout.Row = row;
            surfaceSiteDensityHint.Layout.Column = 3;

            row = row + 1;
            % Electron Density [cite: 3]
            electronDensityLabel = uilabel(grid, 'Text', 'Electron Density (m^-3):');
            electronDensityLabel.Layout.Row = row;
            electronDensityLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.electronDensity = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.electronDensity', evt.Value));
            gui.UIControls.workingConditions.electronDensity.Layout.Row = row;
            gui.UIControls.workingConditions.electronDensity.Layout.Column = 2;
            electronDensityHint = uilabel(grid, 'Text', ''); % Hint
            electronDensityHint.Layout.Row = row;
            electronDensityHint.Layout.Column = 3;

            row = row + 1;
            % Chamber Length [cite: 3]
            chamberLengthLabel = uilabel(grid, 'Text', 'Chamber Length (m):');
            chamberLengthLabel.Layout.Row = row;
            chamberLengthLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.chamberLength = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.chamberLength', evt.Value));
            gui.UIControls.workingConditions.chamberLength.Layout.Row = row;
            gui.UIControls.workingConditions.chamberLength.Layout.Column = 2;
            chamberLengthHint = uilabel(grid, 'Text', ''); % Hint
            chamberLengthHint.Layout.Row = row;
            chamberLengthHint.Layout.Column = 3;

            row = row + 1;
            % Chamber Radius [cite: 3]
            chamberRadiusLabel = uilabel(grid, 'Text', 'Chamber Radius (m):');
            chamberRadiusLabel.Layout.Row = row;
            chamberRadiusLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.chamberRadius = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.chamberRadius', evt.Value));
            gui.UIControls.workingConditions.chamberRadius.Layout.Row = row;
            gui.UIControls.workingConditions.chamberRadius.Layout.Column = 2;
            chamberRadiusHint = uilabel(grid, 'Text', ''); % Hint
            chamberRadiusHint.Layout.Row = row;
            chamberRadiusHint.Layout.Column = 3;

            row = row + 1;
            % Total SCCM Inflow [cite: 3]
            totalSccmInFlowLabel = uilabel(grid, 'Text', 'Total SCCM Inflow (sccm):');
            totalSccmInFlowLabel.Layout.Row = row;
            totalSccmInFlowLabel.Layout.Column = 1;
            gui.UIControls.workingConditions.totalSccmInFlow = uieditfield(grid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.totalSccmInFlow', evt.Value));
            gui.UIControls.workingConditions.totalSccmInFlow.Layout.Row = row;
            gui.UIControls.workingConditions.totalSccmInFlow.Layout.Column = 2;
            totalSccmInFlowHint = uilabel(grid, 'Text', ''); % Hint
            totalSccmInFlowHint.Layout.Row = row;
            totalSccmInFlowHint.Layout.Column = 3;

            row = row + 1;
            % Total SCCM Outflow [cite: 3]
            totalSccmOutFlowLabel = uilabel(grid, 'Text', 'Total SCCM Outflow:');
            totalSccmOutFlowLabel.Layout.Row = row;
            totalSccmOutFlowLabel.Layout.Column = 1;
            
            % Create a container for the dropdown and optional text field
            flowContainer = uigridlayout(grid, [1, 2]);
            flowContainer.Layout.Row = row;
            flowContainer.Layout.Column = 2;
            flowContainer.ColumnWidth = {'fit', '1x'};
            flowContainer.Padding = [0 0 0 0];
            
            gui.UIControls.workingConditions.totalSccmOutFlowType = uidropdown(flowContainer, ...
                'Items', {'ensureIsobaric', 'totalSccmInFlow', 'Number'}, ...
                'Value', 'ensureIsobaric', ...
                'ValueChangedFcn', @(src, evt) gui.handleSccmOutFlowTypeChange(src, evt.Value));
            gui.UIControls.workingConditions.totalSccmOutFlowType.Layout.Column = 1;
            
            gui.UIControls.workingConditions.totalSccmOutFlow = uieditfield(flowContainer, 'numeric', ...
                'Limits', [0, Inf], ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.totalSccmOutFlow', evt.Value), ...
                'Enable', 'off');
            gui.UIControls.workingConditions.totalSccmOutFlow.Layout.Column = 2;
            
            totalSccmOutFlowHint = uilabel(grid, 'Text', ''); % Hint
            totalSccmOutFlowHint.Layout.Row = row;
            totalSccmOutFlowHint.Layout.Column = 3;

            row = row + 1;
            % Discharge Current [cite: 3] - Optional
            dischargeCurrentCheckbox = uicheckbox(grid, 'Text', 'Discharge Current (A):', ...
                'ValueChangedFcn', @(src, evt) gui.toggleDischargeCurrentEnable(evt.Value));
            dischargeCurrentCheckbox.Layout.Row = row;
            dischargeCurrentCheckbox.Layout.Column = 1;
            gui.UIControls.workingConditions.dischargeCurrentCheckbox = dischargeCurrentCheckbox;
            
            gui.UIControls.workingConditions.dischargeCurrent = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.dischargeCurrent', evt.Value));
            gui.UIControls.workingConditions.dischargeCurrent.Layout.Row = row;
            gui.UIControls.workingConditions.dischargeCurrent.Layout.Column = 2;
            dischargeCurrentHint = uilabel(grid, 'Text', ''); % Hint
            dischargeCurrentHint.Layout.Row = row;
            dischargeCurrentHint.Layout.Column = 3;

            row = row + 1;
            % Discharge Power Density [cite: 3] - Optional
            dischargePowerCheckbox = uicheckbox(grid, 'Text', 'Discharge Power Density (W/m³):', ...
                'ValueChangedFcn', @(src, evt) gui.toggleDischargePowerEnable(evt.Value));
            dischargePowerCheckbox.Layout.Row = row;
            dischargePowerCheckbox.Layout.Column = 1;
            gui.UIControls.workingConditions.dischargePowerCheckbox = dischargePowerCheckbox;
            
            gui.UIControls.workingConditions.dischargePowerDensity = uieditfield(grid, 'numeric', ...
                'Limits', [0, Inf], ...
                'Enable', 'off', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'workingConditions.dischargePowerDensity', evt.Value));
            gui.UIControls.workingConditions.dischargePowerDensity.Layout.Row = row;
            gui.UIControls.workingConditions.dischargePowerDensity.Layout.Column = 2;
            dischargePowerHint = uilabel(grid, 'Text', ''); % Hint
            dischargePowerHint.Layout.Row = row;
            dischargePowerHint.Layout.Column = 3;
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
            gui.UIControls.electronKinetics.isOn = uicheckbox(generalGrid, 'Text', 'Enable Electron Kinetics', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.isOn', evt.Value));
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
            shapeParamLabel = uilabel(generalGrid, 'Text', 'Shape Parameter:');
            shapeParamLabel.Layout.Row = row;
            shapeParamLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.shapeParameter = uidropdown(generalGrid, 'Items', {'Maxwellian', 'Druyvesteyn'}, ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.shapeParameter', evt.Value), ...
                'Visible', 'off');
            gui.UIControls.electronKinetics.shapeParameter.Layout.Row = row;
            gui.UIControls.electronKinetics.shapeParameter.Layout.Column = 2;

            row = row + 1;
            % Ionization Operator Type [cite: 5]
            ionizationOperatorLabel = uilabel(generalGrid, 'Text', 'Ionization Operator:');
            ionizationOperatorLabel.Layout.Row = row;
            ionizationOperatorLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.ionizationOperatorType = uidropdown(generalGrid, 'Items', {'conservative', 'oneTakesAll', 'equalSharing', 'usingSDCS'}, 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.ionizationOperatorType', evt.Value));
            gui.UIControls.electronKinetics.ionizationOperatorType.Layout.Row = row;
            gui.UIControls.electronKinetics.ionizationOperatorType.Layout.Column = 2;

            row = row + 1;
            % Growth Model Type [cite: 5]
            growthModelLabel = uilabel(generalGrid, 'Text', 'Growth Model:');
            growthModelLabel.Layout.Row = row;
            growthModelLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.growthModelType = uidropdown(generalGrid, 'Items', {'temporal', 'spatial'}, 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.growthModelType', evt.Value));
            gui.UIControls.electronKinetics.growthModelType.Layout.Row = row;
            gui.UIControls.electronKinetics.growthModelType.Layout.Column = 2;

            row = row + 1;
            % Include e-e Collisions [cite: 5]
            gui.UIControls.electronKinetics.includeEECollisions = uicheckbox(generalGrid, 'Text', 'Include e-e Collisions', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.includeEECollisions', evt.Value));
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
            
            gui.UIControls.electronKinetics.CARgases = uilistbox(carGrid, 'Multiselect', 'on', 'Enable', 'off');
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
            gasPropsGrid.RowHeight = repmat({'fit'}, 1, 6);
            gasPropsGrid.Padding = [10 10 10 10];
            gasPropsGrid.RowSpacing = 5;
            gasPropsGrid.ColumnSpacing = 10;

            % Mass File
            row = 1;
            massFileLabel = uilabel(gasPropsGrid, 'Text', 'Mass File:');
            massFileLabel.Layout.Row = row;
            massFileLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.mass = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.mass', evt.Value));
            gui.UIControls.gasProperties.mass.Layout.Row = row;
            gui.UIControls.gasProperties.mass.Layout.Column = 2;
            massFileBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.mass'));
            massFileBrowseButton.Layout.Row = row;
            massFileBrowseButton.Layout.Column = 3;

            row = row + 1;
            % Harmonic Frequency File
            harmonicFreqLabel = uilabel(gasPropsGrid, 'Text', 'Harmonic Frequency File:');
            harmonicFreqLabel.Layout.Row = row;
            harmonicFreqLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.harmonicFrequency = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.harmonicFrequency', evt.Value));
            gui.UIControls.gasProperties.harmonicFrequency.Layout.Row = row;
            gui.UIControls.gasProperties.harmonicFrequency.Layout.Column = 2;
            harmonicFreqBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.harmonicFrequency'));
            harmonicFreqBrowseButton.Layout.Row = row;
            harmonicFreqBrowseButton.Layout.Column = 3;

            row = row + 1;
            % Anharmonic Frequency File
            anharmonicFreqLabel = uilabel(gasPropsGrid, 'Text', 'Anharmonic Frequency File:');
            anharmonicFreqLabel.Layout.Row = row;
            anharmonicFreqLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.anharmonicFrequency = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.anharmonicFrequency', evt.Value));
            gui.UIControls.gasProperties.anharmonicFrequency.Layout.Row = row;
            gui.UIControls.gasProperties.anharmonicFrequency.Layout.Column = 2;
            anharmonicFreqBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.anharmonicFrequency'));
            anharmonicFreqBrowseButton.Layout.Row = row;
            anharmonicFreqBrowseButton.Layout.Column = 3;

            row = row + 1;
            % Rotational Constant File
            rotationalConstantLabel = uilabel(gasPropsGrid, 'Text', 'Rotational Constant File:');
            rotationalConstantLabel.Layout.Row = row;
            rotationalConstantLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.rotationalConstant = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.rotationalConstant', evt.Value));
            gui.UIControls.gasProperties.rotationalConstant.Layout.Row = row;
            gui.UIControls.gasProperties.rotationalConstant.Layout.Column = 2;
            rotationalConstantBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.rotationalConstant'));
            rotationalConstantBrowseButton.Layout.Row = row;
            rotationalConstantBrowseButton.Layout.Column = 3;

            row = row + 1;
            % Electric Quadrupole Moment File
            electricQuadrupoleMomentLabel = uilabel(gasPropsGrid, 'Text', 'Electric Quadrupole Moment File:');
            electricQuadrupoleMomentLabel.Layout.Row = row;
            electricQuadrupoleMomentLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.electricQuadrupoleMoment = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.electricQuadrupoleMoment', evt.Value));
            gui.UIControls.gasProperties.electricQuadrupoleMoment.Layout.Row = row;
            gui.UIControls.gasProperties.electricQuadrupoleMoment.Layout.Column = 2;
            electricQuadrupoleMomentBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.electricQuadrupoleMoment'));
            electricQuadrupoleMomentBrowseButton.Layout.Row = row;
            electricQuadrupoleMomentBrowseButton.Layout.Column = 3;

            row = row + 1;
            % OPB Parameter File
            opbParameterLabel = uilabel(gasPropsGrid, 'Text', 'OPB Parameter File:');
            opbParameterLabel.Layout.Row = row;
            opbParameterLabel.Layout.Column = 1;
            gui.UIControls.gasProperties.OPBParameter = uieditfield(gasPropsGrid, 'text', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.gasProperties.OPBParameter', evt.Value));
            gui.UIControls.gasProperties.OPBParameter.Layout.Row = row;
            gui.UIControls.gasProperties.OPBParameter.Layout.Column = 2;
            opbParameterBrowseButton = uibutton(gasPropsGrid, 'Text', 'Browse...', 'ButtonPushedFcn', @(src,evt) gui.browseFile('electronKinetics.gasProperties.OPBParameter'));
            opbParameterBrowseButton.Layout.Row = row;
            opbParameterBrowseButton.Layout.Column = 3;

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
                'Items', {});
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
                'Items', {});
            gui.UIControls.electronKinetics.stateProperties.energy.Layout.Row = row;
            gui.UIControls.electronKinetics.stateProperties.energy.Layout.Column = 2;
            energyBtnGrid = uigridlayout(grid, [2, 1]);
            energyBtnGrid.Layout.Row = row;
            energyBtnGrid.Layout.Column = 3;
            energyBtnGrid.RowHeight = {'fit', 'fit'};
            energyBtnGrid.Padding = [0 0 0 0];
            energyAddButton = uibutton(energyBtnGrid, 'Text', 'Add', 'ButtonPushedFcn', @(src,evt) gui.addListItem('electronKinetics.stateProperties.energy', false));
            energyRemoveButton = uibutton(energyBtnGrid, 'Text', 'Remove', 'ButtonPushedFcn', @(src,evt) gui.removeListItem('electronKinetics.stateProperties.energy'));

            % Statistical Weight listbox
            row = row + 1;
            statisticalWeightLabel = uilabel(grid, 'Text', 'Statistical Weight:');
            statisticalWeightLabel.Layout.Row = row;
            statisticalWeightLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.stateProperties.statisticalWeight = uilistbox(grid, ...
                'Multiselect', 'on', ...
                'Items', {});
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
                'Items', {});
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
            % Variable Energy Grid (Coming Soon) - Mockup checkbox
            variableGridCheckbox = uicheckbox(energyGridGrid, 'Text', 'Variable Energy Grid (Soon...)', 'Value', false, 'Enable', 'off');
            variableGridCheckbox.Layout.Row = row;
            variableGridCheckbox.Layout.Column = [1, 2];

            row = row + 1;
            % Max Energy [cite: 10]
            maxEnergyLabel = uilabel(energyGridGrid, 'Text', 'Max Energy (eV):');
            maxEnergyLabel.Layout.Row = row;
            maxEnergyLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.maxEnergy', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.maxEnergy.Layout.Column = 2;

            row = row + 1;
            % Cell Number [cite: 10]
            cellNumberLabel = uilabel(energyGridGrid, 'Text', 'Energy Cell Number:');
            cellNumberLabel.Layout.Row = row;
            cellNumberLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.cellNumber = uieditfield(energyGridGrid, 'numeric', 'Limits', [1, Inf], 'ValueDisplayFormat', '%.0f', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.cellNumber', evt.Value));
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
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueDisplayFormat', '%.0f', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.minEedfDecay.Layout.Column = 2;

            row = row + 1;
            % Max EEDF Decay
            maxEedfDecayLabel = uilabel(energyGridGrid, 'Text', 'Max EEDF Decay:');
            maxEedfDecayLabel.Layout.Row = row;
            maxEedfDecayLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueDisplayFormat', '%.0f', 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.maxEedfDecay.Layout.Column = 2;

            row = row + 1;
            % Update Factor
            updateFactorLabel = uilabel(energyGridGrid, 'Text', 'Update Factor:');
            updateFactorLabel.Layout.Row = row;
            updateFactorLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.energyGrid.smartGrid.updateFactor', evt.Value));
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Layout.Row = row;
            gui.UIControls.electronKinetics.numerics.energyGrid.smartGrid.updateFactor.Layout.Column = 2;

            row = row + 1;
            % Max Power Balance Rel Error [cite: 12]
            maxPowerBalanceLabel = uilabel(energyGridGrid, 'Text', 'Max Power Bal. Rel. Error:');
            maxPowerBalanceLabel.Layout.Row = row;
            maxPowerBalanceLabel.Layout.Column = 1;
            gui.UIControls.electronKinetics.numerics.maxPowerBalanceRelError = uieditfield(energyGridGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.maxPowerBalanceRelError', evt.Value));
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
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.algorithm = uidropdown(nonLinearGrid, 'Items', {'mixingDirectSolutions', 'temporalIntegration'}, 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.algorithm', evt.Value));

            row = row + 1;
            % Mixing Parameter [cite: 12] - Simplified slider with precision input
            uilabel(nonLinearGrid, 'Text', 'Mixing Parameter:');
            mixingParamGrid = uigridlayout(nonLinearGrid, [1, 2]);
            mixingParamGrid.ColumnWidth = {'0.6x', '0.4x'}; % Give more space to text field
            mixingParamGrid.Padding = [0 0 0 0];
            mixingParamGrid.ColumnSpacing = 5;
            
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameter = uislider(mixingParamGrid, 'Limits', [0, 1], 'Value', 0.7, 'MajorTicks', [0, 1], 'MajorTickLabels', {'0', '1'}, 'ValueChangedFcn', @(src, evt) gui.sliderMixingParameterChanged(src, evt.Value));
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.mixingParameterPrecise = uieditfield(mixingParamGrid, 'numeric', 'Limits', [0, 1], 'Value', 0.7, 'ValueChangedFcn', @(src, evt) gui.textMixingParameterChanged(src, evt.Value));
            
            mixingParamGrid.Layout.Row = row;
            mixingParamGrid.Layout.Column = 2;

            row = row + 1;
            % Max EEDF Rel Error [cite: 12]
            uilabel(nonLinearGrid, 'Text', 'Max EEDF Rel. Error:');
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.maxEedfRelError = uieditfield(nonLinearGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.maxEedfRelError', evt.Value));

            row = row + 1;
            % Advanced Parameters Section (Expandable) - Only visible for temporalIntegration
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
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol', evt.Value));
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Layout.Row = 2;
            gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Layout.Column = 2;

                         % Relative Tolerance
             uilabel(contentGrid, 'Text', 'Relative Tolerance:');
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol', evt.Value));
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Layout.Row = 3;
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Layout.Column = 2;

             % Max Step
             uilabel(contentGrid, 'Text', 'Max Step:');
             gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep = uieditfield(contentGrid, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(src, evt) gui.updateField(src, 'electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep', evt.Value));
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
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'output.isOn', evt.Value));
            gui.UIControls.output.isOn.Layout.Row = 5;
            gui.UIControls.output.isOn.Layout.Column = 1;

            % Row 6: Data Format
            dataFormatLabel = uilabel(grid, 'Text', 'Data Format:');
            dataFormatLabel.Layout.Row = 6;
            dataFormatLabel.Layout.Column = 1;
            
            gui.UIControls.output.dataFormat = uidropdown(grid, ...
                'Items', {'txt', 'hdf5', 'hdf5+txt'}, ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'output.dataFormat', evt.Value));
            gui.UIControls.output.dataFormat.Layout.Row = 6;
            gui.UIControls.output.dataFormat.Layout.Column = [2, 3];

            % Row 7: Output Folder
            folderLabel = uilabel(grid, 'Text', 'Output Folder:');
            folderLabel.Layout.Row = 7;
            folderLabel.Layout.Column = 1;
            gui.UIControls.output.folder = uieditfield(grid, 'text', ...
                'ValueChangedFcn', @(src, evt) gui.updateField(src, 'output.folder', evt.Value));
            gui.UIControls.output.folder.Layout.Row = 7;
            gui.UIControls.output.folder.Layout.Column = 2;
            
            browseButton = uibutton(grid, 'Text', 'Browse...', ...
                'ButtonPushedFcn', @gui.browseFolder);
            browseButton.Layout.Row = 7;
            browseButton.Layout.Column = 3;

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
                        'rateCoefficients', 'powerBalance', 'lookUpTable'};
            
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
                                    
                                    dataPath = sprintf('%s.%s.%s.%s', sectionName, subSectionName, subSubSectionName, subSubControlName);
                                    gui.setControlValue(subSubControl, dataPath);
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
                            dataSetNames = {'inputs', 'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTable'};
                            
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
            try
                advancedButton = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton;
                if strcmp(gui.Setup.electronKinetics.numerics.nonLinearRoutines.algorithm, 'temporalIntegration')
                    advancedButton.Visible = 'on';
                else
                    advancedButton.Visible = 'off';
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

            % Configure initial state of SCCM Outflow type
            try
                % Set initial dropdown value based on the current setup value
                currentValue = gui.Setup.workingConditions.totalSccmOutFlow;
                if isnumeric(currentValue)
                    gui.UIControls.workingConditions.totalSccmOutFlowType.Value = 'Number';
                    gui.handleSccmOutFlowTypeChange([], 'Number');
                    % Set the numeric field value
                    gui.UIControls.workingConditions.totalSccmOutFlow.Value = currentValue;
                elseif ischar(currentValue) && ~isempty(str2num(currentValue))
                    gui.UIControls.workingConditions.totalSccmOutFlowType.Value = 'Number';
                    gui.handleSccmOutFlowTypeChange([], 'Number');
                    % Convert string to numeric
                    gui.UIControls.workingConditions.totalSccmOutFlow.Value = str2double(currentValue);
                elseif strcmp(currentValue, 'totalSccmInFlow')
                    gui.UIControls.workingConditions.totalSccmOutFlowType.Value = 'totalSccmInFlow';
                    gui.handleSccmOutFlowTypeChange([], 'totalSccmInFlow');
                else
                    gui.UIControls.workingConditions.totalSccmOutFlowType.Value = 'ensureIsobaric';
                    gui.handleSccmOutFlowTypeChange([], 'ensureIsobaric');
                end
            catch ME
                warning('Error configuring SCCM Outflow type: %s', ME.message);
            end
            
            % Configure data sets checkboxes after all other controls are populated
            try
                if isfield(gui.UIControls, 'output') && isfield(gui.UIControls.output, 'dataSets')
                    if isfield(gui.Setup, 'output') && isfield(gui.Setup.output, 'dataSets')
                        selectedDataSets = gui.Setup.output.dataSets;
                        dataSetNames = {'inputs', 'log', 'eedf', 'swarmParameters', 'rateCoefficients', 'powerBalance', 'lookUpTable'};
                        
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
        end

        function setControlValue(gui, control, dataPath)
            % Helper function to get value from Setup and set control value
            
            % Skip output.dataSets fields completely to avoid warnings
            if contains(dataPath, 'output.dataSets.')
                return;
            end
            
            % Skip UI-only controls that don't exist in Setup
            if contains(dataPath, 'totalSccmOutFlowType') || ...
               contains(dataPath, 'dischargeCurrentCheckbox') || ...
               contains(dataPath, 'dischargePowerCheckbox') || ...
               contains(dataPath, 'LXCatExtraCheckbox') || ...
               contains(dataPath, 'effectivePopCheckbox') || ...
               contains(dataPath, 'CARcheckbox') || ...
               contains(dataPath, 'CARaddButton') || ...
               contains(dataPath, 'CARremoveButton')
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
            elseif strcmp(fieldPath, 'electronKinetics.numerics.nonLinearRoutines.algorithm')
                % Show/hide advanced button based on algorithm selection
                advancedButton = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedButton;
                if strcmp(value, 'temporalIntegration')
                    advancedButton.Visible = 'on';
                else
                    advancedButton.Visible = 'off';
                    % Also hide the advanced panel if algorithm changes
                    advancedPanel = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.advancedPanel;
                    advancedPanel.Visible = 'off';
                    advancedButton.Text = 'Advanced (Optional) ▼';
                end
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
                    gui.Fig.WindowState = oldState;
                    return; % User cancelled
                end
                
                % Restore window state after file dialog
                drawnow;
                pause(0.1);
                gui.Fig.WindowState = 'normal';
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
                prompt = {['Enter new item for ', fieldPath, ':']};
                dlgtitle = 'Add List Item';
                dims = [1 50];
                answer = inputdlg(prompt, dlgtitle, dims);
                if isempty(answer)
                    return; % User cancelled
                end
                newItem = answer{1};
            end

            if ~isempty(newItem) && ~ismember(newItem, currentItems)
                % For gas fractions, validate the format and check if sum equals 1
                if endsWith(fieldPath, 'gasProperties.fraction') || contains(fieldPath, 'gasProperties.fraction')
                    % Validate format: should be "species = number"
                    pattern = '^\s*(\w+)\s*=\s*([0-9]+\.?[0-9]*)\s*$';
                    match = regexp(newItem, pattern, 'tokens');
                    if isempty(match)
                        uialert(gui.Fig, ['Invalid format for gas fraction. Expected format: "species = number"\n' ...
                            'Example: "N2 = 1" or "H2 = 0.5"'], 'Invalid Format');
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
                    
                    % Check if sum equals 1 (with tolerance for floating point)
                    if abs(totalSum - 1.0) > 1e-6
                        warningMsg = sprintf(['Gas fractions must sum to 1.0.\n' ...
                            'Current sum: %.6f\n' ...
                            'Continue anyway?'], totalSum);
                        % Just warn but don't prevent, as user might be entering multiple values
                        % For now, we'll allow it but show a warning
                        uialert(gui.Fig, warningMsg, 'Fraction Sum Warning', 'Icon', 'warning');
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
                
                % Check if sum equals 1 (with tolerance for floating point)
                if abs(totalSum - 1.0) > 1e-6
                    warningMsg = sprintf(['After removal, gas fractions sum to %.6f instead of 1.0.\n' ...
                        'Please ensure all fractions sum to 1.0.'], totalSum);
                    uialert(gui.Fig, warningMsg, 'Fraction Sum Warning', 'Icon', 'warning');
                end
            end
            
            % Update the Setup struct as well
            gui.setNestedField(fieldPath, currentItems);
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
                gui.Fig.WindowState = oldState;
                return;
            end
            
            % Restore window state after file dialog
            drawnow;
            pause(0.1);
            gui.Fig.WindowState = 'normal';
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
            
            control = gui.getNestedField(gui.UIControls, fieldPath);
            control.Value = filePath;
            gui.setNestedField(fieldPath, filePath); % Update setup directly
        end

        function handleSccmOutFlowTypeChange(gui, ~, value)
            % Handle changes to the SCCM Outflow type dropdown
            numField = gui.UIControls.workingConditions.totalSccmOutFlow;
            
            if strcmp(value, 'Number')
                % Enable numeric field for direct numeric input
                numField.Enable = 'on';
                numField.Value = 1; % Set default value for numeric input
            else
                % Disable numeric field
                numField.Enable = 'off';
                if strcmp(value, 'totalSccmInFlow') || strcmp(value, 'ensureIsobaric')
                    % Update the setup struct with the selected value
                    gui.setNestedField('workingConditions.totalSccmOutFlow', value);
                end
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
            if strcmp(eedfType, 'prescribedEedf')
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'on';
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
        end

        function handleEedfTypeChangeForInit(gui, eedfType)
            % Show/hide shape parameter based on EEDF type (for initialization)
            if strcmp(eedfType, 'prescribedEedf')
                gui.UIControls.electronKinetics.shapeParameter.Visible = 'on';
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

        function loadSettings(gui, ~, ~)
            % Placeholder for loading settings from a file
            [file, path] = uigetfile('*.txt;*.setup', 'Load LoKI Setup File');
            if isequal(file, 0) || isequal(path, 0)
                return; % User cancelled
            end
            inputFile = fullfile(path, file);
            try
                % gui.Setup = gui.parseInputFile(inputFile); % Implement this
                gui.populateGUIFromSetup(); % Update UI
                uialert(gui.Fig, ['Settings loaded from ', file], 'Load Successful');
            catch ME
                uialert(gui.Fig, ['Error loading file: ', ME.message], 'Load Error');
            end
        end

        function saveInputFile(gui, ~, ~)
            % Save the current setup to a file
            defaultName = ['LoKI_Setup_', datestr(now,'yyyymmdd_HHMMSS'), '.txt'];
            startPath = gui.Setup.output.folder; % Suggest saving in output folder
            if ~isfolder(startPath)
                startPath = pwd; % Fallback to current directory
            end
            [file, path] = uiputfile('*.txt', 'Save LoKI Input File As', fullfile(startPath, defaultName));

            if isequal(file, 0) || isequal(path, 0)
                uialert(gui.Fig, 'File save cancelled.', 'Cancelled');
                return; % User cancelled
            end

            outputFile = fullfile(path, file);
            try
                gui.generateInputFile(outputFile);
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
            fprintf(fid, '%% LoKI Input File generated by LoKI_GUI on %s %%\n\n', datestr(now));

            % Update totalSccmOutFlow from UI before generating file
            try
                dropdownValue = gui.UIControls.workingConditions.totalSccmOutFlowType.Value;
                if strcmp(dropdownValue, 'Number')
                    numValue = gui.UIControls.workingConditions.totalSccmOutFlow.Value;
                    if isnumeric(numValue)
                        gui.Setup.workingConditions.totalSccmOutFlow = numValue;
                    end
                else
                    gui.Setup.workingConditions.totalSccmOutFlow = dropdownValue;
                end
            catch ME
                warning('Error updating totalSccmOutFlow from UI: %s', ME.message);
            end

            % CARgases will be populated from the listbox items automatically when writing

            % Use recursion or explicit handling for nested structs
            sections = fieldnames(gui.Setup);
            for i = 1:length(sections)
                sectionName = sections{i};
                if strcmp(sectionName, 'electronKinetics')
                    % Special handling for electronKinetics to reorder fields
                    fprintf(fid, '%s:\n', sectionName); % Section header
                    gui.writeElectronKineticsContent(gui.Setup.(sectionName), '  ', fid);
                else
                    fprintf(fid, '%s:\n', sectionName); % Section header
                    gui.writeStructContent(gui.Setup.(sectionName), '  ', fid); % Indent content
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
                try
                    if strcmp(gui.UIControls.electronKinetics.eedfType.Value, 'prescribedEedf')
                        % Convert shape parameter to numeric (1 for Maxwellian, 2 for Druyvesteyn)
                        shapeValue = 1; % default to Maxwellian
                        if strcmp(ek.shapeParameter, 'Druyvesteyn')
                            shapeValue = 2;
                        end
                        fprintf(fid, '%s shapeParameter: %d\n', indent, shapeValue);
                    end
                catch
                    % Skip if not prescribedEedf
                end
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
            
            % LXCat files
            if isfield(ek, 'LXCatFiles')
                fprintf(fid, '%sLXCatFiles:\n', indent);
                for j = 1:length(ek.LXCatFiles)
                    fprintf(fid, '%s  - %s\n', indent, ek.LXCatFiles{j});
                end
            end
            
            % LXCat Extra Files (only if checked)
            if isfield(ek, 'LXCatExtraFiles')
                try
                    if gui.UIControls.electronKinetics.LXCatExtraCheckbox.Value
                        fprintf(fid, '%s LXCatExtraFiles:\n', indent);
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
                        fprintf(fid, '%s effectiveCrossSectionPopulations:\n', indent);
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
                        fprintf(fid, '%s CARgases:\n', indent);
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
            
            % Gas Properties (write before stateProperties)
            if isfield(ek, 'gasProperties')
                fprintf(fid, '%sgasProperties:\n', indent);
                % Write gas properties manually to handle fraction from UI
                fprintf(fid, '%s  mass: %s\n', indent, ek.gasProperties.mass);
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
                fprintf(fid, '%s  harmonicFrequency: %s\n', indent, ek.gasProperties.harmonicFrequency);
                fprintf(fid, '%s  anharmonicFrequency: %s\n', indent, ek.gasProperties.anharmonicFrequency);
                fprintf(fid, '%s  rotationalConstant: %s\n', indent, ek.gasProperties.rotationalConstant);
                fprintf(fid, '%s  electricQuadrupoleMoment: %s\n', indent, ek.gasProperties.electricQuadrupoleMoment);
                fprintf(fid, '%s  OPBParameter: %s\n', indent, ek.gasProperties.OPBParameter);
            end
            
            % State Properties (write normally with writeStructContent)
            if isfield(ek, 'stateProperties')
                fprintf(fid, '%sstateProperties:\n', indent);
                gui.writeStructContent(ek.stateProperties, [indent, '  '], fid);
            end
            
            % Numerics (write normally with writeStructContent)
            if isfield(ek, 'numerics')
                fprintf(fid, '%snumerics:\n', indent);
                gui.writeStructContent(ek.numerics, [indent, '  '], fid);
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
                        % Check if ODE parameters should be written by checking the checkbox state
                        if gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.isOn.Value
                            % If checkbox is checked, write the section with current values
                            fprintf(fid, '%s%s:\n', indent, fieldName);
                            % Write the actual values from the GUI controls with extra indentation
                            absTol = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.AbsTol.Value;
                            relTol = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.RelTol.Value;
                            maxStep = gui.UIControls.electronKinetics.numerics.nonLinearRoutines.odeSetParameters.MaxStep.Value;
                            fprintf(fid, '%s   AbsTol: %s\n', indent, num2str(absTol));
                            fprintf(fid, '%s   RelTol: %s\n', indent, num2str(relTol));
                            fprintf(fid, '%s   MaxStep: %s\n', indent, num2str(maxStep));
                        end
                        % If not enabled, don't write anything
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
                        % If none of the above, don't write empty fields
                    else
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
                    elseif strcmp(fieldName, 'dischargePowerDensity')
                        try
                            if gui.UIControls.workingConditions.dischargePowerCheckbox.Value
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

        function runSimulation(gui, ~, ~)
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
                lokibcl(tempInputFile(7:end)); % Placeholder for success
                % Run the simulation
                
            % catch ME
            %     uialert(gui.Fig, ['Error occurred during simulation: ', ME.message], 'Simulation Runtime Error');
            % end

            % Optional: Clean up temporary file
            % delete(tempInputFile);
            fprintf('--- Simulation Finished ---\n');
         end
    end % methods
end % classdef