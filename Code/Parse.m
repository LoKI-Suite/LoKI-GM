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

classdef Parse
  
  properties (Constant = true)
    
    inputFolder = 'Input';
    commentChar = '%';
    wildCardChar = '*';
    stateRegExp = ['\s*(?<quantity>\d+\s*)?(?<gasName>\w+)\((?<ionCharg>[+-])?'...
      '(?:,)?(?<eleLevel>[\w''.\[\]/+*^|-]+){1}(?:,v=)?(?<vibLevel>(?<=,v=)[\w+*|-]+)?'...
      '(?:,J=)?(?<rotLevel>(?<=,J=)[\d+*|-]+)?\)\s*'];
    electronRegExp = '\s*(?<quantity>\d+)?((?:[eE]\s*)(?=[\s]|[+][^\(])|(?:[eE]$))';
    
  end
  
  methods (Static)
    
    function [setupStruct, unparsed] = setupFile(fileName)
    % setupFile Reads the configuration file for the simulation and returns
    % a setup structure with all the information.
      
      [structArray, unparsed] = file2structArray([Parse.inputFolder filesep fileName]);
      setupStruct = structArray2struct(structArray);
      
    end

    function [setupStruct, unparsed] = setupFileJson(fileName)
    % setupFileJson Reads the configuration file in JSON format for the 
    % simulation and returns a setup structure with all the information.
     
        % open file
        file = [Parse.inputFolder filesep fileName];
        fileID = fopen(file, 'r');
        if fileID == -1
            error('\t Unable to open file: %s\n', file);
        end   
    
        fileStr = fileread(file); % reads file as text
        
        % Use regex to remove key-value pairs with keys starting with '_'
        % after trimming whitespaces
        fileTextCleaned = regexprep(fileStr, '\s*"_\w+":\s*[^,}]*[,]?', ''); 
        % Remove array elements starting with '_'
        fileTextCleaned = regexprep(fileTextCleaned, '\s*"_[^"]*",?', '');

        data = jsondecode(fileTextCleaned); % Using the jsondecode function to parse JSON from string
        
        unparsed = file2structArrayJson(fileID);

        % Transpose all list elements in the data structure
        setupStruct = transposeFields(data);

        fclose(fileID);
      
    end
    
    function LXCatEntryArray = LXCatFiles(fileName)
    % LXCatFiles Reads LXCat files, parse their content and returns an
    % structure array 'LXCatEntryArray' with all the information.
      
      % definition of regular expressions to parse LXCat file
      LXCatRegExp1 = 'PROCESS: (?<reactants>.+?)(?<direction>->|<->) (?<products>.+?), (?<type>\w+)';
      LXCatRegExp2 = 'E = (?<threshold>[\d.Ee+-]+) eV';
      LXCatRegExp3 = ['\[(?<reactants>.+?)(?<direction>->|<->)(?<products>.+?), (?<type>\w+)' ...
        '(?<subtype>, momentum-transfer)?\]'];
      % create a cell array with filenames in case only one file is
      % received as input
      if ischar(fileName)
        fileName = {fileName};
      end
      % create an empty struct array of LXCat entries 
      LXCatEntryArray = struct.empty;
      % loop over different LXCat files that have to be read
      for i = 1:length(fileName)
        % open LXCat file
        fileID = fopen([Parse.inputFolder filesep fileName{i}], 'r');
        if fileID<0
          error(' Unable to open LXCat file: %s\n', fileName{i});
        end
        % parse LXCat file
        while ~feof(fileID)
          description = regexp(fgetl(fileID), LXCatRegExp1, 'names', 'once');
          if ~isempty(description)
            parameter = regexp(fgetl(fileID), LXCatRegExp2, 'names', 'once');
            description = regexp(fgetl(fileID),LXCatRegExp3, 'names', 'once');
            while ~strncmp(fgetl(fileID), '-----', 5)
            end
            rawCrossSection = (fscanf(fileID, '%f', [2 inf]));
            % add LXCat entry information into LXCatEntry struct array
            LXCatEntryArray = addLXCatEntry(description, parameter, rawCrossSection, LXCatEntryArray);
          end
        end
        % close LXCat file
        fclose(fileID);
      end

    end

    function LXCatEntryArray_new = LXCatFilesJson(fileName)
    % LXCatFilesJson Reads LXCat files in JSON format, parse their content 
    % and returns a structure array 'LXCatEntryArray' with all the 
    % information.
      
      % create a cell array with filenames in case only one file is
      % received as input
      if ischar(fileName)
        fileName = {fileName};
      end
      % create an empty struct array of LXCat entries 
      LXCatEntryArray_new = struct.empty;
      % loop over different LXCat files that have to be read
      for i = 1:length(fileName)
        % open LXCat file
        fileID = fopen([Parse.inputFolder filesep fileName{i}], 'r');
        if fileID<0
          error(' Unable to open LXCat file: %s\n', fileName{i});
        end
        % parse LXCat file
        fileStr = fileread([Parse.inputFolder filesep fileName{i}]); % reads file 
        jsonData = jsondecode(fileStr); % Using the jsondecode function to parse JSON from string
        
        % for each LXCat Process
        for p = 1:size(jsonData.processes,1)

            % adapt json data from decoding to be used by the
            % addLXCatEntryJson function
            threshold = jsonData.processes(p,1).info.threshold;
            rawCrossSection = jsonData.processes(p,1).info.data.values';
            description = jsonData.processes(p,1).reaction;
            
            % add LXCat entry information into LXCatEntry struct array
            LXCatEntryArray_new = addLXCatEntryJson(jsonData.states, description, threshold, ...
                rawCrossSection, LXCatEntryArray_new);

        end
        % close LXCat file
        fclose(fileID);
      end

    end

    function gasAndValueArray = gasPropertyFile(fileName)
    % gasPropertyFile Reads a file with gas properties, parse its content
    % and returns an structure array 'gasAndValueArray' with all the
    % information.
    
      regExp = '(?<gasName>\w+)\s*(?<valueStr>.+)';
      fileID = fopen([Parse.inputFolder filesep fileName], 'r');
      if fileID == -1
        error('\t Unable to open file: %s\n', [Parse.inputFolder filesep fileName]);
      end
      
      gasAndValueArray = struct.empty;
      while ~feof(fileID)
        line = removeComments(fgetl(fileID));
        if isempty(line)
          continue
        end
        gasProperty = regexp(line, regExp, 'names', 'once');
        if ~isempty(gasProperty)
          gasAndValueArray(end+1).gasName = gasProperty.gasName;
          gasAndValueArray(end).value = str2num(gasProperty.valueStr);
        end
      end
      
      fclose(fileID);
      
    end
    
    function parsedEntry = gasPropertyEntry(entry)
    % gasPropertyEntry parses an entry (string) of the input file that
    % contains information related to a certain gas property. It returns an
    % structure with the parsed information. 
    
      regExp = ['(?<gasName>\w+)\s*=\s*(?<constant>[\d.eE\s()*/+-]+)?' ...
        '(?<function>\w+)?(?(function)@?)(?<argument>.+)?\s*'];
      parsedEntry = regexp(entry, regExp, 'names');
      if isempty(parsedEntry)
        parsedEntry = struct('fileName', entry);
      elseif ~isempty(parsedEntry.function)
        if ~isempty(parsedEntry.constant)
          error(['Error found when parsing state property entry:\n%s\n' ...
            'Please, fix the problem and run the code again'], entry);
        end
        if isempty(parsedEntry.argument)
          parsedEntry.argument = cell(0);
        else
          parsedEntry.argument = strsplit(parsedEntry.argument, ',');
          for i=1:length(parsedEntry.argument)
            numericArgument = str2num(parsedEntry.argument{i});
            if ~isnan(numericArgument)
              parsedEntry.argument{i} = numericArgument;
            end
          end
        end
      else
        parsedEntry.constant = str2num(parsedEntry.constant);
      end
      
    end
    
    function stateAndValueArray = statePropertyFile(fileName)
    % statePropertyFile Reads a file with state properties, parse its
    % content and returns an structure array 'stateAndValueArray' with all
    % the information.
    
      regExp = [Parse.stateRegExp '\s*(?<valueStr>.+)'];

      fileID = fopen([Parse.inputFolder filesep fileName], 'r');
      if fileID == -1
        error('\t Unable to open file: %s\n', [Parse.inputFolder filesep fileName]);
      end
      
      stateAndValueArray = struct.empty;
      while ~feof(fileID)
        line = removeComments(fgetl(fileID));
        if isempty(line)
          continue
        end
        stateProperty = regexp(line, regExp, 'names', 'once');
        if ~isempty(stateProperty)
          stateAndValueArray(end+1).gasName = stateProperty.gasName;
          stateAndValueArray(end).ionCharg = stateProperty.ionCharg;
          stateAndValueArray(end).eleLevel = stateProperty.eleLevel;
          stateAndValueArray(end).vibLevel = stateProperty.vibLevel;
          stateAndValueArray(end).rotLevel = stateProperty.rotLevel;
          stateAndValueArray(end).value = str2num(stateProperty.valueStr);
        end
      end
      
    end
    
    function parsedEntry = statePropertyEntry(entry)
    % statePropertyEntry parses an entry (string) of the input file that
    % contains information related to a certain state property. It returns
    % an structure with the parsed information. 
    
      regExp = [Parse.stateRegExp '\s*=\s*(?<constant>[\d.eE\s()*/+-]+)?' ...
            '(?<function>\w+)?(?(function)@?)(?<argument>.+)?\s*'];
      parsedEntry = regexp(entry, regExp, 'names');
      if isempty(parsedEntry)
        parsedEntry = struct('fileName', entry);
      elseif ~isempty(parsedEntry.function)
        if ~isempty(parsedEntry.constant)
          error(['Error found when parsing state property entry:\n%s\n' ...
            'Please, fix the problem and run the code again'], entry);
        end
        if isempty(parsedEntry.argument)
          parsedEntry.argument = cell(0);
        else
          parsedEntry.argument = strsplit(parsedEntry.argument, ',');
          for i=1:length(parsedEntry.argument)
            numericArgument = str2num(parsedEntry.argument{i});
            if ~isnan(numericArgument)
              parsedEntry.argument{i} = numericArgument;
            end
          end
        end
      else
        parsedEntry.constant = str2num(parsedEntry.constant);
      end
      
    end
    
  end
  
end

function data = transposeFields(data)
    % find column vectors in structures and transpose to row vectors. Process fields recursively.
    fields = fieldnames(data);  % Get all field names

    for i = 1:numel(fields)
        field = fields{i};
        
        if isstruct(data.(field))  % If the field is a structure, recurse into it
            data.(field) = transposeFields(data.(field));

        else 
          % if the field is not a structure, call str2value to convert the field to a value
          if ischar(data.(field))
            data.(field) = str2value(data.(field));
          end
        
          % if the field is a vector or cell array, transpose to row
          if iscell(data.(field)) || isvector(data.(field))  % Check if it's a vector or cell array
              if size(data.(field), 1) > 1 && size(data.(field), 2) == 1  % If column vector, transpose to row
                  data.(field) = data.(field)';
              end
          end
        end
    end
end

function rawLine = file2structArrayJson(fileID)
% file2structArrayJson Reads an input file and returns the unparsed data 
% for output purposes.
%

  rawLine = cell.empty;
  while ~feof(fileID)
    % erase commas, quotation marks and brackets from the input file
    % match = ["{", "}", "[", "]", '"', ","];
    % match = ["{", "}", "[", "]", '"'];
    % rawLine{end+1} = erase((fgetl(fileID)), match);
    % line = strtrim(rawLine{end});
    line = fgetl(fileID);
    if ~isempty(line)
        rawLine{end+1} = line;
    end
  end
      
end


function [structArray, rawLine] = file2structArray(file)
% file2structArray Reads an input file and and creates an array of input
% structs with all the information.
%
% See also structArray2struct
      
  fileID = fopen(join(file, ''), 'r');
  if fileID == -1
    error('\t Unable to open file: %s\n', file);
  end
  structArray = struct.empty;
  rawLine = cell.empty;
  while ~feof(fileID)
    rawLine{end+1} = removeComments(fgetl(fileID));
    line = strtrim(rawLine{end});
    if isempty(line)
      rawLine = rawLine(1:end-1);
      continue;
    elseif strncmp(line, '-', 1)
      structArray(end).value{end+1} = str2value(line(2:end));
    else
      nameAndValue=strtrim(strsplit(line,':'));
      structArray(end+1).name = nameAndValue{1};
      structArray(end).value = str2value(nameAndValue{2});
      structArray(end).level = regexp(rawLine{end}, '\S', 'once');
    end
  end
  fclose(fileID);
      
end


function structure = structArray2struct(structArray)
% structArray2struct Convert an array of structures into a single structure.
%
% See also file2structArray
    
  i = 1;
  iMax = length(structArray);
  while i<=iMax
    if (i<iMax && structArray(i).level<structArray(i+1).level)
      counter = 1;
      for j = i+2:iMax
        if (structArray(j).level>structArray(i).level)
          counter = counter+1;
        else
          break
        end
      end
      structure.(structArray(i).name) = ...
        structArray2struct(structArray(i+1:i+counter));
      i = i+counter+1;
    else
      structure.(structArray(i).name) = structArray(i).value;
      i = i+1;
    end
  end
      
end

function textLineClean = removeComments(textLine)
% removeComments Return a string without comments.
%
% Comment delimiter is defined in the property commentChar of the Parse class.
% Note: Comment symbol is NOT ignored even if it's inside a string.
  
  idx = regexp(textLine, Parse.commentChar, 'once'); 
  if isempty(idx)
    textLineClean = textLine;
  else
    textLineClean = textLine(1:idx-1);
  end
end

function value = str2value(str)
% str2value Converts a string to a value. Posible values are: numeric,
% logical or string

  [value, is_number_or_logical] = str2num(strtrim(str));
  if ~is_number_or_logical
    value = strtrim(str);
  end
      
end

function LXCatEntryArray = addLXCatEntryJson(jsonStates, description, parameter, rawCrossSection, LXCatEntryArray)
% addLXCatEntryJson analyses the information of a particular LXCat entry in
% JSON format and adds it to the structure array LXCatEntryArray
  
  % The reading of the typeTags is a WIP, which is going to evolve 
  % according to the new LXCat3

  % reaction type and subtype
  % read only the first typeTag
  LXCatEntryArray(end+1).type = description.typeTags{1,1};
  % if typeTag is Electronic, change to Excitation
  % this is going to evolve according to the new LXCat3 - WIP
  if strcmpi(description.typeTags{1,1}, 'Electronic')
    LXCatEntryArray(end).type = 'Excitation';
  end


  % check if it's momentum transfer cross-section
  % This is a WIP, which is going to evolve according to the new LXCat3
  LXCatEntryArray(end).isMomentumTransfer = any(strcmp(description.typeTags, ...
    'MomentumTransfer'));  

  % reversible reaction
  LXCatEntryArray(end).isReverse = description.reversible == 1;

  % add the raw data to the structure
  LXCatEntryArray(end).rawCrossSection = rawCrossSection;

  % get the fieldnames of the json states
  statesFieldnamesStr = fieldnames(jsonStates);
  statesFieldnames = eraseBetween(statesFieldnamesStr,1,1) ;

  % initialize reactantElectrons and productElectrons count
  LXCatEntryArray(end).reactantElectrons = 0;
  LXCatEntryArray(end).productElectrons = 0;

  LXCatEntryArray(end).target = struct.empty;


  % check if the lhs is empty
  if isempty(description.lhs)
      error('A target can not be found in the collision. Please check your LXCat files');
  end
  % loop for each species in the targets (lhs)
  for lhsId = 1:size(description.lhs,1)

      % find the state in the states structure
      lhsState = description.lhs(lhsId).state;
      stateId = find(strcmp(statesFieldnames, lhsState), 1);
      
      if ~isempty(stateId)
          % get the state corresponding to the stateId in the left-hand side
          stateField = jsonStates.(statesFieldnamesStr{stateId});
          if ~isempty(stateField.detailed)
              if isfield(stateField.detailed, 'type') && strcmpi(stateField.detailed.type, 'Electron')
              % if the state is an electron add its count to the reactantElectrons
              % field
                  reactantElectrons = description.lhs(lhsId).count;
                  LXCatEntryArray(end).reactantElectrons = reactantElectrons;
              else
                
                % initialize target fields
                LXCatEntryArray(end).target(end+1).ionCharg = '';
                LXCatEntryArray(end).target(end).gasName = '';
                LXCatEntryArray(end).target(end).quantity = '';
                LXCatEntryArray(end).target(end).eleLevel = '';
                LXCatEntryArray(end).target(end).vibLevel = '';
                LXCatEntryArray(end).target(end).rotLevel = '';

                % if the species is not an electron, add it to the the target structure
                if isfield(stateField.detailed, 'charge') && stateField.detailed.charge ~= 0
                    % add the ion charge to the target structure
                    if stateField.detailed.charge == 1
                        LXCatEntryArray(end).target(end).ionCharg = '+';
                    elseif stateField.detailed.charge == -1
                        LXCatEntryArray(end).target(end).ionCharg = '-';
                    else
                        error('Invalid ion charge.');
                    end
                end
                if isfield(stateField.detailed, 'composition')
                  gasName = '';
                  for cNum = 1:length(stateField.detailed.composition)
                    if isequal(stateField.detailed.composition{cNum}{2},1)
                        gasName = strcat(gasName, stateField.detailed.composition{cNum}{1});
                    else
                        gasName = strcat(gasName, num2str(stateField.detailed.composition{cNum}{1}), ...
                            num2str(stateField.detailed.composition{cNum}{2}));
                    end
                  end
                  LXCatEntryArray(end).target(end).gasName = gasName;
                      
                end
                % fill the target structure with the electronic, vibrational and rotational levels
                [LXCatEntryArray(end).target(end).eleLevel, LXCatEntryArray(end).target(end).vibLevel, ...
                    LXCatEntryArray(end).target(end).rotLevel] = parseStateLevelInfoJson(stateField);

                % Assign quantity to target
                LXCatEntryArray(end).target.quantity = description.lhs(lhsId).count;
    
              end
          end

      end

  end

  %% Products
  LXCatEntryArray(end).productArray = struct.empty;

  % check if the rhs is empty
  if isempty(description.rhs)
    error('A product can not be found in the collision. Please check your LXCat files');
  end

  % loop for each species in the products (rhs)
  for rhsId = 1:size(description.rhs,1)

    % find the state in the states structure
    rhsState = description.rhs(rhsId).state;
    stateId = find(strcmp(statesFieldnames, rhsState), 1);
    
    if ~isempty(stateId)
        stateField = jsonStates.(statesFieldnamesStr{stateId});
        if ~isempty(stateField.detailed)
            if isfield(stateField.detailed, 'type') && strcmpi(stateField.detailed.type, 'Electron')
            % if the state is an electron add its count to the productElectrons
            % field
                productElectrons = description.rhs(rhsId).count;
                LXCatEntryArray(end).productElectrons = productElectrons;
            else
              % if the species is not an electron, add it to the the
              % productArray structure
              
              % initialize product fields
              LXCatEntryArray(end).productArray(end+1).ionCharg = '';
              LXCatEntryArray(end).productArray(end).gasName = '';
              LXCatEntryArray(end).productArray(end).quantity = '';
              LXCatEntryArray(end).productArray(end).eleLevel = '';
              LXCatEntryArray(end).productArray(end).vibLevel = '';
              LXCatEntryArray(end).productArray(end).rotLevel = '';


              if isfield(stateField.detailed, 'charge') && stateField.detailed.charge ~= 0
                if stateField.detailed.charge == 1
                    LXCatEntryArray(end).productArray(end).ionCharg = '+';
                elseif stateField.detailed.charge == -1
                    LXCatEntryArray(end).productArray(end).ionCharg = '-';
                else
                    error('Invalid ion charge.');
                end
              end
              if isfield(stateField.detailed, 'composition')
                  gasName = '';
                  for cNum = 1:length(stateField.detailed.composition)
                      if isequal(stateField.detailed.composition{cNum}{2},1)
                          gasName = strcat(gasName, stateField.detailed.composition{cNum}{1});
                      else
                          gasName = strcat(gasName, num2str(stateField.detailed.composition{cNum}{1}), ...
                              num2str(stateField.detailed.composition{cNum}{2}));
                      end
                  end
                  LXCatEntryArray(end).productArray(end).gasName = gasName;
              end
              % fill the productArray structure
              [LXCatEntryArray(end).productArray(end).eleLevel, LXCatEntryArray(end).productArray(end).vibLevel, ...
                  LXCatEntryArray(end).productArray(end).rotLevel] = parseStateLevelInfoJson(stateField);

              % Assign quantity to product Array
              LXCatEntryArray(end).productArray(end).quantity = description.rhs(rhsId).count;
            end
        end
    end
  end

  % 'parameter' represents the threshold
  if isempty(parameter)
    LXCatEntryArray(end).threshold = 0;
  else
    threshold = parameter;
    if isnan(threshold) 
      error('I can not properly parse ''%s'' as the threshold of reaction ''%s''.\nPlease check your LXCat files', ...
        parameter.threshold, [description.reactants description.direction description.products]);
    elseif threshold<0
      error('Reaction ''%s'' has a negative threshold (''%s'').\nPlease check your LXCat files', ...
        [description.reactants description.direction description.products], parameter.threshold);
    end

    LXCatEntryArray(end).threshold = threshold;
  end

end



function [eleLevel, vibLevel, rotLevel] = parseStateLevelInfoJson(stateField)
% parseStateLevelInfoJson Parses the electronic, vibrational, and rotational levels from the state field's
% serialized information. Reads the summary field of the electronic, vibrational, and rotational levels.
%
% Inputs:
%   stateField - State field from LXCat data in JSON format
%
% Outputs:
%   eleLevel - Electronic level string
%   vibLevel - Vibrational level string  
%   rotLevel - Rotational level string


  eleLevel = '';
  vibLevel = '';
  rotLevel = '';

  % electronic level
  if ~isempty(stateField.serialized) && isfield(stateField.serialized,'electronic')
    if ~isempty(stateField.serialized.electronic) && isfield(stateField.serialized.electronic, 'summary')
      % check if size is 1x1, else join with |
      if size(stateField.serialized.electronic, 1) == 1 && size(stateField.serialized.electronic, 2) == 1
        eleLevel = stateField.serialized.electronic.summary;
      else
        summaries = {stateField.serialized.electronic.summary};
        eleLevel = strjoin(summaries, '|');
      end

      % vibrational level
      if ~isempty(stateField.serialized.electronic) && isfield(stateField.serialized.electronic, 'vibrational')
        if ~isempty(stateField.serialized.electronic.vibrational) && isfield(stateField.serialized.electronic.vibrational, 'summary')
          % check if size is 1x1, else join with |
          if size(stateField.serialized.electronic.vibrational, 1) == 1 && size(stateField.serialized.electronic.vibrational, 2) == 1
            vibLevel = stateField.serialized.electronic.vibrational.summary;
          else
            summaries = {stateField.serialized.electronic.vibrational.summary};
            vibLevel = strjoin(summaries, '|');
          end
          
          % rotational level
          if ~isempty(stateField.serialized.electronic.vibrational) && isfield(stateField.serialized.electronic.vibrational, 'rotational')
            if ~isempty(stateField.serialized.electronic.vibrational.rotational) && ...
                isfield(stateField.serialized.electronic.vibrational.rotational, 'summary')
              % check if size is 1x1, else join with |
              if size(stateField.serialized.electronic.vibrational.rotational, 1) == 1 && size(stateField.serialized.electronic.vibrational.rotational, 2) == 1
                rotLevel = stateField.serialized.electronic.vibrational.rotational.summary;
              else
                summaries = {stateField.serialized.electronic.vibrational.rotational.summary};
                rotLevel = strjoin(summaries, '|');
              end
            else
              error(['There is no summary available for rotational level of state %s . Please check your LXCat files.'], ...
              [stateField.serialized.summary])
            end
          end

        else
          error(['There is no summary available for vibrational level of state %s . Please check your LXCat files.'], ...
          [stateField.serialized.summary])
        end
      end
    else
      error(['There is no summary available for electronic level of state %s . Please check your LXCat files.'], ...
      [stateField.serialized.summary])
    end
    
  end

end


function LXCatEntryArray = addLXCatEntry(description, parameter, rawCrossSection, LXCatEntryArray)
% addLXCatEntry analyses the information of a particular LXCat entry and adds it to the structure array
% LXCatEntryArray.
      
      
  % if type is Electronic, change to Excitation
  % This is a WIP, according to the development of the new LXCat3.0
  if strcmpi(description.type,'Electronic')
    LXCatEntryArray(end+1).type = 'Excitation';
  else
    LXCatEntryArray(end+1).type = description.type;
  end

  if isempty(description.subtype)
    LXCatEntryArray(end).isMomentumTransfer = false;
  else
    LXCatEntryArray(end).isMomentumTransfer = true;
  end
  if strcmp(description.direction, '->')
    LXCatEntryArray(end).isReverse = false;
  elseif strcmp(description.direction, '<->')
    LXCatEntryArray(end).isReverse = true;
  end
  LXCatEntryArray(end).target = regexp(description.reactants, Parse.stateRegExp, 'names', 'once');
  if isempty(LXCatEntryArray(end).target)
    error(['I can not find a target in the collision:\n%s\nthat match the regular expression for a state.\n'...
      'Please check your LXCat files'], [description.reactants description.direction description.products]);
  end
  electronArray = regexp(description.reactants, Parse.electronRegExp, 'names');
  numElectrons = 0;
  for i = 1:length(electronArray)
    if isempty(electronArray(i).quantity)
      numElectrons = numElectrons+1;
    else
      numElectrons = numElectrons+str2double(electronArray(i).quantity);
    end
  end
  LXCatEntryArray(end).reactantElectrons = numElectrons;
  productArray = regexp(description.products, Parse.stateRegExp, 'names');
  LXCatEntryArray(end).productArray = removeDuplicatedStates(productArray);
  if isempty(LXCatEntryArray(end).productArray)
    error(['I can not find a product in the collision:\n%s\nthat match the regular expression for a state.\n'...
      'Please check your LXCat files'], [description.reactants description.direction description.products]);
  end
  electronArray = regexp(description.products, Parse.electronRegExp, 'names');
  numElectrons = 0;
  for i = 1:length(electronArray)
    if isempty(electronArray(i).quantity)
      numElectrons = numElectrons+1;
    else
      numElectrons = numElectrons+str2double(electronArray(i).quantity);
    end
  end
  LXCatEntryArray(end).productElectrons = numElectrons;
  if isempty(parameter)
    LXCatEntryArray(end).threshold = 0;
  else
    threshold = str2double(parameter.threshold);
    if isnan(threshold) 
      error('I can not properly parse ''%s'' as the threshold of reaction ''%s''.\nPlease check your LXCat files', ...
        parameter.threshold, [description.reactants description.direction description.products]);
    elseif threshold<0
      error('Reaction ''%s'' has a negative threshold (''%s'').\nPlease check your LXCat files', ...
        [description.reactants description.direction description.products], parameter.threshold);
    end
    LXCatEntryArray(end).threshold = threshold;
  end
  LXCatEntryArray(end).rawCrossSection = rawCrossSection;
  
end

function newStateArray = removeDuplicatedStates(stateArray)
  
  if isempty(stateArray)
    newStateArray = struct.empty;
    return
  end
  
  numStates = length(stateArray);
  newStateArray = stateArray(1);
  if isempty(newStateArray.quantity)
    newStateArray.quantity = 1;
  else
    newStateArray.quantity = str2double(newStateArray.quantity);
  end
  for i = 2:numStates
    for j = 1:length(newStateArray)
      if strcmp(stateArray(i).gasName, newStateArray(j).gasName) && ...
          strcmp(stateArray(i).ionCharg, newStateArray(j).ionCharg) && ...
          strcmp(stateArray(i).eleLevel, newStateArray(j).eleLevel) && ...
          strcmp(stateArray(i).vibLevel, newStateArray(j).vibLevel) && ...
          strcmp(stateArray(i).rotLevel, newStateArray(j).rotLevel)
        if isempty(stateArray(i).quantity)
          newStateArray(j).quantity = newStateArray(j).quantity+1;
        else
          newStateArray(j).quantity = newStateArray(j).quantity+str2double(stateArray(i).quantity);
        end
        break;
      end
      if j == length(newStateArray)
        newStateArray(end+1) = stateArray(i);
        if isempty(stateArray(i).quantity)
          newStateArray(end).quantity = 1;
        else
          newStateArray(end).quantity = str2double(stateArray(i).quantity);
        end
      end
    end
  end
  
end
