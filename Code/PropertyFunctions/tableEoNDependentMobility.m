function coefficient = tableEoNDependentMobility(state, argumentArray, workCond)
  % tableEoNDependentMobility is a property function that interpolates 
  %  the reduced mobility for the current value of the reduced electric field
  %  from a table with filename specified as the first argument provided to the function in the setup file

  persistent tableArrayState;

  % check if this is the first evaluation for this state
  if isempty(tableArrayState)
    tableArrayState = cell.empty;
    firstEvaluation = true;
  elseif length(tableArrayState) < state.ID
    firstEvaluation = true;
  else
    firstEvaluation = isempty(tableArrayState{state.ID});
  end

  % read the table (only first time) and save it for further evaluations
  if firstEvaluation
    % open file with the table
    fileName = argumentArray{1};
    fileName = ['Input' filesep fileName];
    fileID = fopen(fileName,'r');
    if fileID < 0
      error(['Error found when evaluating mobility of state %s, while using ''tableEfieldDependentMobility''\n'...
        'File ''%s'' could not be opened.'],state.name,fileName);
    end
    % parse file
    tableArrayState{state.ID} = [];
    while ~feof(fileID)
      % get line and eliminate the comments
      line = fgetl(fileID);
      idx = regexp(line, '%', 'once'); 
      if ~isempty(idx)
        line = line(1:idx-1);
      end
      if isempty(line) || all(isspace(line))
        continue;
      end
      tableArrayState{state.ID} = [tableArrayState{state.ID}; str2num(line)];
    end
    fclose(fileID);
  end  

  % retrieve value of reduced electric field
  EoN = workCond.reducedField;

  % interpolate the coefficient from the table
  if EoN <= tableArrayState{state.ID}(1,1)
    coefficient = tableArrayState{state.ID}(1,2);
  elseif EoN >= tableArrayState{state.ID}(end,1)
    coefficient = tableArrayState{state.ID}(end,2);
  else
    coefficient = interp1(tableArrayState{state.ID}(:,1), tableArrayState{state.ID}(:,2), EoN, 'linear');
  end  
  
end
