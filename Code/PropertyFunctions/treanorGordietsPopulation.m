function population = treanorGordietsPopulation(state, argumentArray, workCond)
  % treanorGordietsPopulation is a property function that evaluates a Treanor-Gordiets distribution
  %  with certain temperatures (specified as the first two arguments provided to the function in the setup file)
  
  if ~strcmp(state.type, 'vib')
    error(['Trying to asign treanor gordiets population to non vibrational state '...
      '%s. Check input file', stateAux.name]);
  end
  
  temp0 = argumentArray{1};
  temp1 = argumentArray{2};
  if ~isnumeric(temp0)
    switch temp0
      case 'gasTemperature'
        temp0 = workCond.gasTemperature;
      case 'electronTemperature'
        temp0 = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temp0);
    end
  end
  if ~isnumeric(temp1)
    switch temp1
      case 'gasTemperature'
        temp1 = workCond.gasTemperature;
      case 'electronTemperature'
        temp1 = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temp1);
    end
  end
  
  anharmonicFrequency = state.gas.anharmonicFrequency;
  if isempty(anharmonicFrequency)
    error(['Unable to find %s harmonic frequency to populate state %s and its siblings with function %s.\n' ...
      'Check input file'], state.gas.name, state.name, 'treanorGordietsPopulation');
  end
  E0 = [];
  E1 = [];
  for stateAux = [state state.siblingArray]
    if strcmp(stateAux.vibLevel, '0')
      E0 = stateAux.energy;
    elseif strcmp(stateAux.vibLevel, '1')
      E1 = stateAux.energy;
    end
  end
  if isempty(E0) || isempty(E1)
    error('Unable to find E0 or E1 to populate state %s and its siblings with function %s.\nCheck input file', ...
      state.name, 'treanorGordietsPopulation');
  end
  vLimit = floor(0.5*(1+(E1-E0)*temp0/(Constant.planckReducedInEV*anharmonicFrequency*temp1)));
  
  % evaluate Treanor-Gordiets (Treanor modified above vLimit)
  vLimitPop = [];
  norm = 0;
  for stateAux = [state state.siblingArray]
    v = str2double(stateAux.vibLevel);
    if v <= vLimit
      E = stateAux.energy;
      g = stateAux.statisticalWeight;
      if isempty(E)
        error('Unable to find %s energy for the evaluation of %s function.\nCheck input file', ...
          stateAux.name, 'treanorGordietsPopulation');
      elseif isempty(g)
        error('Unable to find %s statistical weight for the evaluation of %s function.\nCheck input file', ...
          stateAux.name, 'treanorGordietsPopulation');
      end
      if v < vLimit
        stateAux.population = g*exp(-(v*(E1-E0)*(1/temp1-1/temp0)+(E-E0)/temp0)/Constant.boltzmannInEV);
      else
        vLimitPop = g*exp(-(v*(E1-E0)*(1/temp1-1/temp0)+(E-E0)/temp0)/Constant.boltzmannInEV);
        stateAux.population = vLimitPop;
      end
      norm = norm + stateAux.population;
    end
  end
  for stateAux = [state state.siblingArray]
    v = str2double(stateAux.vibLevel);
    if v > vLimit
      stateAux.population = vLimitPop*vLimit/v;
      norm = norm + stateAux.population;
    end
  end
  for stateAux = [state state.siblingArray]
    stateAux.population = stateAux.population/norm;
  end
  
  population = state.population;
  
end
