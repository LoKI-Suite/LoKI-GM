function population = boltzmannPopulationVibrationalCutoff(state, argumentArray, workCond)
  % boltzmannPopulationVibrationalCutoff is a property function that evaluates a Boltzmann distribution with a certain
  % temperature (specified as the first argument provided to the function in the setup file) and up to a certain 
  % vibrational level (specified as the second argument provided to the function in the setup file). This function
  % should be used ONLY to provide the population of vibrational distributions functions, for other distributions check
  % boltzmannPopulation and boltzmannPopulationRotationalLevel.
  
  % obtain temperature of the distribution (either prescribed, i.e. numeric, or found in the working conditions)
  temperature = argumentArray{1};
  if ~isnumeric(temperature)
    switch temperature
      case 'gasTemperature'
        temperature = workCond.gasTemperature;
      case 'electronTemperature'
        temperature = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temperature);
    end
  end
  
  % obtain the cutoff vibrational level for the Boltzmann distribution (this level is included in the distribution)
  vMax = argumentArray{2};
  if ~isnumeric(vMax)
    error(['Error found when evaluating population of state %s.\nMaximum vibrational level ''%s'' should be a ' ...
      'numerical value.\nPlease, fix the problem and run the code again.'], state.name, argumentArray{2});
  end
  
  % error checking (find if state is a vibrational level)
  if ~strcmp(state.type, 'vib')
    error(['Trying to asign ''boltzmannPopulationVibrationalCutoff'' population to non vibrational state %s. ' ...
      'Check input file'], state.name);
  end
  
  % evaluate Boltzmann distribution for the state and its siblings
  norm = 0;
  for stateAux = [state state.siblingArray]
    if isempty(stateAux.energy)
      error(['Unable to find %s energy for the evaluation of ''boltzmannPopulationVibrationalCutoff'' function.\n'...
        'Check input file'], stateAux.name);
    elseif isempty(stateAux.statisticalWeight)
      error(['Unable to find %s statistical weight for the evaluation of ''boltzmannPopulationVibrationalCutoff'' '...
        'function.\nCheck input file'], stateAux.name);
    end
    if str2double(stateAux.vibLevel) > vMax
      stateAux.population = 0;
    else
      stateAux.population = stateAux.statisticalWeight*exp(-stateAux.energy/(Constant.boltzmannInEV*temperature));
      norm = norm + stateAux.population;
    end
  end
  for stateAux = [state state.siblingArray]
    stateAux.population = stateAux.population/norm;
  end
  
  % return population of the current state
  population = state.population;
  
end
