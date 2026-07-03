function population = boltzmannRotationalPopulation_H2(state, argumentArray, workCond)
  % boltzmannRotationalPopulation_H2 is a property function that evaluates the population of rotational state j of H2, 
  %  adopting Boltzmann distributions with a certain temperature 
  %  (specified as the first argument provided to the function in the setup file) 
  %  for the ortho/para configurations, with weights 0.25 and 0.75, respectively
  
  persistent normEven
  persistent normOdd
  persistent iterNorm
  
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
  
  % initialize separate normalization for ortho/para H2 rotational states
  rotdim = length(state.siblingArray);
  if isempty(normEven) && isempty(normOdd)
    normEven = 0;
    normOdd = 0;
    iterNorm = 1;
  else
    iterNorm = iterNorm + 1;  
  end

  % evaluate Boltzmann distribution for H2 rotational states, separating the populations of ortho/para configurations 
  for stateAux = state %[state state.siblingArray]
    if isempty(stateAux.energy)
      error(['Unable to find %s energy for the evaluation of ''boltzmannPopulation'' function.\n'...
        'Check input file'], stateAux.name);
    elseif isempty(stateAux.statisticalWeight)
      error(['Unable to find %s statistical weight for the evaluation of ''boltzmannPopulation'' '...
        'function.\nCheck input file'], stateAux.name);
    end

    if ~strcmp(state.type, 'rot')
      error(['Trying to asign Boltzmann populations to non rotational state %s. Check input file', state.name]);
    end
    J = str2double(stateAux.rotLevel);

    if(rem(J,2) == 0)
      stateAux.population = stateAux.statisticalWeight*exp(-stateAux.energy/(Constant.boltzmannInEV*temperature));
      normEven = normEven + stateAux.population;
    else
      stateAux.population = stateAux.statisticalWeight*exp(-stateAux.energy/(Constant.boltzmannInEV*temperature));
      normOdd = normOdd + stateAux.population;
    end
  end
  
  if (iterNorm == rotdim+1)
    for stateAux = [state state.siblingArray]
        J = str2double(stateAux.rotLevel);
        if(rem(J,2) == 0)
            stateAux.population = 0.25*stateAux.population/normEven;
        else
            stateAux.population = 0.75*stateAux.population/normOdd;
        end
    end
  end
  
  % return population of the current state
  population = state.population;
  
end
