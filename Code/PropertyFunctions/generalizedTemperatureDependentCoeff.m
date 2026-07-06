function coefficient = generalizedTemperatureDependentCoeff(state, argumentArray, workCond)
  % generalizedTemperatureDependentCoeff is a property function that evaluates a generic 
  %  reduced mobility or reduced diffusion coefficient with the form aT^b, 
  %  where the normalization constant a, the temperature T and the power c are specified as 
  %  the first three arguments provided to the function in the setup file
  
  normalizationConstant = argumentArray{1};
  temperature = argumentArray{2};
  power = argumentArray{3};
  if ~isnumeric(temperature)
    switch temperature
      case 'gasTemperature'
        temperature = workCond.gasTemperature;
      case 'electronTemperature'
        temperature = workCond.electronTemperature;
      otherwise
        error(['Error found when evaluating coefficient of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temperature);
    end
  end
  
  coefficient = normalizationConstant*temperature^power;
  
end
