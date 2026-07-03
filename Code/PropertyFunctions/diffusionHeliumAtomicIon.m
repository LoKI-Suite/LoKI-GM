function redDiffCoeff = diffusionHeliumAtomicIon(state, argumentArray, workCond)
   % diffusionHeliumAtomicIon is a property function for the reduced diffusion coefficient of Helium atomic ions 
   %  as described by M.Santos et al, J. Phys. D: Appl. Phys, 47, 265201, 2014 
  
  A = argumentArray{1};
  B = argumentArray{2};
  C = argumentArray{3};
  temperature = argumentArray{4};
  if ~isnumeric(temperature)
    switch temperature
      case 'gasTemperature'
        temperature = workCond.gasTemperature;
      case 'electronTemperature'
        temperature = workCond.electronTemperature;
      otherwise
        error(['Error found when evaluating reduced diffusion coefficient of state %s.\nTemperature ''%s'' not ' ...
          'defined in the working conditions.\nPlease, fix the problem and run the code again.'], ...
          state.name, temperature);
    end
  end
  
  redMob = A/(B*sqrt(temperature)+C);
  redDiffCoeff = redMob*temperature*Constant.boltzmann/Constant.electronCharge;
  
end
