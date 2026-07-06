function thermalConductivity = oxygenAtomThermalConductivity(gas, ~, workCond)
  % oxygenAtomThermalConductivity is a property function that evaluates the atomic oxygen thermal conductivity
  % taken from B J McBride, S Gordon and M A Reno 1993 NASA TM-4513
  %  "Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species"
  
  Tg = workCond.gasTemperature;

  % thermal conductivity for O (SI units)
  if Tg >= 300 && Tg <= 1000
   thermalConductivity = 1e-4*exp(0.73824503*log(Tg) + 0.11221345e2/Tg + 0.31668244e4/((Tg)^2) + 0.17085307e1);
  else
   thermalConductivity = 1e-4*exp(0.79819261*log(Tg) + 0.17970493e3/Tg - 0.52900889e5/((Tg)^2) + 0.11797640e1);
  end

  % change energy units to eV
  thermalConductivity = thermalConductivity/Constant.electronCharge;
  % store value on gas object properties
  gas.thermalConductivity = thermalConductivity;
    
end
