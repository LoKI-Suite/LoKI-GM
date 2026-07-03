function thermalConductivity = oxygenThermalConductivity(gas, ~, workCond)
  % oxygenThermalConductivity is a property function that evaluates the molecular oxygen thermal conductivity
  % taken from B J McBride, S Gordon and M A Reno 1993 NASA TM-4513
  %  "Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species"
  
  Tg = workCond.gasTemperature;

  % thermal conductivity for O2 (SI units)
  if Tg >= 300 && Tg <= 1000
   thermalConductivity = 1e-4*exp(0.81595343*log(Tg) - 0.34366856e2/Tg + 0.22785080e4/((Tg)^2) + 0.10050999E1);
  else
   thermalConductivity = 1e-4*exp(0.80805788*log(Tg) + 0.11982181e3/Tg - 0.47335931e5/((Tg)^2) + 0.95189193);
  end

  % change energy units to eV
  thermalConductivity = thermalConductivity/Constant.electronCharge;
  % store value on gas object properties
  gas.thermalConductivity = thermalConductivity;
    
end
