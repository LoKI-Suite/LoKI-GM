function thermalConductivity = nitrogenThermalConductivity(gas, ~, workCond)
  % nitrogenThermalConductivity is a property function that evaluates the nitrogen thermal conductivity
  % taken from C D Pintassilgo et al 2014 Plasma Sources Sci. Technol. 23 025006
  % https://iopscience.iop.org/article/10.1088/0963-0252/23/2/025006/pdf

  % thermal conductivity for nitrogen (SI units)
  thermalConductivity = (1.717+0.084*workCond.gasTemperature-1.948e-5*workCond.gasTemperature^2)*1e-3;
  % change energy units to eV
  thermalConductivity = thermalConductivity/Constant.electronCharge;
  % store value on gas object properties
  gas.thermalConductivity = thermalConductivity;
  
  
end
