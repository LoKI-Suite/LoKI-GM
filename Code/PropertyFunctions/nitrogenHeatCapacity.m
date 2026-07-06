function heatCapacity = nitrogenHeatCapacity(gas, ~, workCond)
  % nitrogenHeatCapacity is a property function that evaluates the nitrogen heat capacity
  % taken from C D Pintassilgo et al 2014 Plasma Sources Sci. Technol. 23 025006
  % https://iopscience.iop.org/article/10.1088/0963-0252/23/2/025006/pdf
  
  % heat capacity for nitrogen (SI units)
  heatCapacity = 29.1 + 2494.2/(553.4*sqrt(pi/2))*exp(-2*((workCond.gasTemperature-1047.4)/553.4)^2);
  % change energy units to eV
  heatCapacity = heatCapacity/Constant.electronCharge;
  % store value on gas object properties
  gas.heatCapacity = heatCapacity;
  
  
end
