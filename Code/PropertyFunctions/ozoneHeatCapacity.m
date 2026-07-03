function heatCapacity = ozoneHeatCapacity(gas, ~, workCond)
  % ozoneHeatCapacity is a property function that evaluates the ozone heat capacity
  % taken from B J McBride, M J Zehe and S Gordon 2002 NASA/TP-2002-211556
  %  "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
  % The effective range for ozone heat capacity is 200K to 6000K

  Tg = workCond.gasTemperature;
  
  if Tg >= 200 && Tg <= 1000
    a1 = -1.28231451E+04;
    a2 = 5.89821664E+02;
    a3 = -2.54749676E+00;
    a4 = 2.69012153E-02;
    a5 = -3.52825834E-05;
    a6 = 2.31229092E-08;
    a7 = -6.04489327E-12;
    b1 = 1.34836870E+04;
    b2 = 3.85221858E+01;
  else 
    a1 = -3.86966248E+07;
    a2 = 1.02334499E+05;
    a3 = -8.96155160E+01;
    a4 = 3.70614497E-02;
    a5 = -4.13763874E-06;
    a6 = -2.72501859E-10;
    a7 = 5.24818811E-14;
    b1 = -6.51791818E+05;
    b2 = 7.02910952E+02;
  end

  % heat capacity for O3 (SI units)
  heatCapacity = 8.314510*(a1*Tg^(-2)+a2*Tg^(-1)+a3+a4*Tg+a5*Tg^2+a6*Tg^3+a7*Tg^4);

  % change energy units to eV
  heatCapacity = heatCapacity/Constant.electronCharge;
  % store value on gas object properties
  gas.heatCapacity = heatCapacity;
    
end
