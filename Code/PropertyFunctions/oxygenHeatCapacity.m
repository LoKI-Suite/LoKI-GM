function heatCapacity = oxygenHeatCapacity(gas, ~, workCond)
  % oxygenHeatCapacity is a property function that evaluates the molecular oxygen heat capacity
  % taken from B J McBride, M J Zehe and S Gordon 2002 NASA/TP-2002-211556
  %  "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"

  Tg = workCond.gasTemperature;
  
  if Tg >= 200 && Tg <= 1000
    a1 = -3.42556342E+04;
    a2 = 4.84700097E+02;
    a3 = 1.11901096E+00;
    a4 = 4.29388924E-03;
    a5 = -6.83630052E-07;
    a6 = -2.02337270E-09;
    a7 = 1.03904002E-12;
    b1 = -3.39145487E+03;
    b2 = 1.84969947E+01;
  elseif Tg > 1000 && Tg <= 6000
    a1 = -1.03793902E+06;
    a2 = 2.34483028E+03;
    a3 = 1.81973204E+00;
    a4 = 1.26784758E-03;
    a5 = -2.18806799E-07;
    a6 = 2.05371957E-11;
    a7 = -8.19346705E-16;
    b1 = -1.68901093E+04;
    b2 = 1.73871651E+01;
  else 
    a1 = 4.97529430E+08;
    a2 = -2.86610687E+05;
    a3 = 6.69035225E+01;
    a4 = -6.16995902E-03;
    a5 = 3.01639603E-07;
    a6 = -7.42141660E-12;
    a7 = 7.27817577E-17;
    b1 = 2.29355403E+06;
    b2 = -5.53062161E+02;
  end

  % heat capacity for O2 (SI units)
  heatCapacity = 8.314510*(a1*Tg^(-2)+a2*Tg^(-1)+a3+a4*Tg+a5*Tg^2+a6*Tg^3+a7*Tg^4);

  % change energy units to eV
  heatCapacity = heatCapacity/Constant.electronCharge;
  % store value on gas object properties
  gas.heatCapacity = heatCapacity;
    
end
