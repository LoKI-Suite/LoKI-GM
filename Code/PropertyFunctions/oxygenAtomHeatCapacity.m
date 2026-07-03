function heatCapacity = oxygenAtomHeatCapacity(gas, ~, workCond)
  % oxygenAtomHeatCapacity is a property function that evaluates the atomic oxygen heat capacity
  % taken from B J McBride, M J Zehe and S Gordon 2002 NASA/TP-2002-211556
  %  "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"

  Tg = workCond.gasTemperature;
  
  if Tg >= 200 && Tg <= 1000
    a1 = -7.95361130E+03;
    a2 = 1.60717779E+02;
    a3 = 1.96622644E+00;
    a4 = 1.01367031E-03;
    a5 = -1.11041542E-06;
    a6 = 6.51750750E-10;
    a7 = -1.58477925E-13;
    b1 = 2.84036244E+04;
    b2 = 8.40424182E+00;
  elseif Tg > 1000 && Tg <= 6000
    a1 = 2.61902026E+05;
    a2 = -7.29872203E+02;
    a3 = 3.31717727E+00;
    a4 = -4.28133436E-04;
    a5 = 1.03610459E-07;
    a6 = -9.43830433E-12;
    a7 = 2.72503830E-16;
    b1 = 3.39242806E+04;
    b2 = -6.67958535E-01;
  else 
    a1 = 1.77900426E+08;
    a2 = -1.08232826E+05;
    a3 = 2.81077837E+01;
    a4 = -2.97523226E-03;
    a5 = 1.85499753E-07;
    a6 = -5.79623154E-12;
    a7 = 7.19172016E-17;
    b1 = 8.89094263E+05;
    b2 = -2.18172815E+02;
  end

  % heat capacity for O (SI units)
  heatCapacity = 8.314510*(a1*Tg^(-2)+a2*Tg^(-1)+a3+a4*Tg+a5*Tg^2+a6*Tg^3+a7*Tg^4);

  % change energy units to eV
  heatCapacity = heatCapacity/Constant.electronCharge;
  % store value on gas object properties
  gas.heatCapacity = heatCapacity;
    
end
