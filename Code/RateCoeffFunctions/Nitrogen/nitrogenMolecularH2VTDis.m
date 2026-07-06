function [rateCoeff, dependent] = nitrogenMolecularH2VTDis(~, ~, ~, reaction, rateCoeffParams, chemistry)
% nitrogenMolecularH2VTdis evaluates the rate coefficient for N2 dissociation from V-T transitions between nitrogen molecules and H2 described in:
% L.L. Alves, L. Marques, Plasma Sources Sci. Technol.  (2020)
% http://dx.doi.org/10.1088/

  % local copies and definitions of different parameters used in the function
  kb = Constant.boltzmann;
  Tg = chemistry.workCond.gasTemperature;
  uma = Constant.unifiedAtomicMass; 
  hbar = Constant.planckReduced;
  Lmin = 3.5e-11;                                                                   % Minimum distance between N2 molecules during V-V collision (m)
  MN2 = reaction.reactantArray(1).gas.mass;  % Mass of N2 (Kg)
  MH2 = reaction.catalystArray(1).gas.mass;  % Mass of H2 (Kg)

  omega = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
  chi = reaction.reactantArray(1).gas.anharmonicFrequency/omega; % Anharmonicity of the oscillator N2
  NLevels = rateCoeffParams{1};
  
  % identify current V-T reaction by evaluating the vibrational levels ( v + X <-> v-1 + X )
  v = str2double(reaction.reactantArray(1).vibLevel);
  
  if (NLevels-v)~=1  % check the reaction
    error(['nitrogenMolecularH2VT can not evaluate the rate coefficient for the reaction:\n%s\n' ...
      'The reaction does not comply with the following structure:\n v + X <-> v-1 + X;\n'], ...
      reaction.description);
  end
  % intermediate calculations
  uM = 1.0/((1/MH2)+(1/MN2));  %Reduced mass
  Y = (1-2*chi)*(0.5^1.5)*2*pi*omega*Lmin*(uM/(kb*Tg))^0.5;
  if(Y>=0 & Y<=20)
      FY1=0.5*(3-exp(-2*Y/3))*exp(-2*Y/3);
  elseif (Y>20)
      FY1=8*(pi/3)^0.5*Y^(7/3)*exp(-3*Y^(2/3));
  else  % check omega,chi in VT N2 with H2.
      error(['nitrogenMolecularVT with H2 can not evaluate the rate coefficient for the reaction:\n%s\n'], ...
        reaction.description);
  end
  P10 = 5.63d-20*Tg*exp(-80.6*Tg^(-1/3));  %m3s-1
      
  Y = (1-2*chi*NLevels)*(0.5^1.5)*2*pi*omega*Lmin*(uM/(kb*Tg))^0.5;
  if(Y>=0 & Y<=20)
      FY=0.5*(3-exp(-2*Y/3))*exp(-2*Y/3);
  elseif (Y>20)
      FY=8*(pi/3)^0.5*Y^(7/3)*exp(-3*Y^(2/3));
  else  % check VT N2 with H2.
      error(['nitrogenMolecularVT with H2 can not evaluate the rate coefficient for the reaction:\n%s\n'], ...
        reaction.description);
  end

  % evaluate rate coefficient
  rateCoeff = NLevels*(1-chi)/(1-chi*NLevels)*P10*FY/FY1;
  rateCoeff = rateCoeff*exp(-hbar*omega*(1-2*chi*NLevels)/(kb*Tg));
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  
end
