function [rateCoeff, dependent] = nitrogenAtomicHVTDis(~, ~, ~, reaction, rateCoeffParams, chemistry)
% nitrogenAtomicHVTdis evaluates the rate coefficient for N2 dissociation from V-T transitions between nitrogen molecules and H atoms
% described in:
% L.L. Alves, L. Marques, Plasma Sources Sci. Technol. (2020) 
% http://dx.doi.org/10.1088/

  % local copies and definitions of different parameters used in the function
  kb = Constant.boltzmann;
  hbar = Constant.planckReduced;
  Tg = chemistry.workCond.gasTemperature;
  omegaN = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
  chiN = reaction.reactantArray(1).gas.anharmonicFrequency/omegaN; % Anharmonicity of the oscillator N2
  E10N = hbar*omegaN*(1-2*chiN);
  deltaEN = hbar*omegaN*chiN;
  
  % identify current V-T reaction by evaluating the vibrational levels ( N2(X,v) + H <-> N +N + H, v-w<6 ) 
  w = str2double(reaction.reactantArray(1).vibLevel);
  NLevels = rateCoeffParams{1};
  
  % evaluate rate coefficient
  rateCoeff = 0;
  
  %Polyquantum
  if (NLevels-w)>0
      rateCoeff = 4.0d-16*(Tg/300)^0.5/NLevels;  %m3s-1
        
      aux = 0.105*E10N*NLevels*(1-(deltaEN/E10N)*(NLevels-1))/kb;
      if aux<=7500
          rateCoeff = rateCoeff*exp((-7500+aux)/Tg);
      end
  else
      error(['nitrogenAtomicHVT can not evaluate the rate coefficient for the reaction:\n%s\n' ...
      'The reaction does not comply with the following structure:\n v + X <-> v-1 + X;\n'], ...
      reaction.description);
  end
  
  % Monoquantum
  if (NLevels-w)<2
      aux = 1-1.26d-2*(NLevels-1);
      rateCoeff = rateCoeff + (NLevels*(1+(NLevels-1)*deltaEN/E10N))*aux^(8/3)*110*exp((-122.5*aux^0.681*Tg^(-1/3)));      
  end
  
  rateCoeff = rateCoeff*exp(hbar*omegaN*(w-NLevels+chiN*(NLevels*(NLevels+1)-w*(w+1)))/(kb*Tg));
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  
end