function [rateCoeff, dependent] = nitrogenAtomicHVT(~, ~, ~, reaction, ~, chemistry)
% nitrogenAtomicHVT evaluates the rate coefficient for V-T transitions (monoquantum and polyquantum) between nitrogen molecules and H atoms
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
  
  % identify current V-T reaction by evaluating the vibrational levels ( N2(X,v) + H <-> N2(X,w) + H, v-w<6 ) 
  v = str2double(reaction.reactantArray(1).vibLevel);
  w = str2double(reaction.productArray(1).vibLevel);
  
  % evaluate rate coefficient
  rateCoeff = 0;

  %Polyquantum
  if (v-w)>0
      rateCoeff = 4.0d-16*(Tg/300)^0.5/v;  %m3S-1  
        
      aux = 0.105*E10N*v*(1-(deltaEN/E10N)*(v-1))/kb;
      if aux<=7500
          rateCoeff = rateCoeff*exp((-7500+aux)/Tg);
      end
  else
      error(['nitrogenAtomicHVT can not evaluate the rate coefficient for the reaction:\n%s\n' ...
      'The reaction does not comply with the following structure:\n v + X <-> v-1 + X;\n'], ...
      reaction.description);
  end
  
  % Monoquantum   
  if (v-w)<2
      aux = 1-1.26d-2*(v-1);
      rateCoeff = rateCoeff +(v*(1+(v-1)*deltaEN/E10N))*aux^(8/3)*110*exp((-122.5*aux^0.681*Tg^(-1/3)))*1.0E-16;  %m3s-1    
  end
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  
end