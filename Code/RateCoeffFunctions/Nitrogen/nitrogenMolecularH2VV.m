function [rateCoeff, dependent] = nitrogenMolecularH2VV(~, ~, ~, reaction, ~, chemistry)
% nitrogenMolecularVV evaluates the rate coefficient for V-V transitions between nitrogen molecules described in
% J Loureiro and C M Ferreira, J. Phys. D: Appl. Phys. 19 (1986) 17-35
% (CHANGE!!!!)
% http://dx.doi.org/10.1088/0022-3727/19/1/007

  % local copies and definitions of different parameters used in the function
  kb = Constant.boltzmann;
  hbar = Constant.planckReduced;
  Tg = chemistry.workCond.gasTemperature;
  uma = Constant.unifiedAtomicMass; 
  Lmin = 3.0e-11; % Minimum distance between N2 molecules during V-V collision (m)
  dmin=3.34e-10;
  MN2 = reaction.reactantArray(1).gas.mass;  % Mass of N2 (Kg)
  MH2 = reaction.reactantArray(2).gas.mass;  % Mass of H2 (Kg)

  omegaN = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
  chiN = reaction.reactantArray(1).gas.anharmonicFrequency/omegaN; % Anharmonicity of the oscillator N2
  omegaH = reaction.reactantArray(2).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
  chiH = reaction.reactantArray(2).gas.anharmonicFrequency/omegaH; % Anharmonicity of the oscillator N2
  
  % identify current V-V reaction by evaluating the vibrational levels ( v + w-1 -> v-1 + w )
  if reaction.reactantStoiCoeff(1) == 2
    v = str2double(reaction.reactantArray(1).vibLevel);
    w = v+1;
  elseif reaction.productStoiCoeff(1) == 2
    w = str2double(reaction.productArray(1).vibLevel);
    v = w+1;
  else
    levels = zeros(2);
    levels(1,1) = str2double(reaction.reactantArray(1).vibLevel);
    levels(2,1) = str2double(reaction.reactantArray(2).vibLevel);
    levels(1,2) = str2double(reaction.productArray(1).vibLevel);
    levels(2,2) = str2double(reaction.productArray(2).vibLevel);
    if levels(1,2)-levels(1,1)==1 && levels(2,2)-levels(2,1)==-1
      w = levels(1,2);
      v = levels(2,1);
    elseif levels(2,2)-levels(1,1)==1 && levels(1,2)-levels(2,1)==-1
      w = levels(2,2);
      v = levels(2,1);
    elseif levels(1,2)-levels(1,1)==-1 && levels(2,2)-levels(2,1)==1
      v = levels(1,1);
      w = levels(2,2);
    elseif levels(2,2)-levels(1,1)==-1 && levels(1,2)-levels(2,1)==1
      v = levels(1,1);
      w = levels(1,2);
    else  % check for one quanta jump between vibrational levels
      error(['nitrogenMolecularVV can not evaluate the rate coefficient for the reaction:\n%s\n' ...
        'The reaction does not comply with the following structure:\n v + w-1 -> v-1 + w; for w<v+1\n'], ...
        reaction.description);
    end
  end
  
%  if w>v  % check the directionality of the reaction
%    error(['nitrogenMolecularVV can not evaluate the rate coefficient for the reaction:\n%s\n' ...
%      'The reaction does not comply with the following structure:\n v + w-1 -> v-1 + w; for w<v+1\n'], ...
%      reaction.description);
%  end
  
 % intermediate calculations
  uM = 1.0/((1/MH2)+(1/MN2));  %Reduced mass
  uH = 0.25*MH2;
  uN = 0.25*MN2;

  ZNH = pi*dmin*dmin*((8.0d0*kb*Tg)/(pi*uM))^0.5;
  U2N =  (1.0d0/(8.0d0*Lmin*Lmin*uN*omegaN));
  U2H =  (1.0d0/(8.0d0*Lmin*Lmin*uH*omegaH));
  P10= ZNH*U2N*U2H*8.0d0*Lmin*Lmin*uM*kb*Tg;  %*1.0d6  %m3s-1
  
  DE=hbar*(omegaN*(1-(2*chiN*v))-omegaH*(1-(2*chiH*w)));
  Y=abs(DE)*pi*Lmin*(uM/(2*kb*Tg))^0.5/hbar;
  
  if Y<=20
    F=0.5*(3-exp(-Y/1.5))*exp(-Y/1.5);
  else
    F=8*sqrt(pi/3)*Y^(7/3)*exp(-3*Y^(2/3));
  end
   
  % evaluate rate coefficient
  rateCoeff = (v/(1-chiN*v))*(w/(1-chiH*w))*P10*F;
  if (DE>0)
      rateCoeff = rateCoeff*exp(abs(DE)/(kb*Tg));
  end  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  
end
