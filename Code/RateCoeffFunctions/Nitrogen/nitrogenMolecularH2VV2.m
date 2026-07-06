function [rateCoeff, dependent] = nitrogenMolecularH2VV2(~, ~, ~, reaction, ~, chemistry)
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
    if levels(1,2)-levels(1,1)==2 && levels(2,2)-levels(2,1)==-1
      w = levels(1,2);
      v = levels(2,1);
%    elseif levels(2,2)-levels(1,1)==1 && levels(1,2)-levels(2,1)==-1
%      w = levels(2,2);
%      v = levels(2,1);
    elseif levels(1,2)-levels(1,1)==-2 && levels(2,2)-levels(2,1)==1
      v = levels(1,1);
      w = levels(2,2);
%    elseif levels(2,2)-levels(1,1)==-1 && levels(1,2)-levels(2,1)==1
%      v = levels(1,1);
%      w = levels(1,2);
    else  % check for two quanta jump between vibrational levels
      error(['nitrogenMolecularVV can not evaluate the rate coefficient for the reaction:\n%s\n' ...
        'The reaction does not comply with the following structure:\n v + w-1 -> v-2 + w; for w<v+1\n'], ...
        reaction.description);
    end
  end
 % if w>v  % check the directionality of the reaction
 %  error(['nitrogenMolecularVV can not evaluate the rate coefficient for the reaction:\n%s\n' ...
 %     'The reaction does not comply with the following structure:\n v + w-1 -> v-1 + w; for w<v+1\n'], ...
 %     reaction.description);
 % end
  
 % intermediate calculations
 E10H = hbar*omegaH*(1-(2*chiH));
 E10N = hbar*omegaN*(1-(2*chiN));
 E20N = hbar*omegaN*(2-(6*chiN));
 deltaEH = hbar*omegaH*chiH;
 deltaEN = hbar*omegaN*chiN;
 Q0=2.4d-21*(Tg/300.0d0)^1.5d0;    %m3s-1
 
 EH= E10H*(w-1)*(1-(deltaEH/E10H)*(w-2));
 EH= E10H*w*(1-(deltaEH/E10H)*(w-1)) - EH;
% EN= E10N*(v-2)*(1-(deltaEN/E10N)*(v-3))
% EN= E10N*v*(1-(deltaEN/E10N)*(v-1)) - EN
 EN= E10N*v*(1-(deltaEN/E10N)*(v-1));
 EN= E10N*(v+2)*(1-(deltaEN/E10N)*(v+1)) - EN;
 
 dH = w*(1+(w-1)*deltaEH/E10H);
% dN = (v-1)*(1+(v-2)*deltaEN/E10N)
  dN = (v+1)*(1+v*deltaEN/E10N);
 dN2= 0.5d0*dN*(dN+1);
 Denr=EH-EN;
 F = exp(-37/Tg^0.5*abs((EN-EH)/(E10H-E20N)));

 rateCoeff = Q0*dH*dN2*F*(1.5-0.5*F);
 if (Denr<0)
      rateCoeff = rateCoeff*exp(Denr/(kb*Tg)); %tirei abs(denr)
 end             
 
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
  
end
