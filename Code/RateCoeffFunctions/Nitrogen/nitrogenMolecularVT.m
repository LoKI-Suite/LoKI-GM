function [rateCoeff, dependent] = nitrogenMolecularVT(time, densitiesAll, totalGasDensity, reaction, rateCoeffParams, chemistry)
% nitrogenMolecularVT evaluates the rate coefficient for V-T transitions between nitrogen molecules described in:
% L.L. Alves, L. Marques, C.D. Pintassilgo, G. Wattieaux, E. Es-sebbar, J. Berndt, et al.
% Capacitively coupled radio-frequency discharges in nitrogen at low pressures, Plasma Sources Sci. Technol. 21 (2012) 45008.
% http://dx.doi.org/10.1088/0963-0252/21/4/045008

% see also
% V. Guerra, A. Tejero-del-Caz, C. D. Pintassilgo and L. L. Alves
% Modelling N2–O2 plasmas: volume and surface kinetics, Plasma Sources Sci. Technol. 28 073001 (2019)
% https://iopscience.iop.org/article/10.1088/1361-6595/ab252c
  
  % persistent variables for performance reasons
  persistent dependentInfo;
  persistent kb;
  persistent Lmin;
  persistent M;
  persistent omega;
  persistent chi;
  persistent levelsNeedToBeIdentified;
  persistent v;
  
  % --- performance sensitive calculations ---
  
  % initialize variables the first time the classicalAmbipolar function is called
  if isempty(dependentInfo)
    dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
    kb = Constant.boltzmann;
    Lmin = 2e-11;                                                                   % Minimum distance between N2 molecules during V-V collision (m)
    M = reaction.reactantArray(1).gas.mass;                        % Mass of N2 (Kg)
    omega = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
    chi = reaction.reactantArray(1).gas.anharmonicFrequency/omega; % Anharmonicity of the oscillator N2
    reactionArraySize = size(chemistry.reactionArray);
    levelsNeedToBeIdentified = true(reactionArraySize);
    v = zeros(reactionArraySize);
  end
  
  % --- time independent calculations (also performance sensitive) ---
  reactionID = reaction.ID;
  % identify current V-T reaction by evaluating the vibrational levels ( N2(X,v) + N2 <-> N2(X,v-1) + N2) 
  if levelsNeedToBeIdentified(reactionID)
    v(reactionID) = str2double(reaction.reactantArray(1).vibLevel);
    % once identified the levels of the reaction set the boolean variable to false
    levelsNeedToBeIdentified(reactionID) = false;
  end
  
  % --- regular calculations ---
  
  % local copies and definitions of different parameters used in the function
  Tg = chemistry.workCond.gasTemperature;
  
  % intermediate calculations
  Y=0.5*pi*Lmin*omega*sqrt(M/(kb*Tg))*(1-2*chi*v(reactionID));
%   Y=809.708*(1-2*v*chi)/sqrt(Tg); % LoKI1.2.0 expression
  if Y<=20
    F=0.5*(3-exp(-Y/1.5))*exp(-Y/1.5);
  else
    F=8*sqrt(pi/3)*Y^(7/3)*exp(-3*Y^(2/3));      
%     F=8.1845*Y^(7/3)*exp(-3*Y^(2/3)); % LoKI1.2.0 expression
  end
  % Billing corrections
  kor=0.2772*Tg-80.32+35.5*((v(reactionID)-1)/39)^0.8;  
  
  % evaluate rate coefficient
  rateCoeff =v(reactionID)*(1-chi)/(1-chi*v(reactionID))*F*1.07e-18*Tg^(1.5)/kor;
  
  % set function dependencies
  dependent = dependentInfo;
  
end
