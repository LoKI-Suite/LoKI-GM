function [rateCoeff, dependent] = nitrogenAtomicVTDis(~, ~, ~, reaction, rateCoeffParams, chemistry)
% nitrogenAtomicVTDis evaluates the dissociation rate coefficient 
% for V-T transitions between nitrogen molecules and atoms described in:
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
  persistent hbar;
  persistent omega;
  persistent chi;
  persistent levelsNeedToBeIdentified;
  persistent w;
  
  % --- performance sensitive calculations ---
  
  % initialize variables the first time the classicalAmbipolar function is called
  if isempty(dependentInfo)
    dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
    kb = Constant.boltzmann;
    hbar = Constant.planckReduced;
    omega = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
    chi = reaction.reactantArray(1).gas.anharmonicFrequency/omega; % Anharmonicity of the oscillator N2
    reactionArraySize = size(chemistry.reactionArray);
    levelsNeedToBeIdentified = true(reactionArraySize);
    w = zeros(reactionArraySize);
  end
  
  % --- time independent calculations (also performance sensitive) ---
  reactionID = reaction.ID;
  % identify current V-T reaction by evaluating the vibrational levels ( N2(X,w) + N -> 2N + N, NLevels-w<6 ) 
  if levelsNeedToBeIdentified(reactionID)
    w(reactionID) = str2double(reaction.reactantArray(1).vibLevel);
    % once identified the levels of the reaction set the boolean variable to false
    levelsNeedToBeIdentified(reactionID) = false;
  end
  
  % --- regular calculations ---
  
  % local copies and definitions of different parameters used in the function
  Tg = chemistry.workCond.gasTemperature;
  NLevels = rateCoeffParams{1};
  
  % evaluate rate coefficient
  rateCoeff = 0;
  if NLevels-w(reactionID)<6 && NLevels>6
    % Constants for reactive reactions
    A0r=2.2116e-12/(Tg^1.4261);
    A1r=32102.83575e0/(Tg^0.81473);
    A2r=252050.05679e0/(Tg^1.044);
    % next three lines contain old values for DC comparison
%     A0r=26.1e-16*(Tg/500)^1.5; 
%     A1r=251*sqrt(500/Tg);
%     A2r=780*sqrt(500/Tg);
    rateCoeff = rateCoeff + A0r*exp(-A1r/NLevels+A2r/(NLevels*NLevels)); % reactive contribution
    if NLevels>8
      % non-reactive contribution
      A0nr=9.2403e-12/(Tg^1.63471);
      A1nr=18229.37525/(Tg^0.69843);
      A2nr=9892.96628/(Tg^0.43756);
      % next three lines contain old values for DC comparison
%       A0nr=3.16e-16*(Tg/500)^1.5;
%       A1nr=227.2*sqrt(500/Tg);
%       A2nr=625*sqrt(500/Tg);
      rateCoeff = rateCoeff + A0nr*exp(-A1nr/NLevels+A2nr/(NLevels*NLevels)); % non-reactive contribution
    end
    rateCoeff = rateCoeff*...
      exp(hbar*omega*(w(reactionID)-NLevels+chi*(NLevels*(NLevels+1)-w(reactionID)*(w(reactionID)+1)))/(kb*Tg));
  end
  
  % set function dependencies
  dependent = dependentInfo;
  
end