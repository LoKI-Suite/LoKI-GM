function [rateCoeff, dependent] = nitrogenMolecularVVDis(~, ~, ~, reaction, ~, chemistry)
% nitrogenMolecularVVDis evaluates the dissociation rate coefficient 
% for V-V transitions between nitrogen molecules described in:
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
  persistent Lmin;
  persistent M;
  persistent omega;
  persistent chi;
  persistent levelsNeedToBeIdentified;
  persistent w;
  persistent NLevels;
  
  % --- performance sensitive calculations ---
  
  % initialize variables the first time the classicalAmbipolar function is called
  if isempty(dependentInfo)
    dependentInfo = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);
    kb = Constant.boltzmann;
    hbar = Constant.planckReduced;
    Lmin = 2e-11;                                                                   % Minimum distance between N2 molecules during V-V collision (m)
    M = reaction.reactantArray(1).gas.mass;                        % Mass of N2 (Kg)
    omega = reaction.reactantArray(1).gas.harmonicFrequency;       % Harmonic frequency of the oscillator N2 (rad/s)
    chi = reaction.reactantArray(1).gas.anharmonicFrequency/omega; % Anharmonicity of the oscillator N2
    reactionArraySize = size(chemistry.reactionArray);
    levelsNeedToBeIdentified = true(reactionArraySize);
    w = zeros(reactionArraySize);
  end
  
  % --- time independent calculations (also performance sensitive) ---
  reactionID = reaction.ID;
  if levelsNeedToBeIdentified(reactionID)
    % identify current V-V reaction by evaluating the vibrational levels ( w + Nlevels-1 -> w-1 + 2Atoms )
    if reaction.reactantStoiCoeff(1) == 2
      w(reactionID) = str2double(reaction.reactantArray(1).vibLevel);
      if isempty(NLevels)
        NLevels = w(reactionID);
      end
    else
      levels(2) = str2double(reaction.reactantArray(2).vibLevel);
      levels(1) = str2double(reaction.reactantArray(1).vibLevel);
      w(reactionID) = min(levels);
      if isempty(NLevels)
        NLevels = max(levels)+1;
      end
    end
    % --- ERROR CHECKING NEEDS TO BE ADDED ---
    
    % once identified the levels of the reaction set the boolean variable to false
    levelsNeedToBeIdentified(reactionID) = false;
  end
  
  % --- regular calculations ---
  
  % local copy of gas temperature
  Tg = chemistry.workCond.gasTemperature;
  % intermediate calculations
  Y=pi*Lmin*omega*chi*sqrt(M/(kb*Tg))*(NLevels-w(reactionID));
%   Y=9.97*sqrt(1/Tg)*(NLevels-w(reactionID)); % LoKI1.2.0 expression
  if Y<=20
    F=0.5*(3-exp(-Y/1.5))*exp(-Y/1.5);
  else
    F=8*sqrt(pi/3)*Y^(7/3)*exp(-3*Y^(2/3));
%     F=8.1845*Y^(7/3)*exp(-3*Y^(2/3)); % LoKI1.2.0 expression
  end
  % Billing corrections
  if NLevels < 10 % LoKI1.2.0 uses <= here (error check paper in the header)
    kor=39.0625-1.5625*NLevels;
  else
    kor=25.2+24.1*((NLevels-10)/30)^3;
  end
  
  % evaluate rate coefficient
  rateCoeff=(NLevels/(1-chi*NLevels))*(w(reactionID)/(1-chi*w(reactionID)))*(1-chi)^2*F*6.354e-23*Tg^(1.5)/kor*...
    exp(-2*hbar*omega*chi*(NLevels-w(reactionID))/(kb*Tg));
%   rateCoeff=(NLevels/(1-chi*NLevels))*(w(reactionID)/(1-chi*w(reactionID)))*F*6.354e-23*Tg^(1.5)/kor*...
%     exp(-2*hbar*omega*chi*(NLevels-w(reactionID))/(kb*Tg)); % LoKI1.2.0 expression (error, missing (1-chi)^2 factor, check paper in the header)
  
  % set function dependencies
  dependent = dependentInfo;
  
end
