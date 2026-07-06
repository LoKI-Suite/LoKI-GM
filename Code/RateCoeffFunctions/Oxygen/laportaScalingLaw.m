function [rateCoeff, dependent] = laportaScalingLaw(time, densitiesAll, totalGasDensity, reaction, rateCoeffParams, chemistry)
% laportaScalingLaw evaluates the rate coefficient for EID1 and DEA reactions described in:
% Adriana Annusova et al 2018 Plasma Sources Sci. Technol. 27 045006
% https://doi.org/10.1088/1361-6595/aab47d

  %EID1 - electron impact dissociation (R27)
  %DEA - dissociative electron attachment (R31)

  persistent delta;
  persistent IDEID;
  persistent IDDEA;

  %Table A1: Scaling parameter, delta(v), for the processes of EID1 and DEA from equation (13).
           % EID1      DEA       v
  if isempty(delta)
    delta = [0.2872    0.4493    % v =  1
             0.2337    0.3423    % v =  2
             0.1909    0.2659    % v =  3
             0.1664    0.2165    % v =  4
             0.1479    0.1810    % v =  5
             0.1316    0.1546    % v =  6
             0.1181    0.1342    % v =  7
             0.1071    0.1184    % v =  8
             0.0979    0.1057    % v =  9
             0.0901    0.0958    % v = 10
             0.0832    0.0878    % v = 11
             0.0773    0.0809    % v = 12
             0.0721    0.0750    % v = 13
             0.0675    0.0699    % v = 14
             0.0635    0.0654    % v = 15
             0.0598    0.0613    % v = 16
             0.0566    0.0577    % v = 17
             0.0537    0.0545    % v = 18
             0.0511    0.0517    % v = 19
             0.0487    0.0490    % v = 20
             0.0464    0.0466    % v = 21
             0.0445    0.0445    % v = 22
             0.0426    0.0425    % v = 23
             0.0409    0.0407    % v = 24
             0.0393    0.0390    % v = 25
             0.0379    0.0375    % v = 26
             0.0365    0.0361    % v = 27
             0.0352    0.0347    % v = 28
             0.0340    0.0334    % v = 29
             0.0329    0.0322    % v = 30
             0.0318    0.0311    % v = 31
             0.0309    0.0301    % v = 32
             0.0299    0.0290    % v = 33
             0.0289    0.0277    % v = 34
             0.0281    0.0269    % v = 35
             0.0273    0.0258    % v = 36
             0.0264    0.0243    % v = 37
             0.0254    0.0224    % v = 38
             0.0237    0.0174    % v = 39
             0.0228    0.0154    % v = 40
             0.0211    0.0100 ]; % v = 41
  end

  % identify vibrational level of O2 molecule
  v = str2double(reaction.reactantArray(1).vibLevel);

  flagLoopBreak = 0;

  %identify reaction: EID1 or DEA
  if strcmp(rateCoeffParams{1}, 'EID1') %EID1
    deltaParameter = delta(v,1); %get delta(v) value for EID1
    if isempty(IDEID)
      for ID=1:length(chemistry.reactionArray)
        %looking for ID of the EID1 with v=0
        if strcmp(chemistry.reactionArray(ID).description, 'e + O2(X,v=0) -> e + 2O(3P)')
          IDEID = ID;
          flagLoopBreak = 1;
          break;
        end
      end
      if flagLoopBreak == 0
        error('Reaction description given as input argument was not found');
      end
    end
    rxn = chemistry.reactionArray(IDEID);
    rateCoeffV0 = rxn.rateCoeffFuncHandle(time, densitiesAll, totalGasDensity, rxn, rxn.rateCoeffParams, chemistry);
  elseif strcmp(rateCoeffParams{1}, 'DEA') %DEA
    deltaParameter = delta(v,2); %get delta(v) value for DEA
    if isempty(IDDEA)
      for ID=1:length(chemistry.reactionArray)
        %looking for ID of the DEA reaction with v=0
        if strcmp(chemistry.reactionArray(ID).description, 'e + O2(X,v=0) -> O(-,gnd) + O(3P)')
          IDDEA = ID;
          flagLoopBreak = 1;
          break;
        end
      end
      if flagLoopBreak == 0
        error('Reaction description given as input argument was not found');
      end
    end
    rxn = chemistry.reactionArray(IDDEA);
    rateCoeffV0 = rxn.rateCoeffFuncHandle(time, densitiesAll, totalGasDensity, rxn, rxn.rateCoeffParams, chemistry);
  else
    error('Type of reaction not allowed with EID1_DEA function');
  end

  %rateCoeffV0 in cm^3/s or m^3/s?
  rateCoeff = rateCoeffV0/(1-deltaParameter*v);

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);
end
