function [rateCoeff, dependent] = nitrogenEVscaling(~, ~, ~, reaction, ~, chemistry)
% nitrogenEVscaling evaluates the rate coefficients for eV transitions in nitrogen following the scaling law described in
% A. Bourdon, P. Vervisch, J. Thermophys. Heat Transf. 14 (2000) 489-495
% doi:10.2514/2.6571.

% see also
% V. Guerra, A. Tejero-del-Caz, C. D. Pintassilgo and L. L. Alves
% Modelling N2–O2 plasmas: volume and surface kinetics, Plasma Sources Sci. Technol. 28 073001 (2019)
% https://iopscience.iop.org/article/10.1088/1361-6595/ab252c

  persistent equivalentReactionID;

  vini = str2double(reaction.reactantArray.vibLevel);            % obtain initial vibrational level
  vfin = str2double(reaction.productArray.vibLevel);             % obtain final vibrational level
  n = vfin - vini;                                                      % evaluate vibrational quanta jump (positive or negative)
  equivReactionDscrp = sprintf('e + N2(X,v=0) <-> e + N2(X,v=%d)', abs(n));   % evaluate equivalent reaction description

  % local copy of the reaction array
  reactionArray = chemistry.reactionArray;

  % find equivalent reaction
  reactionID = reaction.ID;
  if isempty(equivalentReactionID) || length(equivalentReactionID)<reactionID || equivalentReactionID(reactionID)==0
    for idx=1:length(reactionArray)
      if strcmp(reactionArray(idx).description, equivReactionDscrp)
        equivalentReactionID(reactionID) = idx;
        break;
      elseif idx == length(reactionArray)
        error('Could not find reaction: %s\nNeeded to evaluate ''nitrogenEVscaling'' for reaction %s\n', ...
          equivReactionDscrp, reactionArray(reactionID).description);
      end
    end
  end

  % scale equivalent rate coefficient
  ID = equivalentReactionID(reactionID);
  if n>0   % inelastic collision
    rateCoeff = reactionArray(ID).eedfEquivalent.ineRateCoeff/(1+0.15*vini);
  else     % superelastic collision
    rateCoeff = reactionArray(ID).eedfEquivalent.supRateCoeff/(1+0.15*vfin);
  end

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end
