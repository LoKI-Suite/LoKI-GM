function [rateCoeff, dependent] = eedfCopyv0(~, ~, ~, reaction, ~, ~)
% eedfCopyv0 evaluates the rate coefficient of a reaction involving vibrational states 
%  (either as reactants or products), 
%  by copying the rate coefficient of a similar reaction involving the corresponding parent state
  
  % Check if the reactants are in a vibrational state and define the corresponding parent
  % Identify equivalent targets in the electron kinetics 
  if isempty(reaction.reactantArray.vibLevel) && all( cellfun(@isempty, {reaction.productArray.vibLevel}) )
    error(['No reactants nor products in excited vibrational state in reaction:\n%s.'...
      '\nPlease change the type of reaction or use a different rateCoeff function'], ...
      reaction.descriptionExtended);
  end
  if ~isempty(reaction.reactantArray.vibLevel)
    if isempty(reaction.reactantArray.parent.eedfEquivalent)
      error(['Could not find target:\n%s\nin the electron kinetics module (LXCatFiles or LXCatFilesExtra).'...
        '\nPlease change the target or provide the corresponding cross section'], ...
        reaction.reactantArray.parent);
    end
    reactantEquiv = reaction.reactantArray.parent.eedfEquivalent;
  else
    if isempty(reaction.reactantArray.eedfEquivalent)
      error(['Could not find target:\n%s\nin the electron kinetics module (LXCatFiles or LXCatFilesExtra).'...
        '\nPlease change the target or provide the corresponding cross section'], ...
        reaction.reactantArray);
    end
    reactantEquiv = reaction.reactantArray.eedfEquivalent;
  end
  
  % Check if the products are in a vibrational state and define the corresponding parent
  % Identify equivalent products in the electron kinetics 
  for i = 1:length(reaction.productArray)
    if ~isempty(reaction.productArray(i).vibLevel)
      if isempty(reaction.productArray(i).parent.eedfEquivalent)
        error(['Could not find product:\n%s\nin the electron kinetics module (LXCatFiles or LXCatFilesExtra).'...
          '\nPlease change the product or provide the corresponding cross section'], ...
          reaction.productArray(i).parent);
      end
      productArrayEquiv(i) = reaction.productArray(i).parent.eedfEquivalent;
    else
      if isempty(reaction.productArray(i).eedfEquivalent)
        error(['Could not find product :\n%s\nin the electron kinetics module (LXCatFiles or LXCatFilesExtra).'...
          '\nPlease change the product or provide the corresponding cross section'], ...
          reaction.productArray(i));
      end
      productArrayEquiv(i) = reaction.productArray(i).eedfEquivalent;
    end
  end
  
  % Identify equivalent reaction (with parent states) in the electron kinetics 
  [collisionID, equivColl] = Collision.findEquivalent(reactantEquiv, productArrayEquiv, ...
    reaction.productStoiCoeff, reaction.isReverse);
  if collisionID == -1
    error(['Could not find reaction:\n%s\nin the electron kinetics module (LXCatFiles or LXCatFilesExtra).'...
      '\nPlease change the type of reaction or provide the corresponding cross section'], ...
      reaction.descriptionExtended);
  end
  
  % get rate coefficient (including superelastic rate coefficient if needed)
  rateCoeff = equivColl.ineRateCoeff;
  if reaction.isReverse
    reaction.backRateCoeff = equivColl.supRateCoeff;
  end
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end
