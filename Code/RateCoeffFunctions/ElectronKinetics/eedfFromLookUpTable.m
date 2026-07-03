function [rateCoeff, dependent] = eedfFromLookUpTable(~, ~, ~, reaction, ~, chemistry)
% eedfFromLookUpTable obtains a reaction rate coefficient by interpolating the value from a look-up table of rate 
% coefficients as a function of the reduced electric field or the electron temperature, depending on the local 
% approximation method. The look-up table is defined in the input file and should be consistent with the reaction 
% descriptions in the reaction array. The function also saves the backwards rate coefficient for reversible reactions, 
% if applicable.
  
  persistent rateCoeffValues;
  
  if isempty(rateCoeffValues)
    rateCoeffValues = cell((size(chemistry.reactionArray)));
  end  

  if isempty(rateCoeffValues{reaction.ID})
    descriptionReaction = reaction.description;
    rateCoeffTable = [];
    % remove spaces of description
    descriptionReaction = descriptionReaction(~isspace(descriptionReaction));
    for column = chemistry.lookUpTableRateCoeff
      variableDescription = column.Properties.VariableDescriptions{1};
      if strcmp(variableDescription(1:end-4), descriptionReaction)
        if strcmp(variableDescription(end-2:end), 'ine')
          rateCoeffTable = [rateCoeffTable column.Variables];
        elseif strcmp(variableDescription(end-2:end), 'sup')
          rateCoeffTable = [rateCoeffTable column.Variables];
          break;
        end
      end
    end
    % error if no collision description has been found
    if isempty(rateCoeffTable)
      error('Collision description not found:\n %s',descriptionReaction);
    end
    rateCoeffValues{reaction.ID} = rateCoeffTable;
  end  

  if strcmp(chemistry.lookUpMethod, 'localField')
    % save the forward rate coeff
    rateCoeff = interp1(chemistry.lookUpTableRedFieldValues,rateCoeffValues{reaction.ID}(:,1), ...
      chemistry.workCond.reducedField,'linear');
    if isnan(rateCoeff)
      error(['The reduced electric field value is out of the look-up table range for reaction %d.\n' ...
        'Please check the setup file.'], reaction.description);
    end
    % save the backwards rate coeff, if it has inverse
    if reaction.isReverse
      reaction.backRateCoeff = interp1(chemistry.lookUpTableRedFieldValues,rateCoeffValues{reaction.ID}(:,2), ...
        chemistry.workCond.reducedField,'linear'); 
        if isnan(reaction.backRateCoeff)
          error(['The reduced electric field value is out of the look-up table range for reaction %d.\n' ...
            'Please check the setup file.'], reaction.description);
        end
    end 
  elseif strcmp(chemistry.lookUpMethod, 'localEnergy')
    % save the forward rate coeff
    rateCoeff = interp1(chemistry.lookUpTableEleTempValues,rateCoeffValues{reaction.ID}(:,1), ...
      chemistry.workCond.electronTemperature,'linear');
    if isnan(rateCoeff)
      error(['The electron temperature value is out of the look-up table range for reaction %d.\n' ...
        'Please check the setup file.'], reaction.description);
    end
    % save the backwards rate coeff, if it has inverse
    if reaction.isReverse
      reaction.backRateCoeff = interp1(chemistry.lookUpTableEleTempValues,rateCoeffValues{reaction.ID}(:,2), ...
        chemistry.workCond.electronTemperature,'linear'); 
        if isnan(reaction.backRateCoeff)
          error(['The electron temperature value is out of the look-up table range for reaction %d.\n' ...
            'Please check the setup file.'], reaction.description);
        end
    end  
  else
    error('Invalid look-up method for electron kinetics. Please check the setup file.');
  end  
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);

end
