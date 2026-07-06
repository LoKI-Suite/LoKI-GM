function [rateCoeff, dependent] = pseudoEedfCopy(~, ~, ~, reaction, rateCoeffParams, chemistry)
% pseudoEedfCopy evaluates a rate coefficient of an electron-impact reaction by setting it equal to 
%  the rate coefficient of another reaction using the reaction's description as input parameter

  persistent equivalentCollision;
  
  if isempty(equivalentCollision)
    equivalentCollision = cell((length(chemistry.reactionArray)));
  end  

  reactionID = reaction.ID;

  if isempty(equivalentCollision{reactionID})
    % reconstruct the collision description
    % e.g. when e + O2(X,v=0) -> e + e + O2(+,X) is given as argument
    % the loop concatenates 'e + O2(X' + ',' + 'v=0) -> e + e + O2(+' + ',' + 'X)'
    collisionDescription = rateCoeffParams{1};
    for i=2:length(rateCoeffParams)
      collisionDescription = [collisionDescription ',' rateCoeffParams{i}];
    end
    % remove spaces
    collisionDescription = collisionDescription(~isspace(collisionDescription));
    % find the equivalent collision in the target gas array
    targetGas = reaction.reactantArray(1).gas.eedfEquivalent;
    for collision = [targetGas.collisionArray targetGas.collisionArrayExtra]
      % get the candidate full description (collision,type)
      description = collision.description;
      % remove the type part
      commasPositions = strfind(description,',');
      description = description(1:(commasPositions(end)-1));
      % check if the candidate description coincides
      if strcmp(description, collisionDescription)
        equivalentCollision{reactionID} = collision;
        break;
      end  
    end  
    % check if the equivalent collision was found
    if isempty(equivalentCollision{reactionID})
      error('Could not find collision:\n%s\nNeeded to evaluate ''pseudoEedfCopy'' for reaction\n%s',...
        collisionDescription,reaction.description);
    elseif reaction.isReverse && ~equivalentCollision{reactionID}.isReverse
      error(['Error when using ''pseudoEedf_copy'' for reaction\n%s',...
        '\nThe collision to be copied has no reverse.'],reaction.description);
    end  
  end  

  % save the forward rate coeff
  rateCoeff = equivalentCollision{reactionID}.ineRateCoeff;
  % save the backwards rate coeff, if it has inverse
  if reaction.isReverse
    reaction.backRateCoeff = equivalentCollision{reactionID}.supRateCoeff;
  end  

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', true);
end
