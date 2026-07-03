function EoN = getEoNFromCurrent(I, chemistry)
% getEoNFromCurrent provides the E/N corresponding to a given discharge current.
% This function is used by reduced-field pulse functions based on current imposition

  persistent redFieldValues;
  persistent electronTemperatureValues;
  persistent redMobilityValues;
  persistent driftVelocityValues;
  persistent chamberRadius;

  if I == 0
    EoN = 0;
    return;
  end  
  
  % load values from the look-up table for swarm parameters
  if isempty(redFieldValues) && isempty(electronTemperatureValues)
    if strcmp(chemistry.lookUpMethod, "localField")
      % save local copy of reduced field values
      redFieldValues = chemistry.lookUpTableRedFieldValues;
      % save local copy of drift velocity values
      driftVelocityValues = chemistry.lookUpTableDriftVelocityValues;
    elseif strcmp(chemistry.lookUpMethod, "localEnergy")
      % save local copy of electron temperature values
      electronTemperatureValues = chemistry.lookUpTableEleTempValues;
      % save local copy of reduced mobility values
      redMobilityValues = chemistry.lookUpTableRedMobValues;
    else
      error('Invalid look-up method. Please select either "localField" or "localEnergy".');
    end
    % save local copy of the chamber radius (does not change during the simulation)
    if chemistry.workCond.chamberRadius == 0
      error('Chamber radius cannot be zero for current imposition. Please check the working conditions.');
    else
      chamberRadius = chemistry.workCond.chamberRadius;
    end
  end
  
  % determine drift velocity from the discharge current
  vd = I/(Constant.electronCharge*chemistry.workCond.electronDensity*pi*chamberRadius^2);
  
  if strcmp(chemistry.lookUpMethod, "localField")
    % get the corresponding EoN from linear interpolation
    EoN = interp1(driftVelocityValues, redFieldValues, vd, "linear");
    if isnan(EoN)
      error(sprintf(['The drift velocity (%f m/s) deduced from discharge current is out of the range of values ' ...
        'found in the lookUpTable for the swarm parameters.\nPlease check the provided lookUpTable.'], vd),1); 
    end
  elseif strcmp(chemistry.lookUpMethod, "localEnergy")
    redMob = interp1(electronTemperatureValues, redMobilityValues, chemistry.workCond.electronTemperature, "linear");
    if isnan(redMob)
      error(sprintf(['The electron temperature (%f eV) in the working conditions is out of the range of values ' ...
        'found in the lookUpTable for the swarm parameters.\nPlease check the provided lookUpTable.'], ...
        chemistry.workCond.electronTemperature),1); 
    end
    EoN = vd/redMob*1E21;
  end  
end 