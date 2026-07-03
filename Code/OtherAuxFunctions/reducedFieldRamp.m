function EoN = reducedFieldRamp(time, parameters)
%reducedFieldRamp returns the value of the reduced electric field (in Td) during a ramp pulse 
  
  charTime = parameters{1};
  amplitude = parameters{2};
  if time == 0
    EoN = 0;
  elseif time<charTime
    EoN = amplitude/charTime*time;
  else
    EoN = amplitude;
  end
  
end

