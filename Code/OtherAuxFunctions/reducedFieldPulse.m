function EoN = reducedFieldPulse(time, parameters)
%reducedFieldPulse returns the value of the reduced electric field (in Td) during the pulse 
  
  charTime = parameters{1};
  amplitude = parameters{2};
  EoN = amplitude*sqrt(time/charTime)*exp(-time/charTime);
  
end

