function EoN = reducedFieldFromCurrentFile(time, parameters)
% reducedFieldFromCurrentFile returns the reduced field resulting from the imposition of discharge current [A].
% The fileName of the matrix (time,dischargeCurrent) is given as parameter.

  persistent timeValues;
  persistent currentvalues;

  % read the matrix at the first function evaluation
  if isempty(timeValues)
    fileName = parameters{1};
    fileName = fileName(~isspace(fileName));
    table = readtable(['Input' filesep fileName],'CommentStyle','%');
    timeValues = table{:,1};
    currentvalues = table{:,2};
  end  

  % interpolate current[A] at the given time
  I = interp1(timeValues, currentvalues, time, "linear", 0);
  % calculate the E/N [Td] corresponding to this value of current
  EoN = getEoNFromCurrent(I, parameters{end});
  
end

