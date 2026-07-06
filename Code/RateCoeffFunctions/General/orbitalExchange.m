function [rateCoeff, dependent] = orbitalExchange(~, ~, ~, reaction, rateCoeffParams, chemistry)
  % orbitalExchange evaluates the rate coefficient for quantum orbital-number exchange reactions using the expression
  %  C = a * (Delta u_ij/Tg)^b * exp[-(Delta u_ij)*c/Tg] / gi , Delta u_ij > T_g 
  %  C = a * (Delta u_ij/Tg)^b / gi , Delta u_ij < T_g 
  % M.Santos et al 2014 J. Phys. D: Appl. Phys, 47, 265201
  
  Tg = chemistry.workCond.gasTemperature*Constant.boltzmannInEV; % in eV  
  
  ui = reaction.reactantArray.energy;
  gi = reaction.reactantArray.statisticalWeight;
  uj = reaction.productArray.energy;
  duij = ui-uj;
  
  a = rateCoeffParams{1};
  b = rateCoeffParams{2};
  
  if duij < Tg
   rateCoeff = a * (duij/Tg)^(b) / gi;
  else
   c = rateCoeffParams{3};
   rateCoeff = a * (duij/Tg)^(b) * exp(-duij*c/Tg)/gi;
  end
 
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', true, 'onElectronKinetics', false);

end
