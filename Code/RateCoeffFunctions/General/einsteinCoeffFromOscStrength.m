function [rateCoeff, dependent] = einsteinCoeffFromOscStrength(~, ~, ~, reaction, rateCoeffParams, ~)
% einsteinCoeffFromOscStrength evaluates the Einstein coefficient for a radiative transition using the expression:
%  C = 2 * pi * epsilon_0 * e^4 / (m_e*h^2*c) * (Delta u_ij)^2 * (g_j/g_i) * f_ji

  a = Constant.vacuumPermeability*2*pi*Constant.electronCharge^4/(Constant.electronMass*...
      Constant.planck^2*Constant.speedOfLight);
  fji = rateCoeffParams{1};
  ui = reaction.reactantArray.energy;
  gi = reaction.reactantArray.statisticalWeight;
  uj = reaction.productArray.energy;
  gj = reaction.productArray.statisticalWeight;

  du_ij = ui - uj;

  rateCoeff = a * du_ij^2 * gj/gi * fji;
  
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);

end
