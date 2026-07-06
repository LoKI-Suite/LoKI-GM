function [rateCoeff, dependent] = inFlow(~, ~, ~, reaction, ~, chemistry)
% inFlow evaluates the (generalized) rate coefficient (in m-3 s-1) for the inlet flow of species 
%  assuming conservation of atoms in the gas/plasma mixture

% the in/out flow is automatically handled by the code according to the working conditions
%  totalSccmInFlow (value of the inFlow in sccm)
%  totalSccmOutFlow (value of the outFlow in sccm, or model that defines the outFlow)
%  inFlowFraction
%  inFlowPopulation

  % calculate the persistent variables
  persistent chamberVolume;
  persistent rateCoeffAux;

   if isempty(chamberVolume)
 	 chamberVolume = pi*(chemistry.workCond.chamberRadius^2)*chemistry.workCond.chamberLength;
     rateCoeffAux = chemistry.workCond.totalSccmInFlow*Constant.sccmToParticleRateCoeff;
  end  

  % calculate the inflow rate coefficient (in m-3 s-1)
  %  rate coefficient (m-3 s-1) = Qin * factor * delta / V
  %   Qin = inflow (sccm)
  %   factor = Constant.atmosphereInPa/Constant.boltzmann*1e-6/60/273.15 (cm-3/(s min-1)
  %   delta = inflow_populations
  %   V = volume (m3)
  %   Note: [Qin * factor] (s-1)
  rateCoeff = rateCoeffAux*reaction.productArray.inFlowRelDensity/chamberVolume;
  
  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', false, 'onGasTemperature', false, 'onElectronKinetics', false);
end
