function [rateCoeff, dependent] = outFlow(time, ~, totalGasDensity, ~, rateCoeffParams, chemistry)
% outFlow evaluates the rate coefficient (in s-1) for the outlet flow of species 
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
    rateCoeffAux = Constant.sccmToParticleRateCoeff/chamberVolume;
  end            
  
  % calculate the outflow rate coefficient (in s-1)
  if strcmp(rateCoeffParams{1}, 'ensureIsobaric') && time == 0
    % for flow model 'ensureIsobaric', the rate coefficient is set to zero in the first iteration
    % (to solve the thermal model, where the flow is assumed not to contribute)
    rateCoeff = 0;
  else
    % for other flow models, the rate coefficient varies with the time-dependent totalGasdensity 
    %  rate coefficient (s-1) = Qout * factor / V / N(t)
    %   Qout = outflow (sccm)
    %   factor = Constant.atmosphereInPa/Constant.boltzmann*1e-6/60/273.15 (cm-3/(s min-1)
    %   V = volume (m3)
    %   N(t) = time-dependent total-gas-density (m-3)
    %   Note: [Qout * factor] (s-1)
    %
    % In cases of high dissociation, N(t) tends to increase (at constant gas temperature) 
    %  and the outflow rate coefficient decreases, 
    %  potentially leading to an increase in N(t_final), 
    %  and therefore to a limitation in the time-convergence of species densities
    % In these cases, it is advisable to control/restrict the total simulation time 
    %  or to adopt the flow model 'ensureIsobaric'
    %
    %
    % for flow model 'ensureIsobaric' and time<>0, the rate coefficient is calculated 
    %  using the workCond.totalSccmOutFlow self-consistently obtained by solving the equation 
    %  dp/dt = (dN_chem/dt + dN_outflow/dt) k_B Tg + (N_chem + N_outflow) k_B (dTg/dt) = 0
    %  (note that, strickly, the previous equation is valid in steady-state only)
    % 
    %
    % the flowBarrier factor (between 0-1) affects the outflow of each individual species
    %  for flow models other than 'ensureIsobaric'  
    if isempty(rateCoeffParams{2}) || strcmp(rateCoeffParams{1}, 'ensureIsobaric')
      flowBarrier = 1;  
    else
      flowBarrier = rateCoeffParams{2};
    end  
    rateCoeff = flowBarrier*rateCoeffAux*chemistry.workCond.totalSccmOutFlow/totalGasDensity;
  end 

  % set function dependencies
  dependent = struct('onTime', false, 'onDensities', true, 'onGasTemperature', false, 'onElectronKinetics', false);
  
end