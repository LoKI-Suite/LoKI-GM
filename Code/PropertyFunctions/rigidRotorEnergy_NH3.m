function energy = rigidRotorEnergy_NH3(state, ~, ~)
  % rigidRotorEnergy_NH3 is a property function that evaluates the energy of rotational state J of NH3, 
  %  following the work of Dowling 1968, which shows good agreemnt with the MARVEL database 
  %  (Derzi et al 2015, Furtenbacher et al 2020)  

  % Energy values (eV) for the rotational levels of NH3
  %  (JK)   Marvel database      this function
  %  (00)     0.0                     0.0
  %  (10)   0.0024659               0.00247
  %  (20)   0.0074900               0.00739
  %  (30)   0.0147831               0.01478
  %  (40)   0.0247084               0.02462
  %  (50)   0.0369016               0.03689

  % Dowling J M 1968 
  %  The Rotation-Inversion Spectrum of Ammonia 
  %  J. Mol. Spectrosc. 27 527
  %
  % Al Derzi A R, Furtenbacher T, Tennyson J, Yurchenko S N and Császár A G 2015 
  %  MARVEL analysis of the measured high-resolution spectra of 14NH3 
  %  J. Quant. Spectrosc. Radiat. Transf. 161 117–30
  %
  % Furtenbacher T, Coles P A, Tennyson J, Yurchenko S N, Yu S, Drouin B, Tóbiás R and Császár A G 2020 
  %  Empirical rovibrational energy levels of ammonia up to 7500 cm−1 
  %  J. Quant. Spectrosc. Radiat. Transf. 251 107027–42
  
  if ~strcmp(state.type, 'rot')
    error(['Trying to asign rigid rotor energy to non rotational state %s. Check input file', state.name]);
  elseif isempty(state.gas.rotationalConstant)
    error(['Unable to find rotationalConstant to evaluate the energy of the state %s with function ' ...
      '''rigidRotorEnergy''.\nCheck input file'], state.name);
  end

  t = 6.626070040e-34*299792458/1.6021766208e-19;  % h*c/e
  J = str2double(state.rotLevel);
  energy = (state.gas.rotationalConstant*J*(J+1))-(t*(8.407*10^(-4))*100*(J*(J+1))^2) + (t*(2.38*10^(-7))*100*(J*(J+1))^3);

end
