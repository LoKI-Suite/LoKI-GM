% LoKI-B solves a time and space independent form of the two-term 
% electron Boltzmann equation (EBE), for non-magnetised non-equilibrium 
% low-temperature plasmas excited by DC/HF electric fields from 
% different gases or gas mixtures.
% Copyright (C) 2018 A. Tejero-del-Caz, V. Guerra, D. Goncalves, 
% M. Lino da Silva, L. Marques, N. Pinhao, C. D. Pintassilgo and
% L. L. Alves
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
