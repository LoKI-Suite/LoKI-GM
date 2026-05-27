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

function energy = morseOscillatorEnergy(state, ~, ~)
  % morseOscillatorEnergy is a property function that evaluates the energy of (vibrational) state j of a diatomic molecule, 
  %  adopting the Morse oscillator model
  
  if ~strcmp(state.type, 'vib')
    error('Trying to asign morse oscillator energy to non vibrational state %s. Check input file', state.name);
  elseif isempty(state.gas.harmonicFrequency)
    error(['Unable to find harmonicFrequency to evaluate the energy of the state %s with function ' ...
      '''harmonicOscillatorEnergy''.\nCheck input file'], state.name);
  elseif isempty(state.gas.anharmonicFrequency)
    error(['Unable to find anharmonicFrequency to evaluate the energy of the state %s with function ' ...
      '''morseOscillatorEnergy''.\nCheck input file'], state.name);
  end
  
  v = str2double(state.vibLevel);
  energy = Constant.planckReducedInEV*(state.gas.harmonicFrequency*(v+0.5)-state.gas.anharmonicFrequency*(v+0.5)^2);
  
end
