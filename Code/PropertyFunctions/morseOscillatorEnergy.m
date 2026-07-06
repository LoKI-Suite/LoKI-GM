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
