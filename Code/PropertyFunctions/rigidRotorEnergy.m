function energy = rigidRotorEnergy(state, ~, ~)
  % rigidRotorEnergy is a property function that evaluates the energy of rotational state J of a diatomic molecule, 
  %  adopting the rigid-rotor model
  
  if ~strcmp(state.type, 'rot')
    error(['Trying to asign rigid rotor energy to non rotational state %s. Check input file', state.name]);
  elseif isempty(state.gas.rotationalConstant)
    error(['Unable to find rotationalConstant to evaluate the energy of the state %s with function ' ...
      '''rigidRotorEnergy''.\nCheck input file'], state.name);
  end
  
  J = str2double(state.rotLevel);
  energy = state.gas.rotationalConstant*J*(J+1);
  
end
