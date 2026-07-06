function statisticalWeight = rotationalDegeneracy_H2(state, ~, ~)
  % rotationalDegeneracy_H2 is a property function that evaluates the statistical weight of the intrinsic partition function 
  %  of a rotational state J of the H2 molecule
  % Ridenti M A, Alves L L, Guerra V and Amorim J 2015 Plasma Sources Sci. Technol. 24 035002
  
  if ~strcmp(state.type, 'rot')
    error(['Trying to asign rotational degeneracy to non rotational state %s. Check input file', state.name]);
  end
  
  J = str2double(state.rotLevel);
  statisticalWeight = (2-(-1)^J)*(2*J+1);
  
end
