function statisticalWeight = rotationalDegeneracy_NH3(state, ~, ~)
  % rotationalDegeneracy_NH3 is a property function that evaluates the statistical weight of the intrinsic partition function 
  %  for a rotational state (J,K) of NH3, with K=0

  % The general degeneracy of rotational states is given by g_ns*(2J+1), with g_ns the nuclear spin degeneracy 
  % NH3 appears in ortho and para configurations, depending if K<>3n or K=3n, with n=1,2,3,...
  % 
  % In principle, the (J,0) rotational states should all belong to the ortho configuration, 
  %  therefore characterized by the same g_ns
  % 
  % However, the MARVEL database identifies a number of rotational states (J,K=0) 
  %  with the vibrational ground-state (v1=0, v2=0, v3=0, v4=0)
  %  consistent with the following degeneracies
  %  J = 0          2*(2J+1) = 2
  %  J = 1,2,3,...  (2J+1)
  % The following implements these expressions

  % Al Derzi A R, Furtenbacher T, Tennyson J, Yurchenko S N and Császár A G 2015 
  %  MARVEL analysis of the measured high-resolution spectra of 14NH3 
  %  J. Quant. Spectrosc. Radiat. Transf. 161 117–30
  
  if ~strcmp(state.type, 'rot')
    error(['Trying to asign rotational degeneracy to non rotational state %s. Check input file', state.name]);
  end
  
  J = str2double(state.rotLevel);
  if (J == 0)
    statisticalWeight = 2*(2*J+1);
  else
    statisticalWeight = 2*J+1;
  end    
  
end
