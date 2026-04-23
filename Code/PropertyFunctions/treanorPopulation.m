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

function population = treanorPopulation(state, argumentArray, workCond)
  % treanorPopulation is a property function that evaluates a Treanor distribution
  %  with certain temperatures (specified as the first two arguments provided to the function in the setup file)
  
  temp0 = argumentArray{1};
  temp1 = argumentArray{2};
  if ~isnumeric(temp0)
    switch temp0
      case 'gasTemperature'
        temp0 = workCond.gasTemperature;
      case 'electronTemperature'
        temp0 = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temp0);
    end
  end
  if ~isnumeric(temp1)
    switch temp1
      case 'gasTemperature'
        temp1 = workCond.gasTemperature;
      case 'electronTemperature'
        temp1 = workCond.electronTemperature/Constant.boltzmannInEV;
      otherwise
        error(['Error found when evaluating population of state %s.\nTemperature ''%s'' not defined in the ' ...
          'working conditions.\nPlease, fix the problem and run the code again.'], state.name, temp1);
    end
  end
  
  if ~strcmp(state.type, 'vib')
    error(['Trying to asign treanor population to non vibrational state %s. Check input file', state.name]);
  end
  
  E0 = [];
  E1 = [];
  for stateAux = [state state.siblingArray]
    if strcmp(stateAux.vibLevel, '0')
      E0 = stateAux.energy;
    elseif strcmp(stateAux.vibLevel, '1')
      E1 = stateAux.energy;
    end
  end
  if isempty(E0) || isempty(E1)
    error('Unable to find E0 or E1 to populate state %s and its siblings with function %s.\nCheck input file', ...
      state.name, 'treanorPopulation');
  end
  
  % evaluate Treanor
  norm = 0;
  for stateAux = [state state.siblingArray]
    v = str2double(stateAux.vibLevel);
    E = stateAux.energy;
    g = stateAux.statisticalWeight;
    if isempty(E)
      error('Unable to find %s energy for the evaluation of %s function.\nCheck input file', ...
        stateAux.name, 'treanorPopulation');
    elseif isempty(g)
      error('Unable to find %s statistical weight for the evaluation of %s function.\nCheck input file', ...
        stateAux.name, 'treanorPopulation');
    end
    stateAux.population = g*exp(-(v*(E1-E0)*(1/temp1-1/temp0)+(E-E0)/temp0)/Constant.boltzmannInEV);
    norm = norm + stateAux.population;
  end
  for stateAux = [state state.siblingArray]
    stateAux.population = stateAux.population/norm;
  end
  
  population = state.population;
  
end
