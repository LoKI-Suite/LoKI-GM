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

classdef Collision < handle
  %Collision Class that stores the information of a certain electron collision
  %   Class that stores the information of an electron collision read from
  %   an LXCat file. The information stored here (in particular the cross 
  %   section) is to be used in a Boltzmann solver to obtain the EEDF.

  properties

    ID = -1;    % ID that identifies the collision in the collision 
                %  array of the gas
    
		type = '';  % type of collision. Supported types defined in constant property supportedCollisions

    target = State.empty;       % handle to the target of the collision
    productArray = State.empty; % handle to the products of the collision
    productStoiCoeff = [];      % array of stoichiometric coefficients for the products
    
    isExtra = false;                      % true when the collision is not to be used in the electron kinetics calculations
    isReverse = false;                    % true when super elastic collision is defined
    threshold = 0.0;                      % energy threshold of the collision (eV)
    rawCrossSection = [];                 % integral cross section as read from LXCat file, 2 rows (eV m^2)
    rawMomentumTransferCrossSection = []; % momentum-transfer cross section as read from LXCat file, 2 rows (eV m^2)
    
    energyGrid = Grid.empty;            % handle to the energy grid of the simulation
    energyGridListener;                 % handle to the listener of changes in the energy grid
    crossSection = [];                  % integral cross section interpolated into the grid, 1 row (m^2)
    momentumTransferCrossSection = [];  % momentum-transfer cross section interpolated into the grid, 1 row (m^2)

    
    ineRateCoeff = [];  % inelastic rate coefficient of the collision obtained once the eedf is known (->)
    supRateCoeff = [];  % superelastic rate coefficient of the collision obtained once the eedf is known (if <->)
    
  end

  properties (Constant = true)
    supportedCollisions = {'Effective' 'Elastic' 'Rotational' 'Vibrational' 'Ionization' 'Attachment' 'Excitation'};
  end

  methods (Access = public)

    function collision = Collision(type, isMomentumTransfer, target, productArray, productStoiCoeff, isReverse, ...
        threshold, rawCrossSection, isExtra)
      persistent lastID;
      if isempty(lastID)
        lastID = 0;
      end
      lastID = lastID + 1;
      collision.ID = lastID;
      collision.target = target;
      collision.productArray = productArray;
      collision.productStoiCoeff = productStoiCoeff;
      collision.isExtra = isExtra;
      collision.isReverse = isReverse;
      collision.threshold = threshold;
      collision.type = type;
      if strcmp('Elastic', type) || strcmp('Effective', type) || ~isMomentumTransfer
        collision.rawCrossSection = rawCrossSection;
      else
        collision.rawMomentumTransferCrossSection = rawCrossSection;
      end
      if ~any(strcmp(type, collision.supportedCollisions))
        strSupportedCollisions = '';
        for idx = 1:length(collision.supportedCollisions)-1
          strSupportedCollisions = [strSupportedCollisions collision.supportedCollisions{idx} ', '];
        end
        strSupportedCollisions = [strSupportedCollisions(1:end-2) ' and ' collision.supportedCollisions{end}];
        error(['Error while creating collision ''%s''.\n''%s'' collisions are not supported by LoKI.\nPlease, ' ...
          'check your LXCat files.\nSupported collisions are: %s\n'], collision.description, type, ...
          strSupportedCollisions);
      end
      if isExtra
        target.gas.collisionArrayExtra(end+1) = collision;
        target.collisionArrayExtra(end+1) = collision;
      else
        target.gas.collisionArray(end+1) = collision;
        target.collisionArray(end+1) = collision;
        target.isTarget = true;
      end
      if isReverse
        if length(productArray)==1 && productStoiCoeff==1
          if isExtra
            productArray.collisionArrayExtra(end+1) = collision;
          else
            productArray.collisionArray(end+1) = collision;
            productArray.isTarget = true;
          end
        else
          error(['Error while creating collision ''%s''.\nKlein-Rosseland microreversibility relation valid only ' ...
            'for binary collisions'], collision.description);
        end
      end
      if strcmp(type, 'Effective') && ~strcmp(target.type, 'ele')
        error(['Found ''%s'' collision with ''%s'' as target [%s].\n%s collisions are only allowed for ' ...
          'electronic states, please check LXCat files.\n'], type, target.name, collision.description, type);
      end
    end
    
    function disp(collision)
      fprintf('ID: %d\n', collision.ID);
      fprintf('description: %s\n', collision.description);
      fprintf('threshold: %f eV\n', collision.threshold);
      if isempty(collision.energyGrid) 
        if collision.isReverse
          loglog(collision.rawCrossSection(1,:), collision.rawCrossSection(2,:), '-k', ...
            collision.rawCrossSection(1,:), collision.superElasticCrossSection, '-b');
          legend([collision.type ' (Raw)'], 'Superelastic (Raw)');
        else
          loglog(collision.rawCrossSection(1,:), ...
            collision.rawCrossSection(2,:), '-k');
          legend([collision.type ' (Raw)']);
        end
      else
        if collision.isReverse
          loglog(collision.rawCrossSection(1,:), collision.rawCrossSection(2,:), '-k', ...
            collision.energyGrid.node, collision.crossSection, 'or', ...
            collision.energyGrid.node, collision.superElasticCrossSection, 'xb');
          legend([collision.type ' (Raw)'], [collision.type ' (Interpolated)'], 'Superelastic (Interpolated)');
        else
          loglog(collision.rawCrossSection(1,:), collision.rawCrossSection(2,:), '-k', ...
            collision.energyGrid.node, collision.crossSection, 'or');
          legend([collision.type ' (Raw)'], [collision.type ' (Interpolated)']);
        end
      end
      title(collision.description);
      
    end
    
    function collisionStr = description(collision)
      collisionStr = ['e+' collision.target.name];
      if collision.isReverse
        collisionStr = [collisionStr '<->'];
      else
        collisionStr = [collisionStr '->'];
      end
      if strcmp(collision.type, 'Ionization')
        collisionStr = [collisionStr 'e+e+'];
      elseif ~strcmp(collision.type, 'Attachment')
        collisionStr = [collisionStr 'e+'];
      end
      numberOfProducts = length(collision.productArray);
      for i = 1:numberOfProducts
        product = collision.productArray(i);
        if collision.productStoiCoeff(i) > 1
          stoiCoeffStr = sprintf('%d',collision.productStoiCoeff(i));
        else
          stoiCoeffStr = '';
        end
        collisionStr = [collisionStr stoiCoeffStr product.name];
        if i < numberOfProducts
          collisionStr = [collisionStr '+'];
        end
      end
      collisionStr = [collisionStr ',' collision.type];
    end
    
    function adjustCrossSection(collision, energyGrid)
    % adjustCrossSection interpolates the values of the cross section of a
    % certain collision and stores the adjusted values in the 'crossSection'
    % property of the object 'collision'. The function also stores the
    % energy grid used for the interpolation in the 'energyGrid' property.
    
      % check that 'Elastic' cross section is available for the maximum value of the energy grid
      if strcmp(collision.type, 'Elastic')
        if energyGrid.node(end) > collision.rawCrossSection(1,end)
          error(['''%s'' cross section data is not available for the maximum energy of the simulation (%f eV).\n'...
            'Simulation is not reliable under this conditions.'], collision.description, energyGrid.node(end));
        end
      end
      
      % save energy grid and add the corresponding listener to update the interpolation acording to the grid
      if isempty(collision.energyGrid)
        collision.energyGrid = energyGrid;
        collision.energyGridListener = addlistener(energyGrid, 'updatedMaxEnergy1', @collision.reAdjustCrossSection);
      else
        collision.energyGrid = energyGrid;
        delete(collision.energyGridListener);
        collision.energyGridListener = addlistener(energyGrid, 'updatedMaxEnergy1', @collision.reAdjustCrossSection);
      end
      
      % save interpolated values of the cross section in crossSection property
      [collision.crossSection, collision.momentumTransferCrossSection] = ...
        collision.interpolatedCrossSection(energyGrid.node);
      
    end
        
    function [crossSection, momentumTransferCrossSection] = interpolatedCrossSection(collision, energyValues)
    % interpolatedCrossSection returns the values of the cross section of a
    % "collision" interpolated at certain "energyValues". The interpolation
    % is performed with the interp1 matlab function, with linear
    % interpolation and an extrapolation value of 0.0. 
    
      crossSection = zeros(1,length(energyValues));
      minIndex = 0;
      if strcmp(collision.type, 'Effective') || strcmp(collision.type, 'Elastic')
        minIndex = 1;
      else
        for i = 1:length(energyValues)
          if energyValues(i)>collision.threshold
            minIndex = i;
            break;
          end
        end
      end
      if minIndex > 0
        crossSection(minIndex:end) = interp1(collision.rawCrossSection(1,:), collision.rawCrossSection(2,:), ...
          energyValues(minIndex:end), 'linear', 0.0);
      end

      if ~isempty(collision.rawMomentumTransferCrossSection)
        momentumTransferCrossSection = zeros(1,length(energyValues));
        minIndex = 0;
        if strcmp(collision.type, 'Effective') || strcmp(collision.type, 'Elastic')
          minIndex = 1;
        else
          for i = 1:length(energyValues)
            if energyValues(i)>collision.threshold
              minIndex = i;
              break;
            end
          end
        end
        if minIndex > 0
          momentumTransferCrossSection(minIndex:end) = interp1(collision.rawMomentumTransferCrossSection(1,:), ...
            collision.rawMomentumTransferCrossSection(2,:), energyValues(minIndex:end), 'linear', 0.0);
        end
      else
        momentumTransferCrossSection = [];
      end
      
    end
    
    function [crossSection, momentumTransferCrossSection] = superElasticCrossSection(collision, energyValue)
    % superElasticCrossSection returns the super elastic cross section of a certain collision at the specified 
    % energyValues by using the Klein-Rosseland michroreversivility relation:
    % 
    %    sigma_{f,i}(u) = (g_i/g_f)*(1+threshold/u)*sigma_{i,f}(u+threshold)
    %
    % Note: in case energyValue is not provided the function uses the energy values of the collision raw cross section.
    
      % error checking 
      if ~collision.isReverse
        error('Collision ''%s'' is not defined as bidirectional.', collision.description);
      elseif isempty(collision.target.statisticalWeight)
        error(['The statistical weight of the state ''%s'' is not defined.\nThe super elastic cross section of ' ...
          'the collision ''%s'' can not be evaluated.'], collision.target.name, collision.description);
      elseif isempty(collision.productArray.statisticalWeight)
        error(['The statistical weight of the state ''%s'' is not defined.\nThe super elastic cross section of ' ...
          'the collision ''%s'' can not be evaluated.'], collision.productArray.name, collision.description);
      end
      
      % when no energyValue is given, the function returns the super elastic
      % cross section at: energyGrid.node (if available) or rawCrossSection
      % energy values
      if nargin == 1
        if isempty(collision.energyGrid)
          energyValue = collision.rawCrossSection(1,:);
        else
          energyValue = collision.energyGrid.node;
        end
      end
      
      % Klein-Rosseland microreversivility relation (super elastic cross
      %  section is zero at zero energy)
      if energyValue(1)==0
        if length(energyValue) == 1
          return;
        else
          minIndex = 2;
        end
      else
        minIndex = 1;
      end
      aux = (collision.target.statisticalWeight/collision.productArray.statisticalWeight)*...
        (1+collision.threshold./energyValue(minIndex:end));
      [directCrossSection, directMomentumTransferCrossSection] = ...
        collision.interpolatedCrossSection(energyValue(minIndex:end)+collision.threshold);
      crossSection = zeros(1,length(energyValue));
      crossSection(minIndex:end) = aux.*directCrossSection;
      if isempty(directMomentumTransferCrossSection)
        momentumTransferCrossSection = [];
      else
        momentumTransferCrossSection = zeros(1,length(energyValue));
        momentumTransferCrossSection(minIndex:end) = aux.*directMomentumTransferCrossSection;
      end
      
    end
    
    function [ineRateCoeff, supRateCoeff] = evaluateRateCoeff(collision, eedf)
    % evaluateRateCoeff evaluates the rate coefficient(s) of a certain collision provided an eedf. The function
    % returns the value(s) of the rate coefficient for the inelastic (and superelastic) collision. The value(s) 
    % is (are) also stored in the collision properties.
    
      % evaluate auxiliary variables
      gamma = Constant.gamma;
      lmin = floor(collision.threshold/collision.energyGrid.step);
      cellCrossSection = (collision.crossSection(lmin+1:end-1)+collision.crossSection(lmin+2:end))./2;
      aux = cellCrossSection.*collision.energyGrid.cell(lmin+1:end);
    
      % evaluate inelastic rate coefficient
      ineRateCoeff = gamma*sum(aux.*eedf(lmin+1:end))*collision.energyGrid.step;
      collision.ineRateCoeff = ineRateCoeff;
      
      % evaluate superelastic rate coefficient (if collision is reverse)
      if collision.isReverse
        targetStatWeight = collision.target.statisticalWeight;
        productStatWeight = collision.productArray.statisticalWeight;
        if isempty(targetStatWeight)
          error(['Unable to find %s statistical weight for the evaluation of superelastic rate coefficient of ' ...
            '''%s''.\nCheck input file'], collision.target.name, collision.description);
        elseif isempty(productStatWeight)
          error(['Unable to find %s statistical weight for the evaluation of superelastic rate coefficient of ' ...
            '''%s''.\nCheck input file'], collision.productArray.name, collision.description);
        else
          statWeightRatio = targetStatWeight/productStatWeight;
        end
        supRateCoeff = gamma*statWeightRatio*sum(aux.*eedf(1:end-lmin))*collision.energyGrid.step;
        collision.supRateCoeff = supRateCoeff;
      else
        supRateCoeff = [];
      end
      
    end
    
  end
  
  methods (Access = private)
    
    function reAdjustCrossSection(collision, energyGrid, ~)
    % readjustCrossSection readjust the interpolated cross section of the collision whenever the 'updatedMaxEnergy1'
    % event is trigered by the energyGrid object
      
      % check that 'Elastic' cross section is available for the maximum value of the energy grid
      if strcmp(collision.type, 'Elastic')
        if energyGrid.node(end) > collision.rawCrossSection(1,end)
          error(['''%s'' cross section data is not available for the maximum energy of the simulation (%f eV).\n'...
            'Simulation is not reliable under this conditions.'], collision.description, energyGrid.node(end));
        end
      end
      
      % save new interpolated values of the cross section in crossSection property
      [collision.crossSection, collision.momentumTransferCrossSection] = ...
        collision.interpolatedCrossSection(energyGrid.node);
      
    end
    
  end
  
  methods (Static)
    
    function [collisionArray, collisionID] = add(type, isMomentumTransfer, target, productArray, productStoiCoeff, ...
        isReverse, threshold, rawCrossSection, collisionArray, isExtra)
    % addCollision (needs to be written)
      
      collisionID = Collision.find(type, target, productArray, productStoiCoeff, isReverse, threshold);
      if collisionID == -1
        collisionArray(end+1) = Collision(type, isMomentumTransfer, target, productArray, productStoiCoeff, ...
          isReverse, threshold, rawCrossSection, isExtra);
        collisionID = collisionArray(end).ID;
      elseif isMomentumTransfer && isempty(collisionArray(collisionID).rawMomentumTransferCrossSection)
        collisionArray(collisionID).rawMomentumTransferCrossSection = rawCrossSection;
      elseif ~isMomentumTransfer && isempty(collisionArray(collisionID).rawCrossSection)
        collisionArray(collisionID).rawCrossSection = rawCrossSection;
      else
        warning('Avoiding duplicated electron impact collision:\n\t''%s''\n', collisionArray(collisionID).description);
      end
      
    end
    
    function collisionID = find(type, target, productArray, productStoiCoeff, isReverse, threshold)
    % find (needs to be written)
      
      for collision = [ target.collisionArray target.collisionArrayExtra ]
        if collision.threshold == threshold && strcmp(collision.type, type) && collision.isReverse == isReverse
          numProducts = length(collision.productArray);
          if numProducts == length(productArray)
            equalProducts = true;
            for j = 1:numProducts
              for k = 1:numProducts
                if collision.productArray(j) == productArray(k) && collision.productStoiCoeff(j) == productStoiCoeff(k)
                  break;
                elseif k == numProducts
                  equalProducts = false;
                end
              end
              if ~equalProducts
                break;
              end
            end
            if equalProducts
              collisionID = collision.ID;
              return;
            end
          end
        end
      end
      collisionID = -1;
      
    end
    
    function [collisionID, eqColl] = findEquivalent(target, productArray, productStoiCoeff, isReverse)
    % findEquivalent finds an electron collision from the electron kinetics collisionArray with a prescribed 
    % directionality, target, productArray and productStoiCoeff.
      
      for collision = [ target.collisionArrayExtra target.collisionArray ]
        if collision.isReverse == isReverse
          numProducts = length(collision.productArray);
          if numProducts == length(productArray)
            equalProducts = true;
            for j = 1:numProducts
              for k = 1:numProducts
                if collision.productArray(j) == productArray(k) && collision.productStoiCoeff(j) == productStoiCoeff(k)
                  break;
                elseif k == numProducts
                  equalProducts = false;
                end
              end
              if ~equalProducts
                break;
              end
            end
            if equalProducts
              collisionID = collision.ID;
              eqColl = collision;
              return;
            end
          end
        end
      end
      collisionID = -1;
      eqColl = Collision.empty;
      
    end
    
    function adjustToEnergyGrid(energyGrid, collisionArray)
    
      % adjustic cross sections of each collision to the energy grid
      for collision = collisionArray
        collision.adjustCrossSection(energyGrid);
      end
      
    end
    
  end

end
