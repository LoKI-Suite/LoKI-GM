function out = l1_ee_operator(action, varargin)
% L1_EE_OPERATOR  e-e operator applied to a discrete Maxwellian (L1, D-011; review question T1.4).
%
%   T = l1_ee_operator('study')                   table over variants, grids, T and N
%   g = l1_ee_operator('guard', resFile, variant) check of the copied assembly against a LoKI run
%
% The e-e operator of LoKI is assembled inside the private method Boltzmann.solveEEColl, which always
% iterates from LoKI's own initial EEDF, so it cannot be applied to a chosen EEDF without changing
% Boltzmann.m. This file therefore holds line-by-line copies of the assembly of the auxiliary matrices
% (without the constant alpha) of
%   'branch'  branch 14, Boltzmann.m solveEEColl at 9867ca6 (variable grid, Maxwellian weights)
%   'master'  master, Boltzmann.m solveEEColl at fd96b57 (uniform grid only)
% and 'guard' checks each copy against the code it was taken from: applied to the converged EEDF of an
% L1 maxwell_ee run, the copy must reproduce, up to the constant alpha, the up/downflux vectors stored by
% LoKI (Afinal, Bfinal, read back from the snapshot's eeMatrix). The copies are to be replaced by a call
% to the code once T1.4/P3.1 exposes the e-e operator as a method.
%
% With Ahat = auxA*f, Bhat = auxB*f (both codes, alpha omitted) the operator is
%   (Mee f)_n = (J_(n-1/2) - J_(n+1/2))/du_n,   J_(n+1/2) = Ahat_n f_n - Bhat_(n+1) f_(n+1),
% so particles are conserved by construction. Quantities of the study:
%   detailedBalance  max_n |J_(n+1/2)| / max_n Ahat_n f_n for f_n = C exp(-u_n/T): zero if the discrete
%                    Maxwellian is an exact equilibrium of the operator
%   detailedBalanceTe  diagnostic: the same with Te forced to T in the Maxwellian weights of branch 14
%                    (the code computes Te = (2/3)<u> from the discrete EEDF, which is not exactly T);
%                    equal to detailedBalance for master, whose weights do not depend on Te
%   TeOverT          Te computed by the code from the discrete Maxwellian at T, divided by T
%   particles, energy  |sum du (Mee f)|/sum du |Mee f|, |sum du u (Mee f)|/sum du u |Mee f| for a
%                    non-Maxwellian f = exp(-(u/T)^2)

  switch action
    case 'study'
      out = study();
    case 'guard'
      out = guard(varargin{1}, varargin{2});
    otherwise
      error('l1_ee_operator:action', 'Unknown action "%s".', action);
  end

end

% ---------------------------------------------------------------------------------------------------
function T = study()
  kTg = 300*1.380649e-23/1.602176634e-19;
  temps = [kTg 1];
  cellNumbers = [100 200 400 800 1600];
  rows = {};
  for variant = {'branch', 'master'}
    for gridType = {'uniform', 'geometric'}
      if strcmp(variant{1}, 'master') && strcmp(gridType{1}, 'geometric'), continue; end
      for Tm = temps
        for N = cellNumbers
          [node, cell, du] = makeGrid(gridType{1}, 40*Tm, N);
          r.variant = variant{1}; r.grid = gridType{1}; r.T = Tm; r.N = N;
          f = maxwellian(cell, du, Tm);
          [auxA, auxB, Te] = assemble(variant{1}, node, cell, du, f);
          r.TeOverT = Te/Tm;
          r.detailedBalance = detailedBalance(auxA, auxB, f);
          [auxA, auxB] = assemble(variant{1}, node, cell, du, f, Tm);
          r.detailedBalanceTe = detailedBalance(auxA, auxB, f);
          g = exp(-(cell/Tm).^2); g = g/sum(du.*sqrt(cell).*g);
          [auxA, auxB] = assemble(variant{1}, node, cell, du, g);
          Mg = eeMatrix(auxA, auxB, g, du)*g(:);
          r.particles = abs(sum(du(:).*Mg))/sum(du(:).*abs(Mg));
          r.energy = abs(sum(du(:).*cell(:).*Mg))/sum(du(:).*cell(:).*abs(Mg));
          rows{end+1} = r; %#ok<AGROW>
        end
      end
    end
  end
  T = struct2table([rows{:}]);
  T.pDetailedBalance = NaN(height(T), 1);
  for k = 2:height(T)
    if T.N(k) == 2*T.N(k-1) && strcmp(T.variant{k}, T.variant{k-1}) && strcmp(T.grid{k}, T.grid{k-1})
      T.pDetailedBalance(k) = log2(T.detailedBalance(k-1)/T.detailedBalance(k));
    end
  end
end

function g = guard(resFile, variant)
% max relative spread of (LoKI flux vector)./(copied flux vector) over cells, after removing alpha
  S = load(resFile);
  s = S.res.jobs{end};
  M = s.M.eeMatrix;
  du = s.du(:)'; cell = s.cell(:)'; node = s.node(:)'; f = s.eedf(:)';
  N = numel(du);
  AhatCode = diag(M, -1)'.*du(2:N);       % Mee(n+1,n) = Ahat_n/du_(n+1)
  BhatCode = diag(M, 1)'.*du(1:N-1);      % Mee(n,n+1) = Bhat_(n+1)/du_n
  [auxA, auxB] = assemble(variant, node, cell, du, f);
  Ahat = (auxA*f')'; Bhat = (auxB*f')';
  rA = AhatCode./Ahat(1:N-1);
  rB = BhatCode./Bhat(2:N);
  keep = isfinite(rA) & Ahat(1:N-1) > 1e-300 & isfinite(rB) & Bhat(2:N) > 1e-300;
  r = [rA(keep) rB(keep)];
  g = max(abs(r/median(r) - 1));
end

% ---------------------------------------------------------------------------------------------------
function [node, cell, du] = makeGrid(type, uMax, N)
  switch type
    case 'uniform'
      node = linspace(0, uMax, N + 1);
    case 'geometric'                     % first step uMax/(10N), as in l1_run
      u1 = uMax/(10*N);
      % u1*(a^N - 1)/(a - 1) = uMax, solved for x = a - 1 in log form (a^N overflows for a = 2, N = 1600)
      x = fzero(@(x) log(expm1(N*log1p(x))) - log(x) - log(uMax/u1), [1e-9, 5/N]);
      a = 1 + x;
      node = [0 cumsum(u1*a.^(0:N-1))];
      node(end) = uMax;
  end
  du = diff(node);
  cell = node(1:end-1) + du/2;
end

function f = maxwellian(cell, du, T)
  f = exp(-cell/T);
  f = f/sum(du.*sqrt(cell).*f);
end

function d = detailedBalance(auxA, auxB, f)
  Ahat = (auxA*f')'; Bhat = (auxB*f')';
  up = Ahat(1:end-1).*f(1:end-1);
  J = up - Bhat(2:end).*f(2:end);
  d = max(abs(J))/max(up);
end

function Mee = eeMatrix(auxA, auxB, f, du)
  A = (auxA*f(:))'; B = (auxB*f(:))'; du = du(:)';
  n = numel(du);
  Mee = diag(-(A + B)./du) + diag(A(1:n-1)./du(2:n), -1) + diag(B(2:n)./du(1:n-1), 1);
end

% ---------------------------------------------------------------------------------------------------
function [auxA, auxB, Te] = assemble(variant, energyNode, energyCell, energyStep, eedf, TeForced)
% Copies of the e-e assembly (auxiliary matrices without alpha). Row vectors, as in LoKI.
% TeForced (optional, diagnostic only) replaces the Te computed from eedf in branch 14's weights.
  cellNumber = numel(energyCell);
  switch variant
    case 'branch'   % Boltzmann.m solveEEColl, branch 14 at 9867ca6 (lines 1722-1767)
      weightGrid(1:cellNumber-1) = energyStep(1:cellNumber-1)./...
          (energyStep(2:cellNumber)+energyStep(1:cellNumber-1));
      auxA = zeros(cellNumber);
      auxB = zeros(cellNumber);
      auxMatrix = zeros(cellNumber);
      for k=2:cellNumber
        auxMatrix(k,1:k-1) = energyStep(1:k-1).*((weightGrid(k-1)).*sqrt(energyCell(1:k-1))+...
            (2/3)*2*weightGrid(k-1)./energyStep(k-1).*energyCell(1:k-1).^(3/2));
        auxMatrix(k,k:cellNumber) = (2/3)*energyStep(k:cellNumber).*...
            2*weightGrid(k-1)./energyStep(k-1).*energyNode(k).^(3/2);
      end
      meanEnergy = dot(energyStep.*energyCell.^(3/2),eedf);
      Te = (2/3)*meanEnergy;
      if nargin > 5, Te = TeForced; end
      maxwellianExponent = zeros(cellNumber);
      maxwellianFactor = zeros(cellNumber);
      for k=1:cellNumber-1
        maxwellianExponent(2:cellNumber,k) = 1/Te*0.5*(energyStep(1:cellNumber-1)./...
            weightGrid(1:cellNumber-1)-energyStep(k)./weightGrid(k));
        maxwellianFactor(2:cellNumber,k) = weightGrid(1:cellNumber-1)./weightGrid(k).*...
           energyStep(k)./energyStep(1:cellNumber-1).*exp(maxwellianExponent(2:cellNumber,k)');
        auxB(2:cellNumber,k) = sqrt(maxwellianFactor(2:cellNumber,k).*auxMatrix(2:cellNumber,k).*...
            auxMatrix(k+1,1:cellNumber-1)');
      end
      for k=2:cellNumber
          for j=1:cellNumber-1
            auxA(j,k) = weightGrid(j)/weightGrid(k-1)*energyStep(k-1)/energyStep(j)*auxB(k,j);
          end
      end

    case 'master'   % Boltzmann.m solveEEColl, master at fd96b57 (uniform step); A = (alpha/du)*auxA*f,
                    % Mee not divided by du, i.e. Ahat = du*A/alpha = auxA*f as for branch 14
      if max(abs(energyStep - energyStep(1))) > 1e-12*energyStep(1)
        error('l1_ee_operator:master', 'The master e-e operator is defined on uniform grids only.');
      end
      step = energyStep(1);
      auxA = zeros(cellNumber);
      auxMatrix = zeros(cellNumber);
      auxEnergyArray = -(step/2)*sqrt(energyCell)+(2/3)*energyCell.^(3/2);
      for k=1:cellNumber-1
        auxMatrix(k,2:k) = auxEnergyArray(2:k);
        auxMatrix(k,k+1:cellNumber) = (2/3)*energyNode(k+1)^(3/2);
      end
      for k=1:cellNumber-1
        auxA(k,2:cellNumber) = sqrt(auxMatrix(k,2:cellNumber).*auxMatrix(1:cellNumber-1,k+1)');
      end
      auxB = transpose(auxA);
      meanEnergy = step*dot(energyCell.^(3/2),eedf);
      Te = (2/3)*meanEnergy;
  end
end
