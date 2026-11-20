function cases = l1_cases()
% L1_CASES  Analytic test cases of the L1 harness (variable-energy-grid project, task T5, D-011).
%
% Every case has an exact Maxwellian solution f(u) = C exp(-u/T) of the continuous two-term equation
% (zero energy flux), so the numerical EEDF can be compared with it on any grid:
%
%   constNu_DC       H2 mass, elastic sigma = c/sqrt(u) (constant collision frequency), E/N = 5 Td
%   constNu_HF       same, E/N = 40 Td at 2.45 GHz
%   maxwell_elastic  H2, real elastic cross section, E/N = 0
%   maxwell_CAR      as above + continuous approximation for rotations (CAR)
%   maxwell_ee       as maxwell_elastic + electron-electron collisions
%
% With sigma = c/sqrt(u) (u in eV), u*sigma^2 = c^2 is constant and the field and elastic terms give
%   T = kTg + (E/N)^2 / (6 (m/M) c^2 H),   H = 1 + (omega/N)^2 / (gamma^2 c^2),   gamma = sqrt(2e/me)
% (H = 1 for DC). At E/N = 0 the elastic, CAR and e-e operators all have the Maxwellian at kTg as
% equilibrium, whatever the cross section, so real H2 data are used (D-011).
%
% Fields of each case (compatible with l0_setup_text / l0_run_case):
%   name, gasBlock, ionType, growth, ee, algorithm, excFreq, reducedField  as in l0_cases
%   grid          struct(maxEnergy, cellNumber[, firstEnergyStep]) set per run by l1_run
%   maxEnergy     maximum energy of the case (eV), chosen for f(u_max)/f(0) ~ 1e-17
%   analytic      struct(type, EoNTd, excFreq, c) used by l1_check to evaluate T

  dataDir = '../Output/L1/data';   % LXCat files written by l1_data (paths relative to Code/Input)

  gasProps = {
    '  gasProperties:'
    '    mass: Databases/masses.txt'
    '    fraction:'
    '      - H2 = 1'
    '    harmonicFrequency: Databases/harmonicFrequencies.txt'
    '    anharmonicFrequency: Databases/anharmonicFrequencies.txt'
    '    rotationalConstant: Databases/rotationalConstants.txt'
    '    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt'
    '    OPBParameter: Databases/OPBParameter.txt'
    '  stateProperties:'
    '    population:'
    '      - H2(X) = 1.0'
    };
  constNu = [{'  LXCatFiles:'; ['    - ' dataDir '/H2_constant_nu.txt']}; gasProps];
  elastic = [{'  LXCatFiles:'; ['    - ' dataDir '/H2_elastic.txt']}; gasProps];
  elasticCAR = [{'  LXCatFiles:'; ['    - ' dataDir '/H2_elastic.txt']; '  CARgases:'; '    - H2'}; gasProps];

  c = l1_data('constantNu');   % sigma*sqrt(u) of the constant-nu cross section (m^2 eV^0.5)
  hf = 2.45e9;
  alg = 'mixingDirectSolutions';

  def = {
    % name             gasBlock    ee     excFreq EoN   maxEnergy
    'constNu_DC'       constNu     false  0       5     60
    'constNu_HF'       constNu     false  hf      40    60
    'maxwell_elastic'  elastic     false  0       0     1
    'maxwell_CAR'      elasticCAR  false  0       0     1
    'maxwell_ee'       elastic     true   0       0     1
    };

  for k = size(def, 1):-1:1
    cases(k).name = def{k, 1};
    cases(k).gasBlock = def{k, 2};
    cases(k).ionType = 'conservative';
    cases(k).growth = 'temporal';
    cases(k).ee = def{k, 3};
    cases(k).algorithm = alg;
    cases(k).excFreq = def{k, 4};
    cases(k).reducedField = num2str(def{k, 5});
    cases(k).grid = [];
    cases(k).maxEnergy = def{k, 6};
    if startsWith(def{k, 1}, 'constNu')
      cases(k).analytic = struct('type', 'constantNu', 'EoNTd', def{k, 5}, 'excFreq', def{k, 4}, 'c', c);
    else
      cases(k).analytic = struct('type', 'gasTemperature', 'EoNTd', 0, 'excFreq', 0, 'c', NaN);
    end
  end

end
