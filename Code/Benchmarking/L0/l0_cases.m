function cases = l0_cases()
% L0_CASES  Test configurations of the L0 harness (variable-energy-grid project, task T5).
%
% Each case is a LoKI-B stationary configuration run at E/N = 0.01, 0.1, 1, 10, 100, 1000 Td.
% N2 is run for every valid combination of
%   ionization operator type  conservative, oneTakesAll, equalSharing, usingSDCS
%   growth model              temporal, spatial
%   e-e collisions            off, on
%   field                     DC, HF (2.45 GHz; HF is not allowed with the spatial growth model)
% i.e. 4 x 2 x 2 (DC) + 4 x 1 x 2 (HF) = 24 cases. N2 includes CAR and vibrational superelastic
% collisions. Ar (atomic gas) and O2 (attachment, statistical weights) are added as one case each.
% The temporalIntegration algorithm is left to L4 (D-006).
%
% Fields of each case:
%   name          case identifier (used for file names)
%   gasBlock      setup-file lines from 'LXCatFiles:' to the end of 'stateProperties:'
%   ionType       ionizationOperatorType
%   growth        growthModelType
%   ee            includeEECollisions (logical)
%   algorithm     nonLinearRoutines algorithm
%   excFreq       excitation frequency (Hz)
%   reducedField  reduced-field list (Td), as setup-file text
%   grid          [] (smart grid, 200 cells, as in L0) or struct with fields maxEnergy and
%                 cellNumber for a fixed grid without smart grid (used by l0_refine)

  N2 = {
    '  LXCatFiles:'
    '    - Nitrogen/N2_LXCat.txt'
    '  CARgases:'
    '    - N2'
    '  gasProperties:'
    '    mass: Databases/masses.txt'
    '    fraction:'
    '      - N2 = 1'
    '    harmonicFrequency: Databases/harmonicFrequencies.txt'
    '    anharmonicFrequency: Databases/anharmonicFrequencies.txt'
    '    rotationalConstant: Databases/rotationalConstants.txt'
    '    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt'
    '    OPBParameter: Databases/OPBParameter.txt'
    '  stateProperties:'
    '    energy:'
    '      - N2(X,v=*) = harmonicOscillatorEnergy'
    '    statisticalWeight:'
    '      - N2(X,v=*) = 1.0'
    '    population:'
    '      - N2(X) = 1.0'
    '      - Nitrogen/N2_vibpop.txt'
    };

  Ar = {
    '  LXCatFiles:'
    '    - Argon/Ar_LXCat.txt'
    '  gasProperties:'
    '    mass: Databases/masses.txt'
    '    fraction:'
    '      - Ar = 1'
    '    harmonicFrequency: Databases/harmonicFrequencies.txt'
    '    anharmonicFrequency: Databases/anharmonicFrequencies.txt'
    '    rotationalConstant: Databases/rotationalConstants.txt'
    '    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt'
    '    OPBParameter: Databases/OPBParameter.txt'
    '  stateProperties:'
    '    population:'
    '      - Ar(1S0) = 1.0'
    };

  O2 = {
    '  LXCatFiles:'
    '    - Oxygen/O2_LXCat.txt'
    '  gasProperties:'
    '    mass: Databases/masses.txt'
    '    fraction:'
    '      - O2 = 1'
    '    harmonicFrequency: Databases/harmonicFrequencies.txt'
    '    anharmonicFrequency: Databases/anharmonicFrequencies.txt'
    '    rotationalConstant: Databases/rotationalConstants.txt'
    '    electricQuadrupoleMoment: Databases/quadrupoleMoment.txt'
    '    OPBParameter: Databases/OPBParameter.txt'
    '  stateProperties:'
    '    energy:'
    '      - O2(X,v=*) = harmonicOscillatorEnergy'
    '    statisticalWeight:'
    '      - O2(X) = 3'
    '      - O2(a1Dg) = 2'
    '      - O2(b1Sg+) = 1'
    '      - O2(A3Su+_C3Du_c1Su-) = 10'
    '      - O2(9.97eV)= 3'
    '      - O2(14.7eV) = 1'
    '      - O2(X,v=*) = 3.0'
    '    population:'
    '      - O2(X) = 1.0'
    '      - O2(X,v=0) = 1.0'
    };

  EoN = '[0.01 0.1 1 10 100 1000]';
  alg = 'mixingDirectSolutions';
  hf = 2.45e9;

  def = cell(0, 8);
  for ion = {'conservative', 'oneTakesAll', 'equalSharing', 'usingSDCS'}
    for growth = {'temporal', 'spatial'}
      for ee = [false true]
        for f = [0 hf]
          if f > 0 && strcmp(growth{1}, 'spatial'), continue; end   % not allowed in LoKI-B
          name = ['N2_' ion{1} '_' growth{1}];
          if ee, name = [name '_ee']; end %#ok<AGROW>
          if f > 0, name = [name '_HF']; end %#ok<AGROW>
          def(end+1, :) = {name, N2, ion{1}, growth{1}, ee, alg, f, EoN}; %#ok<AGROW>
        end
      end
    end
  end
  def(end+1, :) = {'Ar_usingSDCS_temporal', Ar, 'usingSDCS', 'temporal', false, alg, 0, EoN};
  def(end+1, :) = {'O2_usingSDCS_spatial_attachment', O2, 'usingSDCS', 'spatial', false, alg, 0, EoN};

  cases = cell2struct(def, {'name','gasBlock','ionType','growth','ee','algorithm','excFreq','reducedField'}, 2);
  [cases.grid] = deal([]);

end
