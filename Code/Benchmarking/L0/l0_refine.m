function out = l0_refine(action, varargin)
% L0_REFINE  Grid-refinement study of the e-e case: branch 14 vs reference code on uniform grids.
%
%   l0_refine('run', mode)                      run the study in the tree of the current folder
%                                               (mode 'ref' in the reference tree, 'uniform' in branch 14)
%   T = l0_refine('report', refFolder, file)    compare (run from the branch-14 Code folder) and write
%                                               a Markdown report
%
% Case N2_conservative_temporal_ee (l0_cases) at fixed grids without smart grid: for each reduced field a fixed
% maxEnergy, and cellNumber = 100, 200, ..., 3200. If both codes discretize the same equation
% consistently, the gap between them must decrease with the cell size (observed order p > 0).

  study.reducedField = [0.1 1 10 100];
  study.maxEnergy    = [2 10 15 60];        % eV, about the smart-grid values of the L0 run
  study.cellNumber   = [100 200 400 800 1600 3200];

  switch action
    case 'run'
      mode = varargin{1};
      base = l0_cases();
      base = base(strcmp({base.name}, 'N2_conservative_temporal_ee'));
      outDir = fullfile('Output', 'L0', ['refine_' mode]);
      for i = 1:numel(study.reducedField)
        for N = study.cellNumber
          c = base;
          c.name = caseName(study.reducedField(i), N);
          c.reducedField = num2str(study.reducedField(i));
          c.grid = struct('maxEnergy', study.maxEnergy(i), 'cellNumber', N);
          l0_run_case(c, mode, outDir);
        end
      end
      out = [];

    case 'report'
      refFolder = varargin{1};
      rows = {};
      for i = 1:numel(study.reducedField)
        EoN = study.reducedField(i);
        prev = [];
        for N = study.cellNumber
          s = load1(fullfile('Output', 'L0', 'refine_uniform', [caseName(EoN, N) '.mat']));
          q = load1(fullfile(refFolder, 'Output', 'L0', 'refine_ref', [caseName(EoN, N) '.mat']));
          r = struct('EoN', EoN, 'maxEnergy', study.maxEnergy(i), 'N', N, 'du', q.du(1), ...
            'eedfGap', NaN, 'meanEnergyGap', NaN, 'eeGainGap', NaN, 'eeMatrixGap', NaN, ...
            'pEedf', NaN, 'pMeanEnergy', NaN, 'pEeGain', NaN, 'eeNetRefBranch', NaN, 'eeNetRefMaster', NaN);
          if ~isempty(s) && ~isempty(q)
            r.eedfGap = max(abs(s.eedf - q.eedf))/max(q.eedf);
            r.meanEnergyGap = abs(meanEnergy(s) - meanEnergy(q))/meanEnergy(q);
            r.eeGainGap = abs(s.power.electronElectronGain - q.power.electronElectronGain)/ ...
              abs(q.power.electronElectronGain);
            r.eeMatrixGap = max(abs(s.M.eeMatrix(:) - q.M.eeMatrix(:)))/max(abs(q.M.eeMatrix(:)));
            r.eeNetRefBranch = abs(s.power.electronElectronNet)/s.power.reference;
            r.eeNetRefMaster = abs(q.power.electronElectronNet)/q.power.reference;
            if ~isempty(prev)
              r.pEedf = log2(prev.eedfGap/r.eedfGap);
              r.pMeanEnergy = log2(prev.meanEnergyGap/r.meanEnergyGap);
              r.pEeGain = log2(prev.eeGainGap/r.eeGainGap);
            end
            prev = r;
          end
          rows{end+1} = r; %#ok<AGROW>
        end
      end
      out = struct2table([rows{:}]);
      disp(out);
      if numel(varargin) > 1
        writeReport(varargin{2}, out, refFolder);
      end
  end

end

function name = caseName(EoN, N)
  name = sprintf('N2_conservative_ee_E%g_N%d', EoN, N);
end

function s = load1(file)
% single-job snapshot of a result file, [] if missing or failed
  s = [];
  if ~isfile(file), return; end
  S = load(file);
  if S.res.ok, s = S.res.jobs{1}; end
end

function m = meanEnergy(s)
  m = sum(s.cell(:).^1.5.*s.eedf(:).*s.du(:));
end

function writeReport(file, T, refFolder)
  fid = fopen(file, 'w');
  cleanup = onCleanup(@() fclose(fid));
  [~, v] = system(sprintf('git -C "%s" rev-parse --short HEAD', refFolder));
  fprintf(fid, '# e-e grid-refinement study (L0)\n\n');
  fprintf(fid, '- Date: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
  fprintf(fid, '- Case: N2_conservative_temporal_ee, uniform grid, no smart grid. Branch 14 (variableGrid: false) vs reference `%s`.\n', strtrim(v));
  fprintf(fid, ['- Gaps between the two codes on the same grid: eedf max|df|/max f; mean energy and e-e gain power, ' ...
    'relative; e-e matrix max|dM|/max|M|. p = log2(gap(N/2)/gap(N)) is the observed order of the gap.\n']);
  fprintf(fid, '- e-e net / ref: |electronElectronNet| / power.reference in each code (energy conservation).\n\n');
  fprintf(fid, ['| E/N (Td) | u_max (eV) | N | du (eV) | eedf gap | p | mean energy gap | p | e-e gain gap | p | ' ...
    'e-e matrix gap | e-e net/ref branch | e-e net/ref master |\n']);
  fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
  for k = 1:height(T)
    fprintf(fid, '| %g | %g | %d | %.3g | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', T.EoN(k), ...
      T.maxEnergy(k), T.N(k), T.du(k), num(T.eedfGap(k)), ord(T.pEedf(k)), num(T.meanEnergyGap(k)), ...
      ord(T.pMeanEnergy(k)), num(T.eeGainGap(k)), ord(T.pEeGain(k)), num(T.eeMatrixGap(k)), ...
      num(T.eeNetRefBranch(k)), num(T.eeNetRefMaster(k)));
  end
end

function s = num(x)
  if isnan(x), s = '-'; else, s = sprintf('%.1e', x); end
end

function s = ord(x)
  if isnan(x), s = '-'; else, s = sprintf('%.2f', x); end
end
