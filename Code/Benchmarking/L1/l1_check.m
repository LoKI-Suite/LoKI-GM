function T = l1_check(refFolder, reportFile)
% L1_CHECK  Compare the L1 runs with the analytic Maxwellians and write a Markdown report.
%
%   T = l1_check(refFolder, reportFile)
%
% Run from the Code folder of the branch-14 tree that holds Output/L1/uniform and Output/L1/veg;
% refFolder is the Code folder of the reference tree (Output/L1/ref). Benchmarking/L1 on the path.
%
% For each case and grid, the analytic f_a(u_n) = C exp(-u_n/T) is sampled at the cell centres and
% normalised like LoKI's EEDF (sum du sqrt(u) f = 1). Error measures:
%   errW      weighted norm  sum du sqrt(u) |f - f_a| / sum du sqrt(u) f_a
%   errMax5   max |f/f_a - 1| over the cells with f_a >= 1e-5 f_a(u_1)   (bulk)
%   errMax15  the same down to 1e-15 (tail)
%   meanE     |<u> - 1.5 T| / 1.5 T
%   p         observed order of errW between N/2 and N
% Acceptance (RESEARCH_PLAN L1): observed order about 2 on smooth grids.

  modes = {'ref', 'uniform', 'veg'};
  folders = {fullfile(refFolder, 'Output', 'L1', 'ref'), fullfile('Output', 'L1', 'uniform'), ...
    fullfile('Output', 'L1', 'veg')};
  cases = l1_cases();
  cellNumbers = [100 200 400 800 1600];

  rows = {};
  for m = 1:numel(modes)
    for k = 1:numel(cases)
      prev = [];
      for N = cellNumbers
        r = struct('mode', modes{m}, 'case', cases(k).name, 'N', N, 'T', NaN, 'errW', NaN, ...
          'errMax5', NaN, 'errMax15', NaN, 'meanE', NaN, 'p', NaN, 'powerBalance', NaN, 'ok', false, ...
          'code', '');
        file = fullfile(folders{m}, sprintf('%s_N%d.mat', cases(k).name, N));
        if isfile(file)
          S = load(file);
          r.code = S.res.codeVersion;
          if S.res.ok
            s = S.res.jobs{end};
            r.ok = true;
            r.T = temperature(cases(k).analytic);
            u = s.cell(:); du = s.du(:); f = s.eedf(:);
            fa = exp(-u/r.T); fa = fa/sum(du.*sqrt(u).*fa);
            r.errW = sum(du.*sqrt(u).*abs(f - fa))/sum(du.*sqrt(u).*fa);
            rel = abs(f./fa - 1);
            r.errMax5 = max(rel(fa >= 1e-5*fa(1)));
            r.errMax15 = max(rel(fa >= 1e-15*fa(1)));
            r.meanE = abs(sum(du.*u.^1.5.*f) - 1.5*r.T)/(1.5*r.T);
            r.powerBalance = abs(s.power.balance)/s.power.reference;
            if ~isempty(prev) && prev.ok
              r.p = log2(prev.errW/r.errW);
            end
          end
        end
        prev = r;
        rows{end+1} = r; %#ok<AGROW>
      end
    end
  end
  T = struct2table([rows{:}]);

  E = l1_ee_operator('study');
  G = table();
  for m = 1:numel(modes)
    file = fullfile(folders{m}, 'maxwell_ee_N400.mat');
    variant = 'branch'; if strcmp(modes{m}, 'ref'), variant = 'master'; end
    g = NaN; if isfile(file), g = l1_ee_operator('guard', file, variant); end
    G = [G; table(modes(m), {variant}, g, 'VariableNames', {'mode', 'copy', 'spread'})]; %#ok<AGROW>
  end

  if nargin > 1
    writeReport(reportFile, T, E, G, refFolder);
  end
end

% ---------------------------------------------------------------------------------------------------
function T = temperature(a)
% analytic temperature (eV); see l1_cases
  kTg = Constant.boltzmannInEV*300;
  switch a.type
    case 'gasTemperature'
      T = kTg;
    case 'constantNu'
      mu = Constant.electronMass/(2*1.008*1.660539040e-27);       % m/M of H2 (Input/Databases/masses.txt)
      N = 133.32/(Constant.boltzmann*300);                           % gas density
      EoN = a.EoNTd*1e-21;
      WoN = 2*pi*a.excFreq/N;
      H = 1 + (WoN/Constant.gamma)^2/a.c^2;
      T = kTg + EoN^2/(6*mu*a.c^2*H);
  end
end

function writeReport(file, T, E, G, refFolder)
  fid = fopen(file, 'w');
  cleanup = onCleanup(@() fclose(fid));
  [~, v] = system(sprintf('git -C "%s" rev-parse --short HEAD', refFolder));
  codes = unique(T.code(~strcmp(T.code, '')));
  fprintf(fid, '# L1 report: analytic Maxwellians\n\n');
  fprintf(fid, '- Date: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
  fprintf(fid, '- Code versions recorded by the runs: %s (reference: `%s`)\n', strjoin(codes, ', '), strtrim(v));
  fprintf(fid, ['- Modes: *ref* reference code, uniform grid; *uniform* branch 14, uniform grid; *veg* branch 14, ' ...
    'geometric grid with first step maxEnergy/(10N).\n']);
  fprintf(fid, ['- errW = sum du sqrt(u)|f-f_a| / sum du sqrt(u) f_a; errMax5/errMax15 = max |f/f_a-1| down to ' ...
    'f_a = 1e-5 / 1e-15 of f_a(u_1); meanE = |<u>-1.5T|/1.5T; p = observed order of errW; power bal. = ' ...
    '|balance|/reference.\n\n']);

  fprintf(fid, '## Solved EEDF vs analytic Maxwellian\n\n');
  fprintf(fid, '| case | mode | N | T (eV) | errW | p | errMax5 | errMax15 | meanE | power bal. |\n');
  fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|\n');
  for k = 1:height(T)
    if ~T.ok(k)
      fprintf(fid, '| %s | %s | %d | run failed or missing | | | | | | |\n', T.case{k}, T.mode{k}, T.N(k));
      continue
    end
    fprintf(fid, '| %s | %s | %d | %.6g | %s | %s | %s | %s | %s | %s |\n', T.case{k}, T.mode{k}, T.N(k), ...
      T.T(k), num(T.errW(k)), ord(T.p(k)), num(T.errMax5(k)), num(T.errMax15(k)), num(T.meanE(k)), ...
      num(T.powerBalance(k)));
  end

  fprintf(fid, '\n## e-e operator on a discrete Maxwellian\n\n');
  fprintf(fid, ['Copies of the assembly of each code (see l1_ee_operator). Grid u_max = 40 T. detailed balance = ' ...
    'max|J|/max(Ahat f) for f = discrete Maxwellian at T; *with Te = T*: diagnostic, Te forced to T in branch 14''s weights; ' ...
    'Te/T = Te computed by the code from the discrete Maxwellian at T; particles/energy for f = exp(-(u/T)^2).\n\n']);
  fprintf(fid, '| code | grid | T (eV) | N | detailed balance | p | detailed balance with Te = T | Te/T - 1 | particles | energy |\n');
  fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|\n');
  for k = 1:height(E)
    fprintf(fid, '| %s | %s | %.4g | %d | %s | %s | %s | %s | %s | %s |\n', E.variant{k}, E.grid{k}, E.T(k), ...
      E.N(k), num(E.detailedBalance(k)), ord(E.pDetailedBalance(k)), num(E.detailedBalanceTe(k)), ...
      num(E.TeOverT(k) - 1), num(E.particles(k)), num(E.energy(k)));
  end

  fprintf(fid, '\n## Guard: copied e-e assembly vs the code (maxwell_ee, N = 400)\n\n');
  fprintf(fid, ['max relative spread of (flux vector stored by LoKI)/(copy applied to the converged EEDF), ' ...
    'after removing alpha. Small values (<= ~1e-8, the EEDF change of the last iteration) mean the copy ' ...
    'is the code.\n\n']);
  fprintf(fid, '| mode | copy | spread |\n|---|---|---|\n');
  for k = 1:height(G)
    fprintf(fid, '| %s | %s | %s |\n', G.mode{k}, G.copy{k}, num(G.spread(k)));
  end
end

function s = num(x)
  if isnan(x), s = '-'; else, s = sprintf('%.1e', x); end
end

function s = ord(x)
  if isnan(x), s = '-'; else, s = sprintf('%.2f', x); end
end
