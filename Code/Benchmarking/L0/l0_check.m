function T = l0_check(refFolder, reportFile)
% L0_CHECK  L0 checks on the results written by l0_run.
%
%   T = l0_check(refFolder)              run from branch-14 Code folder
%   T = l0_check(refFolder, reportFile)  also writes a Markdown report
%
% refFolder: Code folder of the reference tree (e.g. a worktree of master) in which l0_run('ref')
% was run. Results of the current tree are read from Output/L0/uniform and Output/L0/veg.
%
% Checks (per case, worst over all reduced-field jobs):
%   particles  max over operators and source cells m of |sum_n du_n M(n,m)| / sum_n du_n |M(n,m)|,
%              for the particle-conserving operators (field, elastic, CAR, e-e, inelastic+superelastic,
%              conservative ionization and attachment). Matrix rows are per unit energy of cell n.
%   power      relative power balance reported by LoKI-B
%   e-e energy (e-e cases) energy conservation of the e-e operator at the final eedf f:
%              operator  |sum_n du_n u_n (Mee f)_n| / sum_n du_n u_n |(Mee f)_n|
%              power     |electronElectronNet| / power.reference
%              (Mee is the operator of the last e-e iteration, built from the previous iterate, so
%              both are limited by LoKI's e-e convergence criterion |Pee/Pref| < maxPowerBalanceRelError)
%              (power.reference is LoKI's total power gain, Boltzmann.evaluatePower)
%   regression (uniform mode only) branch-14 uniform vs reference code:
%              matrices  max |M - Mref| / max |Mref|, per operator
%              eedf      max |f - fref| / max fref
%              outputs   max relative difference over all numeric fields of power, swarm
%                        parameters and rate coefficients (ignoring values below 1e-30 of the
%                        field's scale), except the explained ones listed in 'known' (reported in
%                        the notes column). Power terms that vanish by construction (electronElectronNet)
%                        or are numerically zero in both codes (below tol.power*power.reference) are
%                        compared as |difference| / power.reference.

  tol.particles = 1e-12;
  tol.power = 1e-9;
  tol.matrix = 5e-13;
  tol.eedf = 1e-9;
  tol.outputs = 1e-8;

  continuousOps = {'fieldMatrix', 'elasticMatrix', 'CARMatrix', 'eeMatrix'};
  discreteOps = {'inelasticMatrix', 'ionizationConservativeMatrix', 'attachmentConservativeMatrix'};
  % residuals checked as identities, not compared between codes
  skipOutputs = {'.power.balance', '.power.relativeBalance'};
  % power terms that vanish by construction (compared against power.reference)
  zeroPower = {'electronElectronNet'};

  % Explained differences (labelled in the notes column; tolerances are not relaxed).
  % Update this list when a difference is resolved or a new one is understood.
  % (redTownsendCoeff, redAttCoeff and driftVelocity depend on D in the spatial growth model)
  known.outputs = {'.swarm.redDiffCoeff', '.swarm.redDiffCoeffEnergy', '.swarm.characEnergy', ...
    '.swarm.redTownsendCoeff', '.swarm.redAttCoeff', '.swarm.driftVelocity'};
  known.outputsNote = ['D, D_energy use node quadrature on branch 14 (fef9b1e), also changing u_k ' ...
    'and, with spatial growth, alpha/N, eta/N, v_d - T1.5/DP-5'];
  known.eeNote = 'e-e operator reformulated on branch 14 (aux. array sign, Maxwellian weights) - T1.4';
  known.eeOps = {'eeMatrix'};   % e-e matrix differences, same reason (own column)
  known.vegDiscNote = 'discrete operators still use the uniform index shift - P3.5';

  cases = l0_cases();
  rows = {};
  for k = 1:numel(cases)
    name = cases(k).name;
    for mode = {'uniform', 'veg'}
      file = fullfile('Output', 'L0', mode{1}, [name '.mat']);
      r = struct('case', name, 'mode', mode{1}, 'code', '', 'status', 'missing', ...
        'partCont', NaN, 'worstCont', '', 'partDisc', NaN, 'worstDisc', '', 'power', NaN, ...
        'eeEnergyOp', NaN, 'eeEnergyPower', NaN, 'matrix', NaN, 'worstMatrix', '', 'eeMatrix', NaN, ...
        'eedf', NaN, 'outputs', NaN, 'worstOutput', '', 'knownOutputs', NaN, 'time', NaN, 'notes', '');
      if ~isfile(file)
        rows{end+1} = r; %#ok<AGROW>
        continue;
      end
      S = load(file); res = S.res;
      r.code = res.codeVersion;
      r.time = res.time;
      if ~res.ok
        r.status = 'run failed';
        r.notes = errorWhere(res.error);
        refFile = fullfile(refFolder, 'Output', 'L0', 'ref', [name '.mat']);
        if strcmp(mode{1}, 'uniform') && isfile(refFile)
          R = load(refFile);
          if ~R.res.ok, r.notes = [r.notes '; reference code fails too']; end
        end
        rows{end+1} = r; %#ok<AGROW>
        continue;
      end

      % --- identities
      [r.partCont, r.partDisc, r.power] = deal(0);
      for j = 1:numel(res.jobs)
        s = res.jobs{j};
        [e, op] = worstDefect(s, continuousOps);
        if e > r.partCont, r.partCont = e; r.worstCont = op; end
        [e, op] = worstDefect(s, discreteOps);
        if e > r.partDisc, r.partDisc = e; r.worstDisc = op; end
        r.power = max(r.power, s.power.relativeBalance);
        if isfield(s.M, 'eeMatrix') && ~isempty(s.M.eeMatrix)
          r.eeEnergyOp = max([r.eeEnergyOp, eeEnergyDefect(s)], [], 'omitnan');
          r.eeEnergyPower = max([r.eeEnergyPower, ...
            abs(s.power.electronElectronNet)/s.power.reference], [], 'omitnan');
        end
      end

      % --- regression against the reference code (uniform grid only)
      if strcmp(mode{1}, 'uniform')
        refFile = fullfile(refFolder, 'Output', 'L0', 'ref', [name '.mat']);
        if isfile(refFile)
          R = load(refFile); ref = R.res;
          if ~ref.ok
            r.status = 'ref failed';
          elseif numel(ref.jobs) ~= numel(res.jobs)
            r.status = 'job count differs';
          else
            [r.matrix, r.eedf, r.outputs, r.knownOutputs] = deal(0);
            if cases(k).ee, r.eeMatrix = 0; end
            for j = 1:numel(res.jobs)
              s = res.jobs{j}; q = ref.jobs{j};
              if numel(s.node) ~= numel(q.node) || max(abs(s.node - q.node)) > 1e-12*max(q.node)
                r.status = 'grids differ'; break;
              end
              for op = fieldnames(q.M)'
                A = q.M.(op{1});
                if isempty(A) || ~isfield(s.M, op{1}), continue; end
                B = s.M.(op{1});
                if isempty(B), B = zeros(size(A)); end
                d = max(abs(B(:) - A(:)))/max(max(abs(A(:))), realmin);
                if ismember(op{1}, known.eeOps)
                  r.eeMatrix = max(r.eeMatrix, d);
                elseif d > r.matrix
                  r.matrix = d; r.worstMatrix = op{1};
                end
              end
              r.eedf = max(r.eedf, max(abs(s.eedf - q.eedf))/max(q.eedf));
              [d, where] = powerDiff(s.power, q.power, skipOutputs, zeroPower, tol.power);
              if d > r.outputs, r.outputs = d; r.worstOutput = where; end
              a = struct('swarm', s.swarm, 'rates', s.rates);
              b = struct('swarm', q.swarm, 'rates', q.rates);
              [d, where] = structDiff(a, b, [skipOutputs, known.outputs]);
              if d > r.outputs, r.outputs = d; r.worstOutput = where; end
              r.knownOutputs = max(r.knownOutputs, structDiff(a, b, {}, known.outputs));
            end
          end
        else
          r.status = 'no ref';
        end
      end

      if strcmp(r.status, 'missing')
        pass = r.partCont <= tol.particles && r.partDisc <= tol.particles && r.power <= tol.power;
        if cases(k).ee
          % energy conservation of e-e judged by |Pee|/Pref (LoKI's own criterion); the operator-level
          % defect is information only (Mee is built from the previous iterate)
          pass = pass && r.eeEnergyPower <= tol.power;
        end
        if strcmp(mode{1}, 'uniform')
          pass = pass && r.matrix <= tol.matrix && r.eedf <= tol.eedf && r.outputs <= tol.outputs;
          if cases(k).ee, pass = pass && r.eeMatrix <= tol.matrix; end
        end
        r.status = ternary(pass, 'PASS', 'FAIL');
      end

      % notes for explained differences
      notes = {};
      if r.partDisc > tol.particles && strcmp(mode{1}, 'veg')
        notes{end+1} = known.vegDiscNote; %#ok<AGROW>
      end
      if r.knownOutputs > tol.outputs
        notes{end+1} = sprintf('%s (max %.1e)', known.outputsNote, r.knownOutputs); %#ok<AGROW>
      end
      if cases(k).ee && strcmp(mode{1}, 'uniform') && ...
          (r.eeMatrix > tol.matrix || r.eedf > tol.eedf || r.outputs > tol.outputs)
        notes{end+1} = known.eeNote; %#ok<AGROW>
      end
      r.notes = strjoin(notes, '; ');
      rows{end+1} = r; %#ok<AGROW>
    end
  end

  T = struct2table([rows{:}]);
  disp(T(:, {'case', 'mode', 'status', 'partCont', 'partDisc', 'power', 'eeEnergyOp', ...
    'eeEnergyPower', 'matrix', 'eeMatrix', 'eedf', 'outputs', 'worstOutput', 'time'}));

  if nargin > 1
    writeReport(reportFile, T, tol, refFolder);
  end

end

function [eMax, worst] = worstDefect(s, ops)
  eMax = 0; worst = '';
  for op = ops
    if isfield(s.M, op{1}) && ~isempty(s.M.(op{1})) && any(s.M.(op{1})(:))
      e = particleDefect(s.M.(op{1}), s.du);
      if e > eMax, eMax = e; worst = op{1}; end
    end
  end
end

function e = particleDefect(M, du)
  du = du(:);
  num = abs(du.'*M);
  den = du.'*abs(M);
  active = den > 0;
  e = max([0, num(active)./den(active)]);
end

function t = errorWhere(report)
% first line of an error report and the first LoKI location ('Error in ...') below it
  lines = strsplit(report, newline);
  t = strtrim(lines{1});
  k = find(startsWith(strtrim(lines), 'Error in'), 1);
  if ~isempty(k), t = sprintf('%s (%s)', t, strtrim(lines{k})); end
  t = strrep(t, '|', '/');
end

function e = eeEnergyDefect(s)
  du = s.du(:); u = s.cell(:); f = s.eedf(:);
  g = s.M.eeMatrix*f;
  e = abs(sum(du.*u.*g))/max(sum(du.*u.*abs(g)), realmin);
end

function [d, where] = powerDiff(p, q, skip, zeroTerms, tolPower)
% Relative differences of the power terms. Terms that vanish by construction (zeroTerms) or are
% numerically zero in both codes (|value| <= tolPower*reference) are compared against the reference
% power of the reference code (LoKI's total power gain), not against themselves.
  d = 0; where = '';
  ref = q.reference;
  lp = leaves(p, '.power'); lq = leaves(q, '.power');   % includes the per-gas terms (.power.gases)
  for k = 1:size(lq, 1)
    name = lq{k, 1};
    if any(endsWith(name, skip)), continue; end
    i = find(strcmp(lp(:, 1), name), 1);
    if isempty(i), continue; end
    x = lp{i, 2}; y = lq{k, 2};
    if ~isequal(size(x), size(y)) || isempty(y), continue; end
    if any(endsWith(name, strcat('.', zeroTerms))) || max(abs([x(:); y(:)])) <= tolPower*ref
      e = max(abs(x(:) - y(:)))/ref;
    else
      e = max(abs(x(:) - y(:))./max(abs(x(:)), abs(y(:))));
    end
    if e > d, d = e; where = name; end
  end
end

function [d, where] = structDiff(a, b, skip, only)
% Largest relative difference over the numeric leaves shared by a and b. Leaves ending in one of
% skip are ignored; if only is given, just the leaves ending in one of only are compared.
  d = 0; where = '';
  la = leaves(a, ''); lb = leaves(b, '');
  for k = 1:size(la, 1)
    if any(endsWith(la{k, 1}, skip)), continue; end
    if nargin > 3 && ~any(endsWith(la{k, 1}, only)), continue; end
    i = find(strcmp(lb(:, 1), la{k, 1}), 1);
    if isempty(i), continue; end
    x = la{k, 2}; y = lb{i, 2};
    if ~isequal(size(x), size(y))
      d = Inf; where = [la{k, 1} ' (size)']; return;
    end
    scale = max(abs(y(:)));
    if scale == 0, scale = max(abs(x(:))); end
    if scale == 0, continue; end
    mask = abs(y(:)) > 1e-30*scale | abs(x(:)) > 1e-30*scale;
    rel = abs(x(mask) - y(mask))./max(abs(y(mask)), abs(x(mask)));
    if ~isempty(rel) && max(rel) > d
      d = max(rel); where = la{k, 1};
    end
  end
end

function L = leaves(v, prefix)
  L = cell(0, 2);
  if isstruct(v)
    for i = 1:numel(v)
      for f = fieldnames(v)'
        p = prefix;
        if numel(v) > 1, p = sprintf('%s(%d)', p, i); end
        L = [L; leaves(v(i).(f{1}), [p '.' f{1}])]; %#ok<AGROW>
      end
    end
  elseif isnumeric(v) || islogical(v)
    L = {prefix, double(v)};
  end
end

function out = ternary(c, a, b)
  if c, out = a; else, out = b; end
end

function writeReport(file, T, tol, refFolder)
  fid = fopen(file, 'w');
  cleanup = onCleanup(@() fclose(fid));
  fprintf(fid, '# L0 report\n\n');
  fprintf(fid, '- Date: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
  fprintf(fid, '- Code under test (as recorded by l0_run): `%s`\n', strjoin(unique(T.code(~cellfun(@isempty, T.code))), '`, `'));
  [~, v] = system(sprintf('git -C "%s" rev-parse --short HEAD', refFolder));
  fprintf(fid, '- Reference code: `%s`\n', strtrim(v));
  fprintf(fid, ['- Tolerances: particles %.0e, power balance %.0e, matrices %.0e, eedf %.0e, ' ...
    'outputs %.0e\n'], tol.particles, tol.power, tol.matrix, tol.eedf, tol.outputs);
  fprintf(fid, '- Columns (worst over the reduced-field jobs of each case):\n');
  fprintf(fid, ['  - *particles cont./disc.*: particle-conservation defect of the continuous (field, elastic, ' ...
    'CAR, e-e) and discrete (inelastic+superelastic, conservative ionization/attachment) operators.\n']);
  fprintf(fid, '  - *power bal.*: relative power balance reported by LoKI.\n');
  fprintf(fid, ['  - *e-e net/ref*: |electronElectronNet| / power.reference (LoKI''s total power gain), the ' ...
    'e-e energy-conservation check; *e-e energy op.*: energy defect of the e-e operator applied to the ' ...
    'final eedf, for information only (the operator is built from the previous iterate).\n']);
  fprintf(fid, '  - *run failed*: the notes give the error and where it occurred.\n');
  fprintf(fid, ['  - *matrices vs ref*: max|M-Mref|/max|Mref| over all operators except e-e, which is in ' ...
    '*e-e matrix vs ref*.\n']);
  fprintf(fid, ['  - *outputs vs ref*: max relative difference over power terms, swarm parameters and rate ' ...
    'coefficients, **excluding** the D-related outputs (D, D_energy, u_k, alpha/N, eta/N, v_d), whose ' ...
    'deviation is given in the notes. Power terms that vanish by construction or are numerically zero ' ...
    'are compared against power.reference.\n\n']);
  fprintf(fid, ['| case | grid | status | particles cont. | particles disc. | power bal. | e-e energy op. | ' ...
    'e-e net/ref | matrices vs ref | e-e matrix vs ref | eedf vs ref | outputs vs ref | worst output | ' ...
    'time (s) | notes |\n']);
  fprintf(fid, '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n');
  for k = 1:height(T)
    fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %.0f | %s |\n', ...
      T.case{k}, T.mode{k}, T.status{k}, num(T.partCont(k)), num(T.partDisc(k)), num(T.power(k)), ...
      num(T.eeEnergyOp(k)), num(T.eeEnergyPower(k)), num(T.matrix(k)), num(T.eeMatrix(k)), ...
      num(T.eedf(k)), num(T.outputs(k)), strrep(T.worstOutput{k}, '|', '/'), T.time(k), T.notes{k});
  end
end

function s = num(x)
  if isnan(x), s = '-'; else, s = sprintf('%.1e', x); end
end
