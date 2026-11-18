function res = l0_run_case(c, mode, outDir)
% L0_RUN_CASE  Run one case (struct from l0_cases) in the code tree of the current folder.
%
%   res = l0_run_case(c, mode, outDir)
%
% Writes the setup file to Output/L0/setups and the result to <outDir>/<c.name>.mat (struct 'res').
% The current folder must be the Code folder of the tree under test (see l0_run).

  if ~isfolder('Input') || ~isfile('loki.m')
    error('l0_run:folder', 'The current folder must be the Code folder of a LoKI-GM tree.');
  end
  setupDir = fullfile('Output', 'L0', 'setups');
  if ~isfolder(setupDir), mkdir(setupDir); end
  if ~isfolder(outDir), mkdir(outDir); end

  setupName = sprintf('%s_%s.in', c.name, mode);
  writelines(l0_setup_text(c, mode), fullfile(setupDir, setupName));
  % Parse reads setup files relative to Code/Input
  setupFile = ['..' filesep setupDir filesep setupName];

  % LoKI-GM keeps object IDs in persistent variables; loki.m clears them before each run
  clear functions %#ok<CLFUNC>

  res = struct('case', c, 'mode', mode, 'ok', false, 'error', '', 'log', '', 'time', NaN, ...
    'jobs', {{}}, 'codeVersion', codeVersion());
  t = tic;
  try
    [res.log, res.jobs] = evalc('solveJobs(setupFile)');
    res.ok = true;
  catch ME
    res.error = getReport(ME, 'extended', 'hyperlinks', 'off');
  end
  res.time = toc(t);
  save(fullfile(outDir, [c.name '.mat']), 'res');
  if res.ok
    fprintf('%-40s %-8s ok     %3d jobs %7.1f s\n', c.name, mode, numel(res.jobs), res.time);
  else
    fprintf('%-40s %-8s FAILED %s\n', c.name, mode, strtok(res.error, newline));
  end

end

function jobs = solveJobs(setupFile)
% Same sequence as loki.m, without output, recording a snapshot after each job.
  setup = Setup(setupFile);
  [electronKinetics, ~] = setup.initializeSimulation();
  jobs = {};
  while setup.currentJobID <= setup.numberOfJobs
    electronKinetics.solve();
    jobs{end+1} = snapshot(electronKinetics); %#ok<AGROW>
    setup.nextJob();
  end
end

function s = snapshot(ek)
  g = ek.energyGrid;
  s.reducedField = ek.workCond.reducedField;
  s.node = g.node;
  s.cell = g.cell;
  variable = isprop(g, 'energyStep') && ~isempty(g.energyStep);
  if variable
    s.du = g.energyStep;
  else
    s.du = g.step*ones(size(g.cell));
  end
  s.eedf = ek.eedf;
  s.power = ek.power;
  s.swarm = ek.swarmParam;
  s.rates = ek.rateCoeffAll;
  s.M = struct();
  for name = {'fieldMatrix', 'elasticMatrix', 'CARMatrix', 'inelasticMatrix', ...
      'ionizationConservativeMatrix', 'ionizationMatrix', 'attachmentConservativeMatrix', 'attachmentMatrix'}
    if isprop(ek, name{1})
      s.M.(name{1}) = ek.(name{1});
    end
  end
  s.M.eeMatrix = eeMatrix(ek, s.du, variable);
end

function Mee = eeMatrix(ek, du, variable)
% e-e operator used in the last iteration of the stationary e-e solver (Boltzmann.solveEEColl),
% rebuilt from the up/downflux vectors stored by the solver, Afinal and Bfinal (Aee and Bee are kept
% only with non-conservative ionization/attachment):
%   branch 14 (energyStep):  Afinal = alpha*Aee*f          Mee(n,:) divided by du_n
%   master (uniform step):   Afinal = alpha/du*Aee*f       Mee(n,:) not divided
% so with A = Afinal (branch) or du*Afinal (master) both are
%   Mee(n,n) = -(A_n+B_n)/du_n,  Mee(n,n-1) = A_(n-1)/du_n,  Mee(n,n+1) = B_(n+1)/du_n.
% Valid for the stationary algorithms only (branch 14's temporalIntegration stores Afinal/du).
% Empty when e-e collisions are off.
  Mee = [];
  if ~ek.includeEECollisions || isempty(ek.Afinal), return; end
  du = du(:);
  A = ek.Afinal(:);
  B = ek.Bfinal(:);
  if ~variable
    A = A*du(1);
    B = B*du(1);
  end
  n = numel(du);
  Mee = diag(-(A + B)./du) + diag(A(1:n-1)./du(2:n), -1) + diag(B(2:n)./du(1:n-1), 1);
end

function v = codeVersion()
  [status, out] = system('git rev-parse --short HEAD');
  if status == 0
    v = strtrim(out);
    [~, dirty] = system('git status --porcelain --untracked-files=no');
    if ~isempty(strtrim(dirty)), v = [v '+modified']; end
  else
    v = 'unknown';
  end
end
