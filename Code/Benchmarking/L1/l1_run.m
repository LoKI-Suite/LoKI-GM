function l1_run(mode, caseNames)
% L1_RUN  Run the L1 analytic cases on a sequence of grids in the code tree of the current folder.
%
%   l1_run(mode)             all cases of l1_cases()
%   l1_run(mode, caseNames)  only the named cases
%
% mode is 'ref' (reference code, uniform grid), 'uniform' or 'veg' (branch 14), as in L0. Call with
% the current folder set to the Code folder of the tree under test, with Benchmarking/L0 and
% Benchmarking/L1 of branch 14 on the path, in a worktree at the commit under test (D-007):
%
%   cd <tree>/Code; addpath('<branch14>/Code/Benchmarking/L0', '<branch14>/Code/Benchmarking/L1'); l1_run('ref')
%
% Each case runs on N = 100, 200, ..., 1600 cells at fixed maxEnergy, without smart grid. On the variable
% grid ('veg') the geometric grid of branch 14 is used with firstEnergyStep = maxEnergy/(10 N), so the
% ratio of the last to the first step is about the same for every N. Results: Output/L1/<mode>/<case>_N<N>.mat
% (struct 'res' from l0_run_case).

  cellNumbers = [100 200 400 800 1600];

  l1_data('write');
  cases = l1_cases();
  if nargin > 1
    cases = cases(ismember({cases.name}, caseNames));
  end
  outDir = fullfile('Output', 'L1', mode);
  for k = 1:numel(cases)
    for N = cellNumbers
      c = cases(k);
      c.grid = struct('maxEnergy', c.maxEnergy, 'cellNumber', N);
      if strcmp(mode, 'veg')
        c.grid.firstEnergyStep = c.maxEnergy/(10*N);
      end
      c.name = sprintf('%s_N%d', cases(k).name, N);
      l0_run_case(c, mode, outDir);
    end
  end

end
