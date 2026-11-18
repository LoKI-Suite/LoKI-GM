function l0_run(mode, caseNames)
% L0_RUN  Run the L0 cases in the code tree of the current folder and store the results.
%
%   l0_run(mode)             runs all cases of l0_cases()
%   l0_run(mode, caseNames)  runs only the named cases (cell array of names)
%
% Must be called with the current folder set to the 'Code' folder of the tree under test, with
% Benchmarking/L0 of branch 14 on the path (the reference tree, e.g. master, does not contain it):
%
%   cd <tree>/Code; addpath('<branch14>/Code/Benchmarking/L0'); l0_run('ref')
%
% Run the tests in a worktree checked out at the commit under test, not in the author's clone (D-007).
%
% mode is 'ref', 'uniform' or 'veg' (see l0_setup_text). For each case, one file
% Output/L0/<mode>/<case>.mat is written (Output/ is git-ignored). It holds a struct 'res' with
% one snapshot per reduced-field job (grid, eedf, power, swarm, rate coefficients, matrices).

  cases = l0_cases();
  if nargin > 1
    cases = cases(ismember({cases.name}, caseNames));
  end
  outDir = fullfile('Output', 'L0', mode);
  for k = 1:numel(cases)
    l0_run_case(cases(k), mode, outDir);
  end

end
