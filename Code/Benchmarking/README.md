# Benchmarking (variable energy grid)

Regression and verification harness for the variable energy grid work (branch `14-variable_energy_grid`).
The test levels L0–L4 are defined in the project's research plan; L0 and L1 are implemented here.

| Level | Folder | Purpose |
|---|---|---|
| L0 | `L0/` | Regression: branch 14 against master on 26 cases (gases, operators, e–e, HF, smart grid, temporal integration) at 0.01–1000 Td, on uniform and variable grids. Compares grids, matrices, EEDF, power balance, swarm parameters and rate coefficients. |
| L1 | `L1/` | Verification against analytic Maxwellians: constant collision frequency (DC and HF), E/N = 0 with real H2 elastic data, + CAR, + e–e; grid refinement N = 100…1600 and observed order. Also applies the e–e operator to a discrete Maxwellian. |

## Rules

- Run the tests in a separate git worktree checked out at the commit under test, never in a working clone, and record the commit in every result.
- Run from the `Code` folder of the tree under test, with the harness folder of branch 14 on the path (the reference tree, e.g. master, does not contain it).
- Results go to `Output/L0/...` and `Output/L1/...` (git-ignored).

## Usage

```matlab
% L0: in each tree (mode 'ref' on master, 'uniform' and 'veg' on branch 14)
cd <tree>/Code; addpath('<branch14>/Code/Benchmarking/L0'); l0_run('ref')
% then, in the branch-14 tree
T = l0_check('<master>/Code', 'L0_report.md');

% L1: same pattern
cd <tree>/Code; addpath('<branch14>/Code/Benchmarking/L0', '<branch14>/Code/Benchmarking/L1');
l1_run('ref')                              % 'uniform' and 'veg' on branch 14 (writes its LXCat files itself)
T = l1_check('<master>/Code', 'L1_report.md');
```

See the help text of each function for details (`help l0_run`, `help l1_check`, ...).
