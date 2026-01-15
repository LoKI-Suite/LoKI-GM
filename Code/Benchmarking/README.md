# 📊 Benchmarking Variable Energy Grid

## 🎯 Purpose

**This benchmarking suite is designed to validate the variable energy grid implementation** in LoKI-B. The benchmarks compare results from:
- **Uniform energy grids** (fixed step size)
- **Variable energy grids** (adaptive step size)

The goal is to ensure that the variable energy grid implementation produces accurate results while potentially using fewer grid points, improving computational efficiency.

---

## 📦 Installation

**⚠️ IMPORTANT:** Before running the benchmarks, you need to copy files to the correct folders in the LoKI-B code:

### 1. Copy Files from `Code/` Folder

Copy all files from the `Code/` folder (in this benchmarking folder) to the main LoKI-B `Code/` folder:

```
LoKI-B/Code/  ← Copy files here
```

**Files to copy:**
- All `.m` files necessary for the benchmarks
- Auxiliary functions and generation scripts

### 2. Copy Files from `Input/` Folder

Copy all files from the `Input/` folder (in this benchmarking folder) to the main LoKI-B `Input/` folder:

```
LoKI-B/Code/Input/  ← Copy files here
```

**Files to copy:**
- Input `.in` files for the benchmarks
- Generated cross sections
- Other necessary configuration files

---

## 🚀 How to Use

**⚠️ IMPORTANT:** All benchmarking scripts are in the `Benchmarking/` folder. Execute them from the project root (the `Code/` folder) so that relative paths work correctly.

```matlab
% Make sure you are in the project root
cd('C:/path/to/LoKI-B/Code')  % adjust the path as needed

% Scripts can be called directly:
test_benchmark_installation.m
```


## 📋 Execution Order

### 1️⃣ **Installation and Verification** (run once)

```matlab
matlab -nodesktop -nosplash -r "test_benchmark_installation; exit"
```

Verifies that all necessary files exist.

---

### 2️⃣ **Generate Cross Sections** (optional, run once or when changing parameters)

If you successfully copied `/Input` folder, this step is not needed, since the files to generate are already there.

(in MATLAB already, run:)

```matlab
% Cross section Maxwellian_const_v (nu = const)
generate_maxwellian_const_v_cross_section

% Cross section Dummy (sigma ≈ 0, for e-e test)
generate_dummy_elastic
```

**Outputs:**
- `Input/Maxwellian_const_v/constant_nu_elastic.txt`
- `Input/Dummy/H2_dummy_elastic.txt`

---

### 3️⃣ **Maxwellian_const_v Diagnosis** (optional, but recommended)

```matlab
diagnose_maxwellian_const_v
```

**Verifies:**
- ✓ Cross section has correct order of magnitude
- ✓ Collision frequency ν is constant
- ✓ Generates verification plots

**Output:** `Input/Maxwellian_const_v/cross_section_diagnostic.png`

---

### 4️⃣ **Run Complete Benchmarks**

```matlab
run_all_benchmarks
```

**Executes 4 types of tests:**

1. **Grid Comparison** (14 tests)
   - Fixed delta u, variable N: 50, 100, 200, 400
   - Fixed N, variable delta u: 5e-4, 1e-3, 5e-3

2. **Maxwellian Elastic** (2 tests)
   - E/N = 0, elastic collisions only
   - Solution: Maxwellian at Tg = 300 K ≈ 0.0259 eV

3. **Maxwellian e-e** (2 tests)
   - E/N = 10 Td, e-e collisions only (gas ≈ 0)
   - Solution: Maxwellian at Te ≈ 1-2 eV

4. **Maxwellian_const_v** (2 tests)
   - E/N = 100 Td, elastic collisions with nu=const
   - Solution: Maxwellian at $T_{eff}$

**Output:** `Output/Output/comprehensive_benchmark/`

---

### 5️⃣ **Analyze Results**

```matlab
Benchmarking/analyze_all_benchmarks
```

**Generates:**
- Convergence plots
- Comparisons with analytical solutions
- Relative errors
- Text report

**Outputs:**
- `Output/comprehensive_benchmark/figures/*.png`
- `Output/comprehensive_benchmark/benchmark_report.txt`

---

## 📊 Test Structure

### Differences between Maxwellian Tests

| Test | E/N | Gas Collisions | e-e Collisions | Solution |
|------|-----|----------------|----------------|----------|
| **Elastic** | 0 Td | **YES** (H2_elastic) | NO | Maxwellian at **Tg ≈ 0.026 eV** |
| **e-e** | 10 Td | **NO** (dummy ≈ 0) | **YES** | Maxwellian at **Te ≈ 1-2 eV** |


---

## 🎯 Expected Results

### Maxwellian Elastic (E/N=0)
- EEDF decays **exponentially** from low energy
- $T_{e}$ ≈ $T_{g}$ ≈ 0.026 eV (300 K)

### Maxwellian e-e (E/N=10 Td)
- EEDF decays **exponentially** but with HIGHER energy
- $T_{e}$ ≈ 1-2 eV (heating by field + e-e collisions)
- **Different** from the elastic test!

### Maxwellian_const_v (E/N=100 Td, nu=const)
- EEDF decays **exponentially** (Maxwellian)
- $T_{eff}$ ≈ 2-3 eV


---

## ⚠️ Troubleshooting

### EEDFs appear as horizontal lines
→ **Problem:** Cross section too large or wrong configuration  
→ **Solution:** Regenerate cross sections with correct parameters

### Error: "Vibrational distribution not normalized"
→ **Problem:** Using `H2_LXCat.txt` (has vibrational) without defining populations  
→ **Solution:** Use `H2_elastic_LXCat.txt` (elastic only)

### Maxwellian tests give identical results
→ **Problem:** Both configured with E/N=0 and elastic collisions  
→ **Solution:** e-e test must have E/N≠0 and dummy cross section

### "Zero found in interval [1, 3]"
→ **This is not an error!** It's the success of `fzero` finding the ratio 'a'

---

## 📁 File Structure

### Main Scripts
- `run_all_benchmarks.m` - Runs all benchmarks
- `analyze_all_benchmarks.m` - Analyzes results

### Test Scripts
- `test_benchmark_installation.m` - Verifies installation
- `test_single_simulation.m` - Single simulation test


### Generation Scripts
- `generate_maxwellian_const_v_cross_section.m` - Cross section nu=const
- `generate_dummy_elastic.m` - Dummy cross section (≈0)

### Analysis and Comparison
- `analyze_all_benchmarks.m` - Results analysis

### Diagnostics
- `test_benchmark_installation.m` - Initial verification
- `diagnose_maxwellian_const_v.m` - Cross section diagnosis

### Analytical Functions
- `analytical_maxwellian.m` - Maxwellian EEDF (in root or Benchmarking/)
- `analytical_maxwellian_const_v.m` - Maxwellian EEDF (nu=const)



---

## 🔬 Critical Parameters

### Grid.m Constraints
```
firstEnergyStep × cellNumber < maxEnergy
```

Valid example:
```yaml
firstEnergyStep: 1e-3
cellNumber: 1000
maxEnergy: 5
# 1e-3 × 1000 = 1 < 5 ✓
```

### Cross Section Values
- **H2 real**: ~1×10⁻¹⁹ m²
- **Maxwellian_const_v**: ~1×10⁻²¹ - 1×10⁻²⁰ m² (smaller, OK)
- **Dummy**: ~1×10⁻³⁰ m² (negligible)

---

**Last updated:** January 15, 2026
