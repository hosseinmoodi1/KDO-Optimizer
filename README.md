# Karma-Dharma Optimizer (KDO) – Complete Implementation

This repository contains the full MATLAB implementation of the **Karma-Dharma Optimizer (KDO)** , a memory-adaptive metaheuristic for high-dimensional continuous optimisation and robust VLSI design. The code includes all benchmark scripts for CEC2017, CEC2022, and four hard electronics engineering problems, as well as ablation study variants and full statistical output generation (Friedman ranks, paired Wilcoxon tests with Holm–Bonferroni correction, effect sizes, convergence curves, boxplots, performance profiles, and CD diagrams).

## Repository Structure

```
KDO-Optimizer/
│
├── Prime_Version/                     # Main algorithm implementations
│   ├── 2017/                          # CEC2017 (D = 10, 30, 50, 100)
│   │   ├── input_data/                # CEC2017 input files
│   │   ├── cec17_func.cpp / .mexw64   # CEC2017 function source and MEX
│   │   ├── main.m                     # Main execution script
│   │   ├── KDO.m, DE.m, L_SHADE.m, …  # All algorithms
│   │   └── CEC2017_Analysis_D*/       # Generated output folders
│   │
│   ├── 2017 - 50/                     # CEC2017 at D = 50 (same structure)
│   ├── 2022/                          # CEC2022 (D = 20)
│   │   ├── input_data/
│   │   ├── cec22_test_func.cpp/.mexw64
│   │   └── main.m
│   │
│   └── Electronics/                   # Hard electronics problems
│       ├── Electro_E*_Hard.m          # Objective functions
│       ├── MAIN_Electronics_Hard.m    # Main electronics script
│       └── Results_Electronics_Hard/  # Output folders
│
├── Ablation_Study/                    # Ablation variants (NoNirvana, NoMemory, NoDharma)
│   ├── 2017/
│   ├── 2022/
│   └── Electronics-New/
│
├── scripts/                           # Optional helper scripts
└── docs/                              # Documentation
```

## Requirements

- **MATLAB R2017b or later
- **Toolboxes:**  
  - `Optimization Toolbox` (for some CEC functions)  
  - `Statistics and Machine Learning Toolbox` (for statistical analyses)  
  - `Parallel Computing Toolbox` (optional, for parallel execution)

> On non‑Windows systems, recompile the MEX files (`cec17_func.mexw64`, `cec22_test_func.mexw64`) using `mex` in MATLAB.

## Running CEC Benchmarks

### 1. Configure `main.m`

The main script is located in `Prime_Version/2017/` or `Prime_Version/2022/`. Adjust the following parameters before running:

%% ==================== CONFIGURATION ====================
benchmark_choice = 'CEC2017';   % or 'CEC2022'
dim = 30;                       % CEC2017: 10, 30, 50, 100 | CEC2022: 10, 20
num_runs = 30;                  % number of independent runs (recommended: 30)
enable_ablation = true;         % true to enable ablation study
enable_parallel = false;        % true if Parallel Toolbox is available
output_format = 'pdf';          % 'pdf', 'eps', or 'both'

% Enable/disable specific outputs
generate_convergence_curves = true;
generate_boxplots = true;
generate_performance_profiles = true;
generate_cd_diagrams = true;
generate_ablation_results = true;
```

### 2. Run the script

matlab
>> main
```

### 3. Output structure

After execution, a folder named `CEC2017_Analysis_D30_yyyy-mm-dd_HH-MM-SS/` is created with the following contents:

| Folder / File | Description |
|:---|:---|
| `Convergence_Curves/` | Convergence plots (PNG + PDF) for each function |
| `Boxplots/` | Boxplot of final results for each function |
| `Performance_Profiles/` | Dolan–Moré performance profiles |
| `CD_Diagrams/` | Critical difference diagrams (Nemenyi) |
| `Tables/LaTeX/` | Full LaTeX tables (mean, std, Friedman ranks, Wilcoxon results) |
| `Summary/` | Plain text summary statistics |
| `Ablation_Results/` | Ablation study tables (if enabled) |
| `Raw_Data/` | MAT file with raw data |

## Running Electronics Problems (Hard Versions)

### 1. Main script

The script is located at `Prime_Version/Electronics/MAIN_Electronics_Hard.m`. Key settings:

matlab
num_runs = 51;                  % runs (increased for statistical power)
population_size = 30;
MaxFEs_per_dim = 10000;
output_format = 'pdf';          % 'pdf', 'eps', or 'both'
save_excel = true;
save_latex = true;
save_png = true;
save_figures = true;

### 2. Run

matlab
>> MAIN_Electronics_Hard

### 3. Outputs

Results are saved in `Results_Electronics_Hard/` with subfolders similar to the CEC benchmark:

- `Convergence_Curves/`
- `Boxplots/`
- `Performance_Profiles/`
- `CD_Diagrams/`
- `Excel_Results/` (complete Excel workbook with all tables)
- `Latex_Tables/` (all LaTeX tables)
- `Summary_Statistics/`
- `Raw_Data/`

## Ablation Study

To enable the ablation study in the CEC benchmark, set `enable_ablation = true` in `main.m`. This runs the following KDO variants:

- `KDO_NoNirvana` – without the Nirvana reset mechanism  
- `KDO_NoMemory`  – without Cosmic Memory  
- `KDO_NoDharma`  – without the Dharma phase  

Results are saved as LaTeX tables and text summaries in the `Ablation_Results/` folder.

For the electronics problems, the ablation study is enabled by default in `MAIN_Electronics_Hard.m`, and the results are included in the corresponding LaTeX tables.

## Changing the Dimension

- **CEC2017:** Supported dimensions are 10, 30, 50, and 100. Set `dim` accordingly.
- **CEC2022:** Maximum dimension is 20 (use `dim = 20`).

> For high dimensions (50 and 100), execution time increases significantly. Consider reducing `num_runs` to 20 or enabling parallel execution.

## Parallel Execution

If you have the Parallel Computing Toolbox, set `enable_parallel = true` to use `parfor` and speed up the benchmarks.

## Recompiling MEX Files (Non‑Windows Systems)

On Linux or macOS, the provided `.mexw64` files will not work. Recompile them with:

matlab
mex -setup C++
mex cec17_func.cpp
mex cec22_test_func.cpp

## Quick Test

To quickly test CEC2017 at D = 10 with only 10 runs:

matlab
% In main.m, set:
benchmark_choice = 'CEC2017';
dim = 10;
num_runs = 10;
enable_ablation = false;
enable_parallel = false;

Then run:

matlab
>> main

## Citation


## License

This project is released under the **MIT License** – see the `LICENSE` file for details.

## Contact

- **Corresponding Author:** Dr. Seyyed Mohammad Razavi  
  Email: [smrazavi@birjand.ac.ir](mailto:smrazavi@birjand.ac.ir)  
- **GitHub:** [hosseinmoodi1/KDO-Optimizer](https://github.com/hosseinmoodi1/KDO-Optimizer)

**Happy Optimizing!**

