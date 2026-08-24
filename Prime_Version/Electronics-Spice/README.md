# Electronics-Spice: ngspice-in-the-Loop Optimization Framework

This repository contains the **MATLAB–ngspice-in-the-loop** evaluation framework for the Karma-Dharma Optimizer (KDO) and other metaheuristic algorithms, applied to four analog/RF and VLSI circuit design problems.

## 🔧 Requirements

- **MATLAB** (R2017b or later)
- **ngspice** (v38 or later) – [Download](https://ngspice.sourceforge.io/)
- **Model Libraries** (not included due to licensing):
  - `ptm180.lib` – 180nm CMOS (for Op-Amp)
  - `ptm65.lib` – 65nm CMOS (for LNA & VCO)
  - `ptm22hp.lib` – 22nm FinFET (for SRAM)

> ⚠️ The model libraries are **not** included. You must obtain them from your foundry or use the open-source PTM models from [PTM website](https://ptm.asu.edu/).

## 📂 Directory Structure
Electronics-Spice/
├── README.md
├── MAIN_Electronics_Hard.m # Main script
├── config.m # Global configuration
├── paths_local.m # Local path definitions
├── .m # All algorithm implementations (KDO, DE, GWO, ...)
├── Electro_E_Hard.m # Circuit cost functions (E1–E4)
├── extract_params_E*.m # Physical parameter extractors
├── Templates/ # SPICE netlist templates
│ ├── template_E1.sp
│ ├── template_E2.sp
│ ├── template_E3.sp
│ └── template_E4.sp
└── Models/ # Model libraries (user-provided)

## 🚀 Quick Start

### 1. Install ngspice
- **Windows**: Download the binary and update `ngspice_exe` path in `paths_local.m`.
- **Linux/macOS**: Install via package manager (`sudo apt install ngspice`).

### 2. Set up model libraries
Place `ptm180.lib`, `ptm65.lib`, and `ptm22hp.lib` in the `Models/` folder.

### 3. Configure paths
Edit `paths_local.m` to point to your `ngspice` executable:

```matlab
if ispc
    paths.ngspice_exe = 'C:\path\to\ngspice_con.exe';
else
    paths.ngspice_exe = 'ngspice';
end

4. Run the benchmark
Open MATLAB, navigate to this folder, and run:

MAIN_Electronics_Hard

This will execute the optimization loop for all circuits and algorithms, generate statistical analyses, and save results in Results_Electronics_Hard/.

📊 Evaluated Circuits
ID	Circuit	Technology	Variables	Key Constraints
E1	Two-Stage Miller Op-Amp	180nm CMOS	18	Gain ≥ 74dB, GBW ≥ 5MHz, PM ≥ 60°, SR ≥ 10V/µs, Power ≤ 5mW
E2	RF CMOS Cascode LNA	65nm CMOS	12	Gain ≥ 15dB, NF ≤ 2.5dB, S11 ≤ -10dB, Yield ≥ 99%, Power ≤ 20mW
E3	LC-VCO with PVT	65nm CMOS	11	PN ≤ -115dBc/Hz, Power ≤ 8mW, TR ≥ 10%, FOM ≤ -180dBc/Hz
E4	6T FinFET SRAM	22nm FinFET	9	SNM ≥ 100mV (3σ), Yield ≥ 99%, CR ≥ 1.2, PR ≤ 0.8
🤖 Algorithms
Algorithm	Type	Reference
KDO	Proposed	Karma-Dharma Optimizer (this work)
DE	Differential Evolution	Storn & Price (1997)
L-SHADE	Adaptive DE	Tanabe & Fukunaga (2014)
CMA-ES	Evolution Strategy	Hansen & Ostermeier (2001)
GBO	Gradient-Based	Ahmadianfar et al. (2020)
RUN	Runge-Kutta	Ahmadianfar et al. (2021)
GWO	Grey Wolf	Mirjalili et al. (2014)
WOA	Whale	Mirjalili & Lewis (2016)
HHO	Harris Hawks	Heidari et al. (2019)
AVOA	African Vultures	Abdollahzadeh et al. (2021)
COA	Coyote	Pierezan & Coelho (2021)
📈 Outputs
After execution, the following results are generated in Results_Electronics_Hard/timestamp/:

Folder	Contents
Convergence_Curves/	Convergence plots (PNG + PDF)
Boxplots/	Box plots of fitness distributions
CD_Diagrams/	Critical Difference diagrams
Performance_Profiles/	Dolan–Moré performance profiles
Excel_Results/	Full results in Excel (.xlsx)
Latex_Tables/	LaTeX tables for publications
Raw_Data/	MAT files and physical parameters
⚙️ Configuration
In MAIN_Electronics_Hard.m, you can adjust:
num_runs = 5;            % Number of independent runs (fast testing)
MaxFEs_per_dim = 100;    % Max function evaluations per dimension
population_size = 10;    % Population size per algorithm
For publication-quality results, set num_runs = 51 (analytical) or num_runs = 5 (ngspice) as per the paper.

🔬 SPICE Evaluation Details
Physics-driven: Each candidate solution is compiled into a SPICE netlist and simulated in batch mode.

Direct vs. derived metrics: Some metrics (e.g., gain, power) are parsed directly from SPICE logs; others (e.g., slew rate, SNM) are derived analytically.

Failure handling: Simulation failures are trapped and penalized with distance-dependent penalties.

📖 Citation
If you use this code, please cite our paper:


📬 Contact
Dr. Seyyed Mohammad Razavi
Email: smrazavi@birjand.ac.ir
Department of Electronics, University of Birjand, Iran

📝 License
This code is provided for research purposes. Please cite the paper when using it in your work.