%% MAIN_Electronics.m - Comprehensive Metaheuristic Algorithm Comparison for Electronics Engineering Problems
% Version: 7.1 - FAST TRACK WITH PHYSICAL SPACE FIX
% Author: Hossein Moodi
% Date: 2026-08-21
% Description: This script compares various metaheuristic algorithms on 
%              electronics engineering optimization problems with full
%              statistical analysis (Paired Wilcoxon Signed-Rank, Cliff's Delta,
%              Holm-Bonferroni Correction, Friedman, CD Diagram, 
%              Performance Profiles) as required for Q1 publications.
%              NEW: Extracts physical parameters for SPICE validation.
%              FAST VERSION: num_runs=5, MaxFEs_per_dim=100
%              FIX: KDO now works in physical space (lb, ub) directly.

%% Clear workspace and close figures
clear all;
close all;
clc;
warning('off', 'all');

% Start timer for total execution time
total_timer = tic;

%% User Configuration
% =========================================================================
benchmark_choice = 'Electronics';

% *** Number of independent runs  ***
num_runs = 15;

% Algorithm parameters
population_size = 10;      % Increased for better exploration
% *** MaxFEs per dimension 
MaxFEs_per_dim = 100;

% Analysis levels
analysis_levels = {'Full', 'Partial', 'Behavioral'};

% Results saving options
save_excel   = true;
save_latex   = true;
save_png     = true;      % PNG for preview
save_mat     = true;
save_figures = true;

% Output format for journal submission (Nature/Scientific Reports)
% Options: 'pdf' (recommended), 'eps', 'both'
output_format = 'pdf';     % <-- NEW: PDF format for journal

% Statistical significance level
alpha = 0.05;

%% Define Electronics Problems (HARD VERSIONS)
% =========================================================================
problem_names = {
    'E1: Two-Stage Miller Op-Amp (18D) - HARD';
    'E2: RF CMOS Cascode LNA 2.4GHz (12D) - HARD with Noise';
    'E3: Low Phase-Noise LC-VCO with PVT (11D) - HARD';
    'E4: Robust 6T FinFET SRAM Cell (9D) - HARD';
};

problem_descriptions = {
    'Sizing of two-stage op-amp (18 design variables)';
    '2.4 GHz cascode LNA optimization with process variation, yield requirement 99%';
    'LC-VCO phase-noise minimization with full PVT variation, tighter PN constraint';
    'FinFET SRAM cell robust optimization with Nmc=200, yield requirement 99%';
};
% Problem IDs for different analysis levels
problem_ids.Full       = [1, 2, 3, 4];
problem_ids.Partial    = [1, 4];        % OpAmp + SRAM
problem_ids.Behavioral = [2, 3];        % LNA + VCO


%% Define Algorithms
% =========================================================================
algorithms = {
   
    %'CMA_ES',       @CMA_ES,       'Covariance Matrix Adaptation ES';
    %'AVOA',         @AVOA,         'African Vultures Optimization Algorithm';
    %'COA',          @COA,          'Crayfish Optimization Algorithm';
    %'RUN',          @RUN,          'Runge-Kutta Optimizer';
    %'L_SHADE',      @L_SHADE,      'Success-History DE';
    %'GBO',          @GBO,          'Gradient-Based Optimizer';
    'KDO',          @KDO,          'Karma Dharma Optimizer (Proposed)';
    %'DE',           @DE,           'Differential Evolution (F=0.8, CR=0.9)';    
    'GWO',          @GWO,          'Grey Wolf Optimizer';
    %'WOA',          @WOA,          'Whale Optimization Algorithm';
    %'HHO',          @HHO,          'Harris Hawks Optimization';
    
    
};
num_algorithms = size(algorithms, 1);

%% Setup Results Directory
% =========================================================================
timestamp      = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
results_folder = fullfile('Results_Electronics_Hard', timestamp);

subfolders = {
    'Convergence_Curves';
    'Boxplots';
    'Performance_Profiles';
    'CD_Diagrams';
    'Excel_Results';
    'Latex_Tables';
    'Summary_Statistics';
    'Raw_Data';
    'Figures';
};

for i = 1:length(subfolders)
    if ~exist(fullfile(results_folder, subfolders{i}), 'dir')
        mkdir(fullfile(results_folder, subfolders{i}));
    end
end

% Start diary for logging
diary_file = fullfile(results_folder, 'execution_log.txt');
diary(diary_file);
fprintf('Execution started at: %s\n', datestr(now));
fprintf('Results folder: %s\n', results_folder);
fprintf('Output format for journal: %s\n', upper(output_format));
fprintf('VERSION: 7.1 - FAST TRACK WITH PHYSICAL SPACE FIX\n\n');

%% Display Configuration
% =========================================================================
fprintf('\n========================================\n');
fprintf('CMOS ELECTRONICS OPTIMIZATION BENCHMARK SUITE - HARD VERSION (FAST)\n');
fprintf('========================================\n\n');

fprintf('Configuration:\n');
fprintf('  Benchmark type: %s (HARD VERSION)\n', benchmark_choice);
fprintf('  Number of runs: %d (Reduced for fast testing)\n', num_runs);
fprintf('  Population size: %d\n', population_size);
fprintf('  Max FEs per dimension: %d (Reduced from 500)\n', MaxFEs_per_dim);
fprintf('  Statistical significance level: alpha = %.2f\n', alpha);
fprintf('  Output format: %s (PNG + %s)\n', upper(output_format), upper(output_format));
fprintf('\n');

fprintf('Electronics Problems (HARD VERSIONS):\n');
for i = 1:length(problem_names)
    fprintf('  %s\n', problem_names{i});
    fprintf('    %s\n', problem_descriptions{i});
end
fprintf('\n');

fprintf('Algorithms (%d total):\n', num_algorithms);
for i = 1:num_algorithms
    fprintf('  %2d. %-10s - %s\n', i, algorithms{i,1}, algorithms{i,3});
end
fprintf('\n');

%% Initialize Results Storage
% =========================================================================
fprintf('Initializing results storage...\n');

current_level    = analysis_levels{1};
current_problems = problem_ids.(current_level);
num_problems     = length(current_problems);

% Main result arrays
best_values        = zeros(num_algorithms, num_problems, num_runs);
convergence_curves = cell(num_algorithms, num_problems, num_runs);
execution_times    = zeros(num_algorithms, num_problems, num_runs);
auc_values         = zeros(num_algorithms, num_problems, num_runs);

% NEW: Store best positions for physical parameter extraction
best_positions     = cell(num_algorithms, num_problems, num_runs);

% Statistical summary arrays
median_values = zeros(num_algorithms, num_problems);
mean_values   = zeros(num_algorithms, num_problems);
std_values    = zeros(num_algorithms, num_problems);
min_values    = zeros(num_algorithms, num_problems);
max_values    = zeros(num_algorithms, num_problems);

% Ranking arrays
ranks          = zeros(num_algorithms, num_problems);
friedman_ranks = zeros(num_algorithms, num_problems);

%% Main Optimization Loop
% =========================================================================
fprintf('\n========================================\n');
fprintf('STARTING OPTIMIZATION PROCESS (HARD PROBLEMS - FAST)\n');
fprintf('========================================\n\n');

total_iterations   = num_algorithms * num_problems * num_runs;
current_iteration  = 0;

for algo_idx = 1:num_algorithms
    algo_name   = algorithms{algo_idx, 1};
    algo_handle = algorithms{algo_idx, 2};

    fprintf('\n--- Algorithm: %s ---\n', algo_name);

    for prob_idx = 1:num_problems
        problem_id   = current_problems(prob_idx);
        problem_name = problem_names{problem_id};

        % Set problem-specific parameters (HARD VERSIONS) - CORRECTED BOUNDS
        switch problem_id
            case 1   % E1: Two-Stage Miller Op-Amp (18D)
                fobj = @Electro_E1_Hard;
                dim = 18;
                lb = [1e-6 * ones(1,8), 200e-9 * ones(1,8), 1e-6, 0.1e-12];
                ub = [100e-6 * ones(1,8), 1e-6 * ones(1,8), 500e-6, 10e-12];
            case 2   % E2: RF CMOS Cascode LNA (12D) - Hard VERSION
                fobj = @Electro_E2_Hard;
                dim  = 12;
                lb = [0.5e-6, 0.5e-6, 0.5e-6, 0.5e-6, 0.1e-9, 0.1e-9, 0.1e-9, ...
                      10e-15, 10e-15, 10e-15, 0.1e-3, 0.5];
                ub = [50e-6, 50e-6, 50e-6, 50e-6, 10e-9, 10e-9, 10e-9, ...
                      500e-15, 500e-15, 500e-15, 5e-3, 1.8];

            case 3   % E3: Low Phase-Noise LC-VCO with PVT (11D) - Hard VERSION
                fobj = @Electro_E3_Hard;
                dim  = 11;
                lb = [0.5e-9, 50e-15, 10e-15, 1e-6, 1e-6, 0.5e-3, 100, 5e-15, 0, 0.1, 50e6];
                ub = [20e-9, 2e-12, 500e-15, 15e-6, 25e-6, 10e-3, 10e3, 200e-15, 1.2, 1.2, 10e9];
            case 4   % E4: Robust 6T FinFET SRAM (9D) - Hard VERSION
                fobj = @Electro_E4_Hard;
                dim  = 9;
                lb = [1, 1, 1, 1, 14e-9, 14e-9, 14e-9, 0.6, 0.8e-9];
                ub = [8, 8, 8, 8, 40e-9, 40e-9, 40e-9, 1.0, 1.5e-9];
        end

        MaxFEs         = MaxFEs_per_dim * dim;
        max_iterations = ceil(MaxFEs / population_size);

        fprintf('  Problem: %s (dim=%d, MaxFEs=%d, MaxIter=%d)\n', ...
            problem_name, dim, MaxFEs, max_iterations);
        fprintf('    Run: ');

        for run_idx = 1:num_runs
            current_iteration = current_iteration + 1;

            if mod(run_idx, 5) == 0
                fprintf('%d ', run_idx);
            end

            % Reproducible seed
            rng(run_idx * 1000 + algo_idx * 100 + prob_idx * 10);

            try
                t_start = tic;
                [Best_score, Best_pos, Convergence_curve] = algo_handle( ...
                    population_size, max_iterations, lb, ub, dim, fobj, MaxFEs);
                exec_time = toc(t_start);

                best_values(algo_idx, prob_idx, run_idx)       = Best_score;
                convergence_curves{algo_idx, prob_idx, run_idx} = Convergence_curve;
                execution_times(algo_idx, prob_idx, run_idx)    = exec_time;
                best_positions{algo_idx, prob_idx, run_idx}     = Best_pos;

                % Normalized AUC
                if ~isempty(Convergence_curve)
                    cv_min = min(Convergence_curve);
                    cv_max = max(Convergence_curve);
                    if cv_max > cv_min
                        norm_curve = (Convergence_curve - cv_min) / (cv_max - cv_min + eps);
                    else
                        norm_curve = zeros(size(Convergence_curve));
                    end
                    auc_values(algo_idx, prob_idx, run_idx) = trapz( ...
                        linspace(0, 1, length(norm_curve)), norm_curve);
                end

            catch ME
                fprintf('\n    Error in %s - %s - Run %d: %s\n', ...
                    algo_name, problem_name, run_idx, ME.message);
                best_values(algo_idx, prob_idx, run_idx)       = inf;
                convergence_curves{algo_idx, prob_idx, run_idx} = [];
                execution_times(algo_idx, prob_idx, run_idx)    = inf;
                auc_values(algo_idx, prob_idx, run_idx)         = 0;
                best_positions{algo_idx, prob_idx, run_idx}     = [];
            end
        end

        fprintf('\n');

        % Per-combination statistics
        valid_vals = squeeze(best_values(algo_idx, prob_idx, :));
        valid_vals = valid_vals(isfinite(valid_vals));

        if ~isempty(valid_vals)
            median_values(algo_idx, prob_idx) = median(valid_vals);
            mean_values(algo_idx, prob_idx)   = mean(valid_vals);
            std_values(algo_idx, prob_idx)    = std(valid_vals);
            min_values(algo_idx, prob_idx)    = min(valid_vals);
            max_values(algo_idx, prob_idx)    = max(valid_vals);

            fprintf('    Results: Min=%.4e, Median=%.4e, Mean=%.4e, Std=%.4e\n', ...
                min_values(algo_idx, prob_idx), ...
                median_values(algo_idx, prob_idx), ...
                mean_values(algo_idx, prob_idx), ...
                std_values(algo_idx, prob_idx));
        else
            fprintf('    Warning: All runs failed for this combination\n');
        end
    end  % end of algo_idx loop
end  % end of prob_idx loop

fprintf('\n========================================\n');
fprintf('OPTIMIZATION COMPLETED\n');
fprintf('========================================\n\n');

%% ==================== EXTRACT PHYSICAL PARAMETERS ====================
% This section extracts physical circuit parameters from the best solutions
% for SPICE validation and comparison with MATLAB models.
% =========================================================================
fprintf('\n========================================\n');
fprintf('EXTRACTING PHYSICAL PARAMETERS FOR ALL CIRCUITS\n');
fprintf('========================================\n\n');

% Define parameter extraction functions for each problem
% NOTE: These functions (extract_params_E1.m, etc.) must be in the MATLAB path
% If they are not present, the extraction will fail with a clear error message.
% Required files:
%   extract_params_E1.m - extracts Gain, GBW, PM, SR, Power, ICMR, Output Swing from Op-Amp solution
%   extract_params_E2.m - extracts Gain, NF, S11, S22, Power, Yield from LNA solution
%   extract_params_E3.m - extracts f_osc, PN, Power, TR, FOM from LC-VCO solution
%   extract_params_E4.m - extracts SNM_mean, Yield, I_leak, Area, CR, PR from SRAM solution

extract_funcs = {
    @extract_params_E1;   % Op-Amp
    @extract_params_E2;   % LNA
    @extract_params_E3;   % LC-VCO
    @extract_params_E4;   % SRAM
};

% Define parameter names and labels for each problem
param_defs = cell(4, 1);

% E1: Op-Amp
param_defs{1}.names = {'Gain_dB', 'GBW_Hz', 'PM_deg', 'SR_Vus', 'Power_W', 'ICMR_V', 'Output_swing_V'};
param_defs{1}.labels = {'Gain (dB)', 'GBW (MHz)', 'PM (deg)', 'SR (V/us)', 'Power (mW)', 'ICMR (V)', 'Output Swing (V)'};
param_defs{1}.scales = [1, 1e-6, 1, 1, 1e3, 1, 1];
param_defs{1}.format = '%.2f';

% E2: LNA
param_defs{2}.names = {'Gain_dB', 'NF_dB', 'S11_dB', 'S22_dB', 'Power_W', 'Yield'};
param_defs{2}.labels = {'Gain (dB)', 'NF (dB)', 'S11 (dB)', 'S22 (dB)', 'Power (mW)', 'Yield (%)'};
param_defs{2}.scales = [1, 1, 1, 1, 1e3, 100];
param_defs{2}.format = '%.2f';

% E3: LC-VCO
param_defs{3}.names = {'f_osc_Hz', 'PN_dBc', 'Power_W', 'TR_percent', 'FOM_dB'};
param_defs{3}.labels = {'f_osc (GHz)', 'PN (dBc/Hz)', 'Power (mW)', 'TR (%)', 'FOM (dB)'};
param_defs{3}.scales = [1e-9, 1, 1e3, 1, 1];
param_defs{3}.format = '%.2f';

% E4: SRAM
param_defs{4}.names = {'SNM_mean_V', 'Yield', 'I_leak_A', 'Area_um2', 'CR', 'PR'};
param_defs{4}.labels = {'SNM (mV)', 'Yield (%)', 'I_leak (nA)', 'Area (um^2)', 'CR', 'PR'};
param_defs{4}.scales = [1e3, 100, 1e9, 1, 1, 1];
param_defs{4}.format = '%.2f';

% Storage for physical parameters
physical_results = cell(num_algorithms, num_problems);

for prob_idx = 1:num_problems
    problem_id = current_problems(prob_idx);
    fprintf('Processing Problem P%d (%s)...\n', prob_idx, problem_names{problem_id});
    
    % Select the appropriate extraction function
    extract_func = extract_funcs{problem_id};
    pdef = param_defs{problem_id};
    
    % Store extracted parameters for each algorithm
    param_mean = zeros(num_algorithms, length(pdef.names));
    param_std  = zeros(num_algorithms, length(pdef.names));
    
    for algo_idx = 1:num_algorithms
        % Collect all valid runs
        param_vals = [];
        valid_runs = 0;
        failed_runs = 0;
        
        for run_idx = 1:num_runs
            best_pos = best_positions{algo_idx, prob_idx, run_idx};
            if ~isempty(best_pos) && all(isfinite(best_pos))
                try
                    p = extract_func(best_pos);
                    % Convert struct to numeric row vector
                    p_vec = zeros(1, length(pdef.names));
                    for k = 1:length(pdef.names)
                        p_vec(k) = p.(pdef.names{k});
                    end
                    param_vals = [param_vals; p_vec];
                    valid_runs = valid_runs + 1;
                catch ME
                    failed_runs = failed_runs + 1;
                    % Log the error with full context for debugging
                    fprintf('WARNING: Extraction failed for Algorithm=%s, Problem=P%d, Run=%d\n', ...
                        algorithms{algo_idx,1}, problem_id, run_idx);
                    fprintf('  Error: %s\n', ME.message);
                    % Continue with next run (do NOT silently skip)
                end
            else
                failed_runs = failed_runs + 1;
                fprintf('WARNING: Empty or invalid best_pos for Algorithm=%s, Problem=P%d, Run=%d\n', ...
                    algorithms{algo_idx,1}, problem_id, run_idx);
            end
        end
        
        if ~isempty(param_vals)
            param_mean(algo_idx, :) = mean(param_vals, 1);
            param_std(algo_idx, :)  = std(param_vals, 0, 1);
            fprintf('  Algorithm %s: %d valid runs, %d failed runs\n', ...
                algorithms{algo_idx,1}, valid_runs, failed_runs);
        else
            param_mean(algo_idx, :) = NaN;
            param_std(algo_idx, :)  = NaN;
            fprintf('  Algorithm %s: NO valid runs (all %d failed)\n', ...
                algorithms{algo_idx,1}, failed_runs);
        end
        
        physical_results{algo_idx, prob_idx} = struct(...
            'mean', param_mean(algo_idx, :), ...
            'std',  param_std(algo_idx, :));
    end
    
    % Display table for this problem
    fprintf('\n--- %s ---\n', problem_names{problem_id});
    fprintf('%-12s', 'Algorithm');
    for k = 1:length(pdef.labels)
        fprintf('  %-14s', pdef.labels{k});
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 12 + 16*length(pdef.labels)));
    
    % Find best for each parameter (optional)
    for algo_idx = 1:num_algorithms
        if ~isnan(param_mean(algo_idx, 1))
            % Bold if KDO
            if strcmp(algorithms{algo_idx,1}, 'KDO')
                fprintf('\\textbf{%-10s}', algorithms{algo_idx,1});
            else
                fprintf('%-12s', algorithms{algo_idx,1});
            end
            for k = 1:length(pdef.names)
                val = param_mean(algo_idx, k) * pdef.scales(k);
                fprintf('  %14.4f', val);
            end
            fprintf('\n');
        end
    end
    fprintf('\n');
end

% Save physical parameters to a separate MAT file
phys_file = fullfile(results_folder, 'Raw_Data', 'Physical_Parameters.mat');
save(phys_file, 'physical_results', 'param_defs', 'problem_names', 'algorithms');
fprintf('Physical parameters saved to: %s\n', phys_file);

% ==================== GENERATE SPICE COMPARISON TABLE ====================
% This table compares MATLAB model predictions with SPICE simulations
% (to be filled after SPICE simulations are run)
% =========================================================================
fprintf('\n========================================\n');
fprintf('SPICE VALIDATION TABLE (TO BE COMPLETED)\n');
fprintf('========================================\n\n');

fprintf('The following table compares MATLAB model predictions with SPICE simulations.\n');
fprintf('Fill the SPICE columns after running HSPICE simulations.\n\n');

% Find KDO and L-SHADE indices
kdo_idx = find(strcmp(algorithms(:,1), 'KDO'));
lshade_idx = find(strcmp(algorithms(:,1), 'L_SHADE'));
de_idx = find(strcmp(algorithms(:,1), 'DE'));

if ~isempty(kdo_idx) && ~isempty(lshade_idx)
    % For Op-Amp (Problem 1)
    prob_idx = 1;
    problem_id = current_problems(prob_idx);
    pdef = param_defs{problem_id};
    
    fprintf('=== SPICE VALIDATION FOR %s ===\n', problem_names{problem_id});
    fprintf('%-20s', 'Parameter');
    fprintf('  %-15s', 'KDO (MATLAB)');
    fprintf('  %-15s', 'KDO (SPICE)');
    fprintf('  %-15s', 'Error (%)');
    fprintf('  %-15s\n', 'L-SHADE (SPICE)');
    fprintf('%s\n', repmat('-', 1, 100));
    
    % Extract Op-Amp parameters
    kdo_matlab = physical_results{kdo_idx, prob_idx}.mean;
    lshade_spice = physical_results{lshade_idx, prob_idx}.mean;  % Will be replaced with SPICE data
    
    % !!! IMPORTANT: Replace this section with actual SPICE results !!!
    % For now, using NaN to indicate missing data (prevents accidental use of dummy data)
    kdo_spice = NaN(size(kdo_matlab));
    fprintf('WARNING: SPICE validation data is placeholder. Replace with actual results.\n');
    
    % Display table
    for k = 1:length(pdef.names)
        param_label = pdef.labels{k};
        scale = pdef.scales(k);
        matlab_val = kdo_matlab(k) * scale;
        spice_val = kdo_spice(k) * scale;
        error_pct = abs(matlab_val - spice_val) / max(abs(matlab_val), 1e-9) * 100;
        lshade_val = lshade_spice(k) * scale;
        
        if isnan(spice_val)
            fprintf('%-20s', param_label);
            fprintf('  %15.4f', matlab_val);
            fprintf('  %15s', 'N/A');
            fprintf('  %14s', 'N/A');
            fprintf('  %15.4f\n', lshade_val);
        else
            fprintf('%-20s', param_label);
            fprintf('  %15.4f', matlab_val);
            fprintf('  %15.4f', spice_val);
            fprintf('  %14.2f%%', error_pct);
            fprintf('  %15.4f\n', lshade_val);
        end
    end
end

fprintf('\nNOTE: The SPICE column currently shows NaN to indicate missing data.\n');
fprintf('Replace these with actual HSPICE simulation results for final validation.\n\n');

%% Calculate Basic Rankings (Median-based)
% =========================================================================
fprintf('Calculating basic rankings (median-based)...\n');

for prob_idx = 1:num_problems
    prob_medians = median_values(:, prob_idx);
    [sorted_vals, sort_idx] = sort(prob_medians);

    current_rank = 1;
    ranks(sort_idx(1), prob_idx) = 1;

    for i = 2:num_algorithms
        if abs(sorted_vals(i) - sorted_vals(i-1)) < 1e-10
            ranks(sort_idx(i), prob_idx) = current_rank;
        else
            current_rank = i;
            ranks(sort_idx(i), prob_idx) = current_rank;
        end
    end
end

avg_ranks = mean(ranks, 2);

%% ==================== ADVANCED PAIRED STATISTICAL ANALYSIS ====================
% *** UPGRADED: Replaces old Wilcoxon section with paired tests, effect size, and Holm correction ***
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('ADVANCED PAIRED STATISTICAL ANALYSIS');
fprintf('\n%s\n', repmat('=', 1, 80));

% Find proposed algorithm index
proposed_idx = find(strcmp(algorithms(:,1), 'KDO'), 1);
if isempty(proposed_idx), proposed_idx = 1; end

% Pre-allocate storage for new stats
P_raw       = zeros(num_algorithms-1, num_problems); % Excluding KDO
P_holm      = zeros(num_algorithms-1, num_problems);
EffectSize  = zeros(num_algorithms-1, num_problems);

fprintf('\nComputing paired tests, effect sizes, and Holm correction vs %s...\n', algorithms{proposed_idx,1});

% Define other_indices outside the problem loop
other_indices = setdiff(1:num_algorithms, proposed_idx);

% Accumulators for overall W/T/L
overall_wins = 0;
overall_ties = 0;
overall_losses = 0;

for prob_idx = 1:num_problems
    fprintf('  Problem %d: %s\n', prob_idx, problem_names{current_problems(prob_idx)});
    kdo_results = squeeze(best_values(proposed_idx, prob_idx, :));
    
    raw_p_vec = zeros(length(other_indices), 1);
    
    for a = 1:length(other_indices)
        comp_idx = other_indices(a);
        comp_results = squeeze(best_values(comp_idx, prob_idx, :));
        
        valid_mask = isfinite(kdo_results) & isfinite(comp_results);
        kdo_valid = kdo_results(valid_mask);
        comp_valid = comp_results(valid_mask);
        
        if length(kdo_valid) >= 5
            % --- 1. PAIRED Wilcoxon Signed-Rank Test (NOT Rank-Sum!) ---
            [p, ~] = signrank(kdo_valid, comp_valid, 'alpha', alpha, 'method', 'approximate');
            P_raw(a, prob_idx) = p;
            
            % --- 2. Effect Size: Cliff's Delta ---
            EffectSize(a, prob_idx) = cliffs_delta(kdo_valid, comp_valid);
        else
            P_raw(a, prob_idx) = NaN;
            EffectSize(a, prob_idx) = NaN;
        end
    end
    
    % --- 3. Holm-Bonferroni Correction per problem ---
    p_holm_vec = holm_bonferroni(P_raw(~isnan(P_raw(:, prob_idx)), prob_idx));
    P_holm(~isnan(P_raw(:, prob_idx)), prob_idx) = p_holm_vec;
    
    % --- Accumulate W/T/L for this problem ---
    for a = 1:length(other_indices)
        comp_idx = other_indices(a);
        p_adj = P_holm(a, prob_idx);
        if isnan(p_adj)
            continue;
        elseif p_adj >= alpha
            overall_ties = overall_ties + 1;
        else
            kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
            comp_mean = mean(squeeze(best_values(comp_idx, prob_idx, :)), 'omitnan');
            if kdo_mean < comp_mean
                overall_wins = overall_wins + 1;
            else
                overall_losses = overall_losses + 1;
            end
        end
    end
end

% --- 4. Overall Summary Statistics (Win/Tie/Loss based on Holm-corrected p) ---
fprintf('\n--- Summary of Paired Wilcoxon Signed-Rank Test vs KDO ---\n');
fprintf('Overall KDO performance: %d Wins / %d Ties / %d Losses (Holm-corrected)\n', ...
    overall_wins, overall_ties, overall_losses);

% Also print per-algorithm summary
fprintf('\nDetailed per-algorithm summary:\n');
for a = 1:length(other_indices)
    comp_idx = other_indices(a);
    w = 0; t = 0; l = 0;
    for prob_idx = 1:num_problems
        p_adj = P_holm(a, prob_idx);
        if isnan(p_adj)
            continue;
        elseif p_adj >= alpha
            t = t + 1;
        else
            kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
            comp_mean = mean(squeeze(best_values(comp_idx, prob_idx, :)), 'omitnan');
            if kdo_mean < comp_mean
                w = w + 1;
            else
                l = l + 1;
            end
        end
    end
    fprintf('  vs %-10s : %d W / %d T / %d L\n', algorithms{comp_idx,1}, w, t, l);
end

%% Friedman Test (Q1 Standard)
% =========================================================================
fprintf('\n========================================\n');
fprintf('FRIEDMAN TEST (Q1 Standard)\n');
fprintf('========================================\n\n');

friedman_data = zeros(num_runs, num_algorithms, num_problems);
for prob_idx = 1:num_problems
    for algo_idx = 1:num_algorithms
        friedman_data(:, algo_idx, prob_idx) = ...
            squeeze(best_values(algo_idx, prob_idx, :));
    end
end

friedman_stats = struct();
for prob_idx = 1:num_problems
    fprintf('\nProblem %d: %s\n', prob_idx, ...
        problem_names{current_problems(prob_idx)});

    prob_data = squeeze(friedman_data(:, :, prob_idx));
    % --- FIX: Remove rows that are completely NaN ---
    valid_rows = ~all(isnan(prob_data), 2);
    prob_data = prob_data(valid_rows, :);
    if size(prob_data, 1) < 2
        fprintf('  Not enough valid data for Friedman test on P%d\n', prob_idx);
        friedman_stats(prob_idx).p_value = NaN;
        friedman_stats(prob_idx).table = [];
        friedman_stats(prob_idx).ranks = NaN(1, num_algorithms);
        continue;
    end
    % -------------------------------------------------

    [p_friedman, tbl, stats] = friedman(prob_data, 1, 'off');

    friedman_stats(prob_idx).p_value = p_friedman;
    friedman_stats(prob_idx).table   = tbl;
    friedman_stats(prob_idx).ranks   = stats.meanranks;

    fprintf('  Friedman test p-value: %.6f\n', p_friedman);
    if p_friedman < alpha
        fprintf('  Significant differences exist among algorithms (p < %.2f)\n', alpha);
    else
        fprintf('  No significant differences (p >= %.2f)\n', alpha);
    end

    fprintf('  Average ranks (Friedman):\n');
    for algo_idx = 1:num_algorithms
        fprintf('    %-10s: %.3f\n', algorithms{algo_idx,1}, ...
            stats.meanranks(algo_idx));
        friedman_ranks(algo_idx, prob_idx) = stats.meanranks(algo_idx);
    end
end

% Overall Friedman test
fprintf('\n========================================\n');
fprintf('OVERALL FRIEDMAN TEST (All Problems)\n');
fprintf('========================================\n\n');

overall_data = zeros(num_runs * num_problems, num_algorithms);
for algo_idx = 1:num_algorithms
    overall_data(:, algo_idx) = reshape(best_values(algo_idx, :, :), [], 1);
end

% --- FIX: Remove rows that are completely NaN ---
valid_rows = ~all(isnan(overall_data), 2);
overall_data = overall_data(valid_rows, :);

if size(overall_data, 1) < 2
    fprintf('Not enough valid data for overall Friedman test (all runs failed).\n');
    p_overall = NaN;
    tbl_overall = [];
    stats_overall.meanranks = NaN(1, num_algorithms);
else
    [p_overall, tbl_overall, stats_overall] = friedman(overall_data, 1, 'off');
    fprintf('Overall Friedman test p-value: %.6f\n', p_overall);
    if p_overall < alpha
        fprintf('Significant differences exist among algorithms overall (p < %.2f)\n', alpha);
    else
        fprintf('No significant differences overall (p >= %.2f)\n', alpha);
    end

    fprintf('\nOverall average ranks (Friedman):\n');
    for algo_idx = 1:num_algorithms
        fprintf('  %-10s: %.3f\n', algorithms{algo_idx,1}, ...
            stats_overall.meanranks(algo_idx));
    end
end

%% Critical Difference (CD) Diagram
% =========================================================================
fprintf('\n========================================\n');
fprintf('CRITICAL DIFFERENCE DIAGRAM\n');
fprintf('========================================\n\n');

% Nemenyi post-hoc test critical difference
q_alpha_table = [0, 1.960, 2.344, 2.569, 2.728, ...
                 2.850, 2.949, 3.031, 3.102, 3.164, 3.219];

k_algs = num_algorithms;
N_cd   = num_problems;

if k_algs <= length(q_alpha_table)
    q_alpha = q_alpha_table(k_algs);
else
    q_alpha = 3.5;
end

CD = q_alpha * sqrt(k_algs * (k_algs + 1) / (6 * N_cd));
fprintf('Critical Difference (CD) = %.4f\n', CD);
fprintf('  k = %d algorithms, N = %d problems\n', k_algs, N_cd);
fprintf('  q_alpha(%.2f) = %.3f\n', alpha, q_alpha);

% Compute overall Friedman ranks per algorithm
overall_friedman_ranks = stats_overall.meanranks;

fprintf('\nOverall Friedman average ranks:\n');
for i = 1:num_algorithms
    fprintf('  %-10s: %.3f\n', algorithms{i,1}, overall_friedman_ranks(i));
end

% --- Plot CD Diagram (JOURNAL FORMAT: PNG + PDF)---
fig_cd = figure('Name', 'CD Diagram', ...
    'Position', [100 100 900 500], ...
    'Color', 'white');

ax_cd = axes('Parent', fig_cd);

% Sort by rank
[sorted_fr, sort_fr_idx] = sort(overall_friedman_ranks);
sorted_names = algorithms(sort_fr_idx, 1);

% Layout parameters
x_left   = 1;
x_right  = k_algs;
y_top    = k_algs + 2;
y_bottom = 0;

hold(ax_cd, 'on');

% Axis line
plot(ax_cd, [x_left x_right], [y_top y_top], 'k-', 'LineWidth', 2);

% Tick marks
tick_step = 0.5;
tick_vals = x_left:tick_step:x_right;
for tv = tick_vals
    if mod((tv - x_left) / tick_step, 2) == 0
        tick_len = 0.3;
        lw_t = 1.5;
    else
        tick_len = 0.15;
        lw_t = 1.0;
    end
    plot(ax_cd, [tv tv], [y_top y_top+tick_len], 'k-', 'LineWidth', lw_t);
    if mod((tv - x_left) / tick_step, 2) == 0
        text(ax_cd, tv, y_top + 0.55, sprintf('%.1f', tv), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end
end

% Color palette for algorithms
colors_cd = lines(k_algs);

% Draw algorithm markers and connecting lines
half = ceil(k_algs / 2);
for i = 1:k_algs
    algo_rank = sorted_fr(i);
    algo_clr  = colors_cd(i, :);

    if i <= half
        y_line = y_top - 1 - (i-1) * 0.8;
        x_end  = x_left - 0.2;
        h_align = 'right';
    else
        y_line = y_top - 1 - (i - half - 1) * 0.8;
        x_end  = x_right + 0.2;
        h_align = 'left';
    end

    % Vertical drop line
    plot(ax_cd, [algo_rank algo_rank], [y_top y_line + 0.05], ...
        '-', 'Color', algo_clr, 'LineWidth', 1.5);

    % Horizontal connector
    plot(ax_cd, [algo_rank x_end], [y_line y_line], ...
        '-', 'Color', algo_clr, 'LineWidth', 1.5);

    % Marker
    plot(ax_cd, algo_rank, y_top, 'o', ...
        'MarkerFaceColor', algo_clr, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8);

    % Label
    if i == find(sort_fr_idx == proposed_idx)
        label_text = sprintf('\\bf{%s} (%.2f)', sorted_names{i}, algo_rank);
    else
        label_text = sprintf(' %s (%.2f) ', sorted_names{i}, algo_rank);
    end
    
    text(ax_cd, x_end, y_line, label_text, ...
        'HorizontalAlignment', h_align, ...
        'FontSize', 9, 'Color', algo_clr, 'FontWeight', 'bold');
end

% Draw CD bar
proposed_rank_cd = overall_friedman_ranks(proposed_idx);
cd_x_start = proposed_rank_cd;
cd_x_end   = proposed_rank_cd + CD;
cd_y       = y_bottom + 0.5;

plot(ax_cd, [cd_x_start cd_x_end], [cd_y cd_y], ...
    'k-', 'LineWidth', 3);
plot(ax_cd, [cd_x_start cd_x_start], ...
    [cd_y - 0.15 cd_y + 0.15], 'k-', 'LineWidth', 2);
plot(ax_cd, [cd_x_end cd_x_end], ...
    [cd_y - 0.15 cd_y + 0.15], 'k-', 'LineWidth', 2);
text(ax_cd, (cd_x_start + cd_x_end) / 2, cd_y + 0.3, ...
    sprintf('CD = %.3f', CD), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');

% Mark significantly different from proposed
for i = 1:k_algs
    j = sort_fr_idx(i);
    if j ~= proposed_idx
        rank_diff = abs(overall_friedman_ranks(j) - ...
                        overall_friedman_ranks(proposed_idx));
        if rank_diff > CD
            plot(ax_cd, sorted_fr(i), y_top, 'r*', 'MarkerSize', 12);
        end
    end
end

title(ax_cd, sprintf('Critical Difference Diagram (Nemenyi, \\alpha=%.2f)', alpha), ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax_cd, 'Average Rank (lower is better)', 'FontSize', 11);

axis(ax_cd, 'off');
set(ax_cd, 'XLim', [x_left - 2, x_right + 2], ...
           'YLim', [y_bottom - 0.5, y_top + 1.5]);

hold(ax_cd, 'off');

if save_figures
    cd_path_base = fullfile(results_folder, 'CD_Diagrams', 'CD_Diagram');
    save_figure_journal_safe(fig_cd, cd_path_base, output_format, save_png);
    fprintf('CD Diagram saved (PNG + %s)\n', upper(output_format));
end

%% Performance Profiles (Dolan & More, 2002)
% =========================================================================
fprintf('\n========================================\n');
fprintf('PERFORMANCE PROFILES (Dolan & More)\n');
fprintf('========================================\n\n');

% Performance ratio matrix
perf_matrix = zeros(num_algorithms, num_problems);

for prob_idx = 1:num_problems
    prob_medians_pp = median_values(:, prob_idx);
    min_perf        = min(prob_medians_pp(isfinite(prob_medians_pp)));

    for algo_idx = 1:num_algorithms
        if isfinite(prob_medians_pp(algo_idx)) && ...
           isfinite(min_perf) && min_perf > 0
            perf_matrix(algo_idx, prob_idx) = ...
                prob_medians_pp(algo_idx) / min_perf;
        elseif prob_medians_pp(algo_idx) == 0 && min_perf == 0
            perf_matrix(algo_idx, prob_idx) = 1;
        else
            perf_matrix(algo_idx, prob_idx) = inf;
        end
    end
end

% Performance profile function
tau_max  = 10;
tau_vals = linspace(1, tau_max, 500);
rho      = zeros(num_algorithms, length(tau_vals));

for algo_idx = 1:num_algorithms
    for t_idx = 1:length(tau_vals)
        tau = tau_vals(t_idx);
        rho(algo_idx, t_idx) = sum(perf_matrix(algo_idx, :) <= tau) / num_problems;
    end
end

% --- Plot Performance Profiles (JOURNAL FORMAT: PNG + PDF)---
fig_pp = figure('Name', 'Performance Profiles', ...
    'Position', [100 100 900 600], ...
    'Color', 'white');

ax_pp = axes('Parent', fig_pp);
hold(ax_pp, 'on');

line_styles = {'-', '--', '-.', ':', '-', '--', '-.', ':'};
colors_pp   = lines(num_algorithms);

pp_handles = zeros(num_algorithms, 1);

for algo_idx = 1:num_algorithms
    pp_handles(algo_idx) = stairs(ax_pp, ...
        tau_vals, rho(algo_idx, :), ...
        'Color',     colors_pp(algo_idx, :), ...
        'LineStyle', line_styles{mod(algo_idx-1, length(line_styles)) + 1}, ...
        'LineWidth', 2.0, ...
        'DisplayName', algorithms{algo_idx, 1});

    if algo_idx == proposed_idx
        set(pp_handles(algo_idx), 'LineWidth', 3.5);
    end
end

% Reference lines
plot(ax_pp, [1 1], [0 1], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(ax_pp, [1 tau_max], [1 1], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

hold(ax_pp, 'off');

xlim(ax_pp, [1, tau_max]);
ylim(ax_pp, [0, 1.05]);
grid(ax_pp, 'on');
box(ax_pp, 'on');

xlabel(ax_pp, '\tau (Performance Ratio)', 'FontSize', 12);
ylabel(ax_pp, '\rho_s(\tau) - Fraction of Problems Solved', 'FontSize', 12);
title(ax_pp, 'Performance Profiles (Dolan & More, 2002)', ...
    'FontSize', 14, 'FontWeight', 'bold');

lgd_pp = legend(ax_pp, 'show');
lgd_pp.Location = 'southeast';
lgd_pp.FontSize = 10;
lgd_pp.Box      = 'on';

set(ax_pp, 'FontSize', 11);

if save_figures
    pp_path_base = fullfile(results_folder, 'Performance_Profiles', 'Performance_Profiles');
    save_figure_journal_safe(fig_pp, pp_path_base, output_format, save_png);
    fprintf('Performance Profiles saved (PNG + %s)\n', upper(output_format));
end

%% Convergence Curves (JOURNAL FORMAT: PNG + PDF)
% =========================================================================
fprintf('\n========================================\n');
fprintf('CONVERGENCE CURVES\n');
fprintf('========================================\n\n');

colors_conv = lines(num_algorithms);

for prob_idx = 1:num_problems
    try
        fig_conv = figure('Name', ...
            sprintf('Convergence - Problem %d', prob_idx), ...
            'Position', [50 50 900 550], ...
            'Color', 'white');

        ax_conv = axes('Parent', fig_conv);
        hold(ax_conv, 'on');

        conv_handles = zeros(num_algorithms, 1);

        for algo_idx = 1:num_algorithms
            all_curves   = {};
            max_len      = 0;

            % Collect valid curves
            for run_idx = 1:num_runs
                cv = convergence_curves{algo_idx, prob_idx, run_idx};
                if ~isempty(cv) && isvector(cv) && numel(cv) > 1
                    all_curves{end+1} = cv(:)';
                    max_len = max(max_len, length(cv));
                end
            end

            if isempty(all_curves)
                continue;
            end

            % Pad curves to equal length
            padded = NaN(length(all_curves), max_len);
            for c = 1:length(all_curves)
                padded(c, 1:length(all_curves{c})) = all_curves{c};
            end

            % Mean convergence
            mean_conv = zeros(1, max_len);
            for col = 1:max_len
                col_data = padded(:, col);
                col_data = col_data(isfinite(col_data));
                if ~isempty(col_data)
                    mean_conv(col) = mean(col_data);
                else
                    mean_conv(col) = NaN;
                end
            end

            std_conv = zeros(1, max_len);
            for col = 1:max_len
                col_data = padded(:, col);
                col_data = col_data(isfinite(col_data));
                if length(col_data) > 1
                    std_conv(col) = std(col_data);
                else
                    std_conv(col) = 0;
                end
            end

            % Remove leading NaNs
            valid_idx = find(isfinite(mean_conv), 1, 'first');
            if isempty(valid_idx)
                continue;
            end
            
            mean_conv = mean_conv(valid_idx:end);
            std_conv  = std_conv(valid_idx:end);
            x_iter    = (valid_idx:max_len) - valid_idx + 1;

            ls     = line_styles{mod(algo_idx-1, length(line_styles)) + 1};
            clr    = colors_conv(algo_idx, :);
            lw     = 1.8;
            if algo_idx == proposed_idx
                lw = 3.0;
            end

            % Shaded std region
            if length(x_iter) > 2
                x_fill  = [x_iter, fliplr(x_iter)];
                y_upper = mean_conv + std_conv;
                y_lower = mean_conv - std_conv;
                y_fill  = [y_upper, fliplr(y_lower)];

                valid_fill = isfinite(y_fill);
                if sum(valid_fill) > 2
                    fill(ax_conv, ...
                        x_fill(valid_fill), y_fill(valid_fill), clr, ...
                        'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
                        'HandleVisibility', 'off');
                end
            end

            conv_handles(algo_idx) = plot(ax_conv, x_iter, mean_conv, ...
                'Color', clr, 'LineStyle', ls, 'LineWidth', lw, ...
                'DisplayName', algorithms{algo_idx, 1});
        end

        hold(ax_conv, 'off');

        set(ax_conv, 'YScale', 'log', 'FontSize', 11);
        grid(ax_conv, 'on');
        box(ax_conv, 'on');

        xlabel(ax_conv, 'Iteration', 'FontSize', 12);
        ylabel(ax_conv, 'Best Fitness (log scale)', 'FontSize', 12);
        title(ax_conv, ...
            sprintf('Convergence Curves - %s', ...
                problem_names{current_problems(prob_idx)}), ...
            'FontSize', 13, 'FontWeight', 'bold');

        valid_h = conv_handles(conv_handles ~= 0);
        if ~isempty(valid_h)
            lgd_conv = legend(ax_conv, valid_h, ...
                algorithms(conv_handles ~= 0, 1));
            lgd_conv.Location = 'northeast';
            lgd_conv.FontSize = 9;
            lgd_conv.Box      = 'on';
        end
        
        drawnow; % Force render

        if save_figures
            conv_path_base = fullfile(results_folder, 'Convergence_Curves', ...
                sprintf('Convergence_P%d', prob_idx));
            save_figure_journal_safe(fig_conv, conv_path_base, output_format, save_png);
            fprintf('Convergence curve P%d saved (PNG + %s)\n', prob_idx, upper(output_format));
        end
        
    catch ME_conv
        fprintf('Warning: Could not save convergence curve for P%d: %s\n', prob_idx, ME_conv.message);
    end
end

%% Box Plots (JOURNAL FORMAT: PNG + PDF)
% =========================================================================
fprintf('\n========================================\n');
fprintf('BOX PLOTS\n');
fprintf('========================================\n\n');

colors_box = lines(num_algorithms);

for prob_idx = 1:num_problems
    try
        fig_box = figure('Name', ...
            sprintf('Boxplot - Problem %d', prob_idx), ...
            'Position', [50 50 1000 500], ...
            'Color', 'white');

        ax_box = axes('Parent', fig_box);

        % Collect data
        box_data  = zeros(num_runs, num_algorithms);
        box_mask  = false(num_runs, num_algorithms);

        for algo_idx = 1:num_algorithms
            vals = squeeze(best_values(algo_idx, prob_idx, :));
            finite_mask = isfinite(vals);
            box_data(finite_mask, algo_idx) = vals(finite_mask);
            box_mask(:, algo_idx)           = finite_mask;
        end

        % Build cell array for boxplot
        all_vals  = [];
        all_grps  = [];

        for algo_idx = 1:num_algorithms
            v = box_data(box_mask(:, algo_idx), algo_idx);
            all_vals = [all_vals; v];
            all_grps = [all_grps; repmat(algo_idx, length(v), 1)];
        end

        if ~isempty(all_vals)
            boxplot(ax_box, all_vals, all_grps, ...
                'Labels',   algorithms(:, 1)', ...
                'Colors',   colors_box, ...
                'Symbol',   'r+', ...
                'Widths',   0.6, ...
                'Notch',    'off');

            % Color boxes
            box_patches = findobj(ax_box, 'Tag', 'Box');
            for bp = 1:length(box_patches)
                c_idx = num_algorithms - bp + 1;
                if c_idx >= 1 && c_idx <= num_algorithms
                    try
                        patch(get(box_patches(bp), 'XData'), ...
                              get(box_patches(bp), 'YData'), ...
                              colors_box(c_idx, :), ...
                              'FaceAlpha', 0.4, ...
                              'Parent', ax_box);
                    catch
                        % Skip if patch fails
                    end
                end
            end

            set(ax_box, 'YScale', 'log', 'FontSize', 10);
            grid(ax_box, 'on');
            box(ax_box, 'on');
            xtickangle(ax_box, 30);

            ylabel(ax_box, 'Best Fitness (log scale)', 'FontSize', 12);
            title(ax_box, ...
                sprintf('Box Plot - %s', ...
                    problem_names{current_problems(prob_idx)}), ...
                'FontSize', 13, 'FontWeight', 'bold');

            % Highlight proposed
            ax_box.XTickLabel{proposed_idx} = ...
                sprintf('\\bf{%s}', algorithms{proposed_idx,1});
        end
        
        drawnow; % Force render

        if save_figures
            box_path_base = fullfile(results_folder, 'Boxplots', ...
                sprintf('Boxplot_P%d', prob_idx));
            save_figure_journal_safe(fig_box, box_path_base, output_format, save_png);
            fprintf('Boxplot P%d saved (PNG + %s)\n', prob_idx, upper(output_format));
        end
        
    catch ME_box
        fprintf('Warning: Could not save boxplot for P%d: %s\n', prob_idx, ME_box.message);
    end
end

%% Runtime Analysis
% =========================================================================
fprintf('\n========================================\n');
fprintf('RUNTIME ANALYSIS\n');
fprintf('========================================\n\n');

% Mean runtime per algorithm per problem
mean_runtime = zeros(num_algorithms, num_problems);
std_runtime  = zeros(num_algorithms, num_problems);

for algo_idx = 1:num_algorithms
    for prob_idx = 1:num_problems
        rt_vals = squeeze(execution_times(algo_idx, prob_idx, :));
        rt_vals = rt_vals(isfinite(rt_vals));
        if ~isempty(rt_vals)
            mean_runtime(algo_idx, prob_idx) = mean(rt_vals);
            std_runtime(algo_idx, prob_idx)  = std(rt_vals);
        end
    end
end

overall_mean_rt = zeros(num_algorithms, 1);
for algo_idx = 1:num_algorithms
    all_rt = squeeze(execution_times(algo_idx, :, :));
    all_rt = all_rt(:);
    all_rt = all_rt(isfinite(all_rt));
    if ~isempty(all_rt)
        overall_mean_rt(algo_idx) = mean(all_rt);
    end
end

fprintf('Mean execution time per algorithm (seconds):\n');
fprintf('%-12s', 'Algorithm');
for prob_idx = 1:num_problems
    fprintf('  P%-8d', prob_idx);
end
fprintf('  Overall\n');
fprintf('%s\n', repmat('-', 1, 12 + (num_problems+1)*11));

for algo_idx = 1:num_algorithms
    fprintf('%-12s', algorithms{algo_idx,1});
    for prob_idx = 1:num_problems
        fprintf('  %8.3f  ', mean_runtime(algo_idx, prob_idx));
    end
    fprintf('  %8.3f\n', overall_mean_rt(algo_idx));
end

% --- Runtime Bar Chart (JOURNAL FORMAT: PNG + PDF)---
try
    fig_rt = figure('Name', 'Runtime Analysis', ...
        'Position', [100 100 900 500], ...
        'Color', 'white');

    ax_rt = axes('Parent', fig_rt);

    b_rt = bar(ax_rt, mean_runtime', 'grouped');

    for algo_idx = 1:num_algorithms
        b_rt(algo_idx).FaceColor = colors_conv(algo_idx, :);
        b_rt(algo_idx).FaceAlpha = 0.8;
    end

    set(ax_rt, 'XTickLabel', ...
        arrayfun(@(x) sprintf('P%d', x), 1:num_problems, ...
        'UniformOutput', false), 'FontSize', 11);
    grid(ax_rt, 'on');
    box(ax_rt, 'on');

    xlabel(ax_rt, 'Problem', 'FontSize', 12);
    ylabel(ax_rt, 'Mean Execution Time (s)', 'FontSize', 12);
    title(ax_rt, 'Runtime Comparison (Hard Problems - Fast)', ...
        'FontSize', 14, 'FontWeight', 'bold');

    lgd_rt = legend(ax_rt, algorithms(:,1)');
    lgd_rt.Location  = 'northeast';
    lgd_rt.FontSize  = 9;
    lgd_rt.Box       = 'on';
    
    drawnow; % Force render

    if save_figures
        rt_path_base = fullfile(results_folder, 'Summary_Statistics', 'Runtime_Comparison');
        save_figure_journal_safe(fig_rt, rt_path_base, output_format, save_png);
        fprintf('Runtime chart saved (PNG + %s)\n', upper(output_format));
    end
    
catch ME_rt
    fprintf('Warning: Could not save runtime chart: %s\n', ME_rt.message);
end

%% AUC Analysis (Convergence Speed)
% =========================================================================
fprintf('\n========================================\n');
fprintf('AUC ANALYSIS (Convergence Speed)\n');
fprintf('========================================\n\n');

mean_auc = zeros(num_algorithms, num_problems);
std_auc  = zeros(num_algorithms, num_problems);

for algo_idx = 1:num_algorithms
    for prob_idx = 1:num_problems
        auc_v = squeeze(auc_values(algo_idx, prob_idx, :));
        auc_v = auc_v(isfinite(auc_v));
        if ~isempty(auc_v)
            mean_auc(algo_idx, prob_idx) = mean(auc_v);
            std_auc(algo_idx, prob_idx)  = std(auc_v);
        end
    end
end

fprintf('Mean AUC (lower = faster convergence):\n');
fprintf('%-12s', 'Algorithm');
for prob_idx = 1:num_problems
    fprintf('  P%-8d', prob_idx);
end
fprintf('\n%s\n', repmat('-', 1, 12 + num_problems*11));

for algo_idx = 1:num_algorithms
    fprintf('%-12s', algorithms{algo_idx,1});
    for prob_idx = 1:num_problems
        fprintf('  %8.4f  ', mean_auc(algo_idx, prob_idx));
    end
    fprintf('\n');
end

% --- AUC Heatmap (JOURNAL FORMAT: PNG + PDF)---
try
    fig_auc = figure('Name', 'AUC Heatmap', ...
        'Position', [100 100 700 500], ...
        'Color', 'white');

    ax_auc = axes('Parent', fig_auc);

    imagesc(ax_auc, mean_auc);
    colormap(ax_auc, 'cool');
    colorbar(ax_auc);

    for i = 1:num_algorithms
        for j = 1:num_problems
            if mean_auc(i, j) > 0
                text(ax_auc, j, i, sprintf('%.3f', mean_auc(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'middle', ...
                    'FontSize', 9, 'FontWeight', 'bold', ...
                    'Color', 'k');
            end
        end
    end

    set(ax_auc, ...
        'XTick', 1:num_problems, ...
        'XTickLabel', arrayfun(@(x) sprintf('P%d', x), ...
            1:num_problems, 'UniformOutput', false), ...
        'YTick', 1:num_algorithms, ...
        'YTickLabel', algorithms(:, 1), ...
        'FontSize', 10);

    xlabel(ax_auc, 'Problem', 'FontSize', 12);
    ylabel(ax_auc, 'Algorithm', 'FontSize', 12);
    title(ax_auc, 'AUC Heatmap (lower = faster) - Hard Problems (Fast)', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow; % Force render

    if save_figures
        auc_path_base = fullfile(results_folder, 'Summary_Statistics', 'AUC_Heatmap');
        save_figure_journal_safe(fig_auc, auc_path_base, output_format, save_png);
        fprintf('AUC heatmap saved (PNG + %s)\n', upper(output_format));
    end
    
catch ME_auc
    fprintf('Warning: Could not save AUC heatmap: %s\n', ME_auc.message);
end

%% Rank Heatmap (JOURNAL FORMAT: PNG + PDF)
% =========================================================================
fprintf('\n========================================\n');
fprintf('RANK HEATMAP\n');
fprintf('========================================\n\n');

try
    fig_rank = figure('Name', 'Rank Heatmap', ...
        'Position', [100 100 700 500], ...
        'Color', 'white');

    ax_rank = axes('Parent', fig_rank);

    imagesc(ax_rank, friedman_ranks);
    colormap(ax_rank, 'summer');
    colorbar(ax_rank);

    for i = 1:num_algorithms
        for j = 1:num_problems
            if friedman_ranks(i, j) > 0
                text(ax_rank, j, i, sprintf('%.2f', friedman_ranks(i, j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment',   'middle', ...
                    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end

    set(ax_rank, ...
        'XTick', 1:num_problems, ...
        'XTickLabel', arrayfun(@(x) sprintf('P%d', x), ...
            1:num_problems, 'UniformOutput', false), ...
        'YTick', 1:num_algorithms, ...
        'YTickLabel', algorithms(:, 1), ...
        'FontSize', 10);

    xlabel(ax_rank, 'Problem (Hard Versions - Fast)', 'FontSize', 12);
    ylabel(ax_rank, 'Algorithm', 'FontSize', 12);
    title(ax_rank, 'Friedman Rank Heatmap (lower is better) - Fast', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    drawnow; % Force render

    if save_figures
        rank_path_base = fullfile(results_folder, 'Summary_Statistics', 'Rank_Heatmap');
        save_figure_journal_safe(fig_rank, rank_path_base, output_format, save_png);
        fprintf('Rank heatmap saved (PNG + %s)\n', upper(output_format));
    end
    
catch ME_rank
    fprintf('Warning: Could not save rank heatmap: %s\n', ME_rank.message);
end

%% Excel Results - Compatible with all MATLAB versions
% =========================================================================
fprintf('\n========================================\n');
fprintf('SAVING EXCEL RESULTS\n');
fprintf('========================================\n\n');

excel_filename = fullfile(results_folder, 'Excel_Results', ...
    'Electronics_Results_Hard_Fast.xlsx');

try
    % --- Sheet 1: Main Results Table ---
    main_headers = {'Algorithm', 'Problem', ...
        'Min', 'Max', 'Mean', 'Median', 'Std', ...
        'Rank', 'Friedman_Rank'};

    main_data = {};
    row_idx = 1;

    for prob_idx = 1:num_problems
        for algo_idx = 1:num_algorithms
            main_data{row_idx, 1} = algorithms{algo_idx, 1};
            main_data{row_idx, 2} = problem_names{current_problems(prob_idx)};
            main_data{row_idx, 3} = min_values(algo_idx, prob_idx);
            main_data{row_idx, 4} = max_values(algo_idx, prob_idx);
            main_data{row_idx, 5} = mean_values(algo_idx, prob_idx);
            main_data{row_idx, 6} = median_values(algo_idx, prob_idx);
            main_data{row_idx, 7} = std_values(algo_idx, prob_idx);
            main_data{row_idx, 8} = ranks(algo_idx, prob_idx);
            main_data{row_idx, 9} = friedman_ranks(algo_idx, prob_idx);
            row_idx = row_idx + 1;
        end
    end

    writecell_safe([main_headers; main_data], excel_filename, 'Main_Results');
    fprintf('Sheet "Main_Results" saved.\n');

    % --- Sheet 2: Summary Statistics ---
    summary_headers = {'Algorithm', ...
        'Avg_Rank', 'Friedman_Avg_Rank', ...
        'Overall_Mean_Time_s', 'Overall_Mean_AUC'};

    summary_data = {};

    for algo_idx = 1:num_algorithms
        rt_all = execution_times(algo_idx, :, :);
        rt_all = rt_all(:);
        rt_all = rt_all(isfinite(rt_all));
        if ~isempty(rt_all)
            mean_rt_overall = mean(rt_all);
        else
            mean_rt_overall = NaN;
        end

        auc_all = auc_values(algo_idx, :, :);
        auc_all = auc_all(:);
        auc_all = auc_all(isfinite(auc_all));
        if ~isempty(auc_all)
            mean_auc_overall = mean(auc_all);
        else
            mean_auc_overall = NaN;
        end

        summary_data{algo_idx, 1} = algorithms{algo_idx, 1};
        summary_data{algo_idx, 2} = avg_ranks(algo_idx);
        summary_data{algo_idx, 3} = overall_friedman_ranks(algo_idx);
        summary_data{algo_idx, 4} = mean_rt_overall;
        summary_data{algo_idx, 5} = mean_auc_overall;
    end

    writecell_safe([summary_headers; summary_data], excel_filename, 'Summary');
    fprintf('Sheet "Summary" saved.\n');

    % --- Sheet 3: Advanced Wilcoxon Results vs Proposed ---
    wilcoxon_headers = {'Algorithm'};
    for prob_idx = 1:num_problems
        wilcoxon_headers{end+1} = sprintf('P%d_Adjusted_p', prob_idx);
        wilcoxon_headers{end+1} = sprintf('P%d_Effect_Size', prob_idx);
        wilcoxon_headers{end+1} = sprintf('P%d_Result', prob_idx);
    end

    wilcoxon_data = {};
    row_w = 1;

    for a = 1:length(other_indices)
        comp_idx = other_indices(a);
        wilcoxon_data{row_w, 1} = algorithms{comp_idx, 1};
        col_w = 2;

        for prob_idx = 1:num_problems
            p_adj = P_holm(a, prob_idx);
            d = EffectSize(a, prob_idx);
            
            if isnan(p_adj)
                result_str = 'N/A';
            elseif p_adj >= alpha
                result_str = '=';
            else
                kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
                comp_mean = mean(squeeze(best_values(comp_idx, prob_idx, :)), 'omitnan');
                if kdo_mean < comp_mean, result_str = '+'; else, result_str = '-'; end
            end
            
            wilcoxon_data{row_w, col_w}   = p_adj;
            wilcoxon_data{row_w, col_w+1} = d;
            wilcoxon_data{row_w, col_w+2} = result_str;
            col_w = col_w + 3;
        end
        row_w = row_w + 1;
    end

    writecell_safe([wilcoxon_headers; wilcoxon_data], ...
        excel_filename, 'Wilcoxon_Advanced');
    fprintf('Sheet "Wilcoxon_Advanced" saved.\n');

    % --- Sheet 4: Friedman Test Results ---
    friedman_headers = {'Problem', 'p_value', 'Significant'};
    friedman_export = {};

    for prob_idx = 1:num_problems
        p_f = friedman_stats(prob_idx).p_value;
        friedman_export{prob_idx, 1} = ...
            problem_names{current_problems(prob_idx)};
        friedman_export{prob_idx, 2} = p_f;
        if p_f < alpha
            friedman_export{prob_idx, 3} = 'Yes';
        else
            friedman_export{prob_idx, 3} = 'No';
        end
    end

    writecell_safe([friedman_headers; friedman_export], ...
        excel_filename, 'Friedman_Tests');
    fprintf('Sheet "Friedman_Tests" saved.\n');

    % --- Sheet 5: Raw Best Values ---
    raw_headers = {'Algorithm', 'Problem'};
    for r = 1:num_runs
        raw_headers{end+1} = sprintf('Run_%d', r);
    end

    raw_data = {};
    row_r = 1;

    for prob_idx = 1:num_problems
        for algo_idx = 1:num_algorithms
            raw_data{row_r, 1} = algorithms{algo_idx, 1};
            raw_data{row_r, 2} = ...
                problem_names{current_problems(prob_idx)};

            for run_idx = 1:num_runs
                raw_data{row_r, 2 + run_idx} = ...
                    best_values(algo_idx, prob_idx, run_idx);
            end
            row_r = row_r + 1;
        end
    end

    writecell_safe([raw_headers; raw_data], excel_filename, 'Raw_Data');
    fprintf('Sheet "Raw_Data" saved.\n');

    % --- Sheet 6: Runtime Data ---
    rt_headers = {'Algorithm', 'Problem', 'Mean_Time_s', 'Std_Time_s'};
    rt_data = {};
    row_rt = 1;

    for prob_idx = 1:num_problems
        for algo_idx = 1:num_algorithms
            rt_data{row_rt, 1} = algorithms{algo_idx, 1};
            rt_data{row_rt, 2} = ...
                problem_names{current_problems(prob_idx)};
            rt_data{row_rt, 3} = mean_runtime(algo_idx, prob_idx);
            rt_data{row_rt, 4} = std_runtime(algo_idx, prob_idx);
            row_rt = row_rt + 1;
        end
    end

    writecell_safe([rt_headers; rt_data], excel_filename, 'Runtime');
    fprintf('Sheet "Runtime" saved.\n');

    % --- Sheet 7: Physical Parameters ---
    phys_headers = {'Algorithm', 'Problem'};
    % Add parameter names for each problem
    % This is simplified - in practice you'd create separate sheets per problem
    
    writecell_safe({'Physical parameters saved in Physical_Parameters.mat'}, ...
        excel_filename, 'Physical_Params');
    fprintf('Sheet "Physical_Params" saved (reference).\n');

    fprintf('\nAll Excel sheets saved to:\n  %s\n', excel_filename);

catch ME_excel
    fprintf('Warning: Excel export failed: %s\n', ME_excel.message);
    fprintf('Saving results as CSV instead...\n');

    csv_file = fullfile(results_folder, 'Excel_Results', ...
        'Summary_Results.csv');
    fid = fopen(csv_file, 'w');
    if fid ~= -1
        fprintf(fid, 'Algorithm,Avg_Rank,Friedman_Avg_Rank\n');
        for algo_idx = 1:num_algorithms
            fprintf(fid, '%s,%.4f,%.4f\n', ...
                algorithms{algo_idx,1}, ...
                avg_ranks(algo_idx), ...
                overall_friedman_ranks(algo_idx));
        end
        fclose(fid);
        fprintf('CSV fallback saved: %s\n', csv_file);
    end
end

%% LaTeX Tables Generation (Complete)
% =========================================================================
fprintf('\n========================================\n');
fprintf('GENERATING LATEX TABLES\n');
fprintf('========================================\n\n');

latex_dir = fullfile(results_folder, 'Latex_Tables');
if ~exist(latex_dir, 'dir')
    mkdir(latex_dir);
end

% ==================== TABLE 1: Main Results (Mean ± Std) ====================
fprintf('  Generating Table 1: Main Results...\n');

latex_file_main = fullfile(latex_dir, 'Table1_Main_Results.tex');
fid_main = fopen(latex_file_main, 'w');

if fid_main ~= -1
    fprintf(fid_main, '%%%% Table 1: Main Optimization Results (Mean ± Std)\n');
    fprintf(fid_main, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_main, '\\begin{table}[htbp]\n');
    fprintf(fid_main, '\\centering\n');
    fprintf(fid_main, '\\caption{Main Optimization Results: Mean $\\pm$ Std over %d independent runs (HARD VERSION - FAST)}\n', num_runs);
    fprintf(fid_main, '\\label{tab:main_results}\n');
    fprintf(fid_main, '\\small\n');
    fprintf(fid_main, '\\setlength{\\tabcolsep}{4pt}\n');
    fprintf(fid_main, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_main, 'c');
    end
    fprintf(fid_main, '}\n');
    fprintf(fid_main, '\\toprule\n');
    
    fprintf(fid_main, '\\textbf{Algorithm}');
    for prob_idx = 1:num_problems
        fprintf(fid_main, ' & \\textbf{P%d}', prob_idx);
    end
    fprintf(fid_main, ' \\\\\n');
    fprintf(fid_main, '\\midrule\n');
    
    % Find best per problem
    best_per_prob = zeros(1, num_problems);
    for prob_idx = 1:num_problems
        [~, best_per_prob(prob_idx)] = min(mean_values(:, prob_idx));
    end
    
    for algo_idx = 1:num_algorithms
        if algo_idx == proposed_idx
            fprintf(fid_main, '\\textbf{%s}', algorithms{algo_idx, 1});
        else
            fprintf(fid_main, '%s', algorithms{algo_idx, 1});
        end
        
        for prob_idx = 1:num_problems
            if mean_values(algo_idx, prob_idx) > 0
                if algo_idx == best_per_prob(prob_idx)
                    fprintf(fid_main, ' & \\textbf{$%.2e \\pm %.2e$}', ...
                        mean_values(algo_idx, prob_idx), std_values(algo_idx, prob_idx));
                else
                    fprintf(fid_main, ' & $%.2e \\pm %.2e$', ...
                        mean_values(algo_idx, prob_idx), std_values(algo_idx, prob_idx));
                end
            else
                fprintf(fid_main, ' & --');
            end
        end
        fprintf(fid_main, ' \\\\\n');
    end
    
    fprintf(fid_main, '\\bottomrule\n');
    fprintf(fid_main, '\\end{tabular}\n');
    fprintf(fid_main, '\\end{table}\n');
    fclose(fid_main);
    fprintf('    Saved: %s\n', latex_file_main);
end

% ==================== TABLE 2: Best, Median, Worst ====================
fprintf('  Generating Table 2: Best/Median/Worst Results...\n');

latex_file_bmw = fullfile(latex_dir, 'Table2_Best_Median_Worst.tex');
fid_bmw = fopen(latex_file_bmw, 'w');

if fid_bmw ~= -1
    fprintf(fid_bmw, '%%%% Table 2: Best, Median, and Worst Results\n');
    fprintf(fid_bmw, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_bmw, '\\begin{sidewaystable}[htbp]\n');
    fprintf(fid_bmw, '\\centering\n');
    fprintf(fid_bmw, '\\caption{Best, Median, and Worst Fitness Values over %d Runs (HARD VERSION - FAST)}\n', num_runs);
    fprintf(fid_bmw, '\\label{tab:best_median_worst}\n');
    fprintf(fid_bmw, '\\tiny\n');
    fprintf(fid_bmw, '\\setlength{\\tabcolsep}{2pt}\n');
    fprintf(fid_bmw, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_bmw, 'ccc');
    end
    fprintf(fid_bmw, '}\n');
    fprintf(fid_bmw, '\\toprule\n');
    
    fprintf(fid_bmw, '\\multirow{2}{*}{\\textbf{Algorithm}}');
    for prob_idx = 1:num_problems
        fprintf(fid_bmw, ' & \\multicolumn{3}{c}{\\textbf{P%d}}', prob_idx);
    end
    fprintf(fid_bmw, ' \\\\\n');
    
    fprintf(fid_bmw, '\\cmidrule(lr){2-4}');
    for prob_idx = 2:num_problems
        fprintf(fid_bmw, '\\cmidrule(lr){%d-%d}', (prob_idx-1)*3+2, prob_idx*3+1);
    end
    fprintf(fid_bmw, '\n');
    
    fprintf(fid_bmw, ' & \\textbf{Best} & \\textbf{Median} & \\textbf{Worst}');
    for prob_idx = 2:num_problems
        fprintf(fid_bmw, ' & \\textbf{Best} & \\textbf{Median} & \\textbf{Worst}');
    end
    fprintf(fid_bmw, ' \\\\\n');
    fprintf(fid_bmw, '\\midrule\n');
    
    for algo_idx = 1:num_algorithms
        if algo_idx == proposed_idx
            fprintf(fid_bmw, '\\textbf{%s}', algorithms{algo_idx, 1});
        else
            fprintf(fid_bmw, '%s', algorithms{algo_idx, 1});
        end
        
        for prob_idx = 1:num_problems
            fprintf(fid_bmw, ' & %.2e & %.2e & %.2e', ...
                min_values(algo_idx, prob_idx), ...
                median_values(algo_idx, prob_idx), ...
                max_values(algo_idx, prob_idx));
        end
        fprintf(fid_bmw, ' \\\\\n');
    end
    
    fprintf(fid_bmw, '\\bottomrule\n');
    fprintf(fid_bmw, '\\end{tabular}\n');
    fprintf(fid_bmw, '\\end{sidewaystable}\n');
    fclose(fid_bmw);
    fprintf('    Saved: %s\n', latex_file_bmw);
end

% ==================== TABLE 3: Friedman Ranks ====================
fprintf('  Generating Table 3: Friedman Ranks...\n');

latex_file_friedman = fullfile(latex_dir, 'Table3_Friedman_Ranks.tex');
fid_friedman = fopen(latex_file_friedman, 'w');

if fid_friedman ~= -1
    fprintf(fid_friedman, '%%%% Table 3: Friedman Average Ranks\n');
    fprintf(fid_friedman, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_friedman, '\\begin{table}[htbp]\n');
    fprintf(fid_friedman, '\\centering\n');
    fprintf(fid_friedman, '\\caption{Friedman Average Ranks (lower is better) - HARD VERSION (FAST)}\n');
    fprintf(fid_friedman, '\\label{tab:friedman_ranks}\n');
    fprintf(fid_friedman, '\\small\n');
    fprintf(fid_friedman, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_friedman, 'c');
    end
    fprintf(fid_friedman, 'c}\n');
    fprintf(fid_friedman, '\\toprule\n');
    
    fprintf(fid_friedman, '\\textbf{Algorithm}');
    for prob_idx = 1:num_problems
        fprintf(fid_friedman, ' & \\textbf{P%d}', prob_idx);
    end
    fprintf(fid_friedman, ' & \\textbf{Overall} \\\\\n');
    fprintf(fid_friedman, '\\midrule\n');
    
    [~, rank_sort_idx] = sort(overall_friedman_ranks);
    
    for si = 1:num_algorithms
        algo_idx = rank_sort_idx(si);
        
        if algo_idx == proposed_idx
            fprintf(fid_friedman, '\\textbf{%s}', algorithms{algo_idx, 1});
        else
            fprintf(fid_friedman, '%s', algorithms{algo_idx, 1});
        end
        
        for prob_idx = 1:num_problems
            fprintf(fid_friedman, ' & %.2f', friedman_ranks(algo_idx, prob_idx));
        end
        if algo_idx == proposed_idx
            fprintf(fid_friedman, ' & \\textbf{%.2f}', overall_friedman_ranks(algo_idx));
        else
            fprintf(fid_friedman, ' & %.2f', overall_friedman_ranks(algo_idx));
        end
        fprintf(fid_friedman, ' \\\\\n');
    end
    
    fprintf(fid_friedman, '\\bottomrule\n');
    fprintf(fid_friedman, '\\end{tabular}\n');
    fprintf(fid_friedman, '\\end{table}\n');
    fclose(fid_friedman);
    fprintf('    Saved: %s\n', latex_file_friedman);
end

% ==================== TABLE 4: Advanced Wilcoxon Test Results ====================
fprintf('  Generating Table 4: Wilcoxon Test Results...\n');

latex_file_wilcoxon = fullfile(latex_dir, 'Table4_Wilcoxon_Advanced.tex');
fid_wilcoxon = fopen(latex_file_wilcoxon, 'w');

if fid_wilcoxon ~= -1
    fprintf(fid_wilcoxon, '%%%% Table 4: Wilcoxon Signed-Rank Test Results vs Proposed Algorithm (Holm-Bonferroni Corrected)\n');
    fprintf(fid_wilcoxon, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_wilcoxon, '\\begin{table}[htbp]\n');
    fprintf(fid_wilcoxon, '\\centering\n');
    fprintf(fid_wilcoxon, '\\caption{Paired Wilcoxon Signed-Rank Test Results (Holm-Bonferroni corrected, $\\alpha=%.2f$) vs. \\textbf{%s}. +: %s better, =: no difference, -: %s worse. Effect size (Cliff''s $\\delta$) in brackets. (HARD VERSION - FAST)}\n', ...
        alpha, algorithms{proposed_idx,1}, algorithms{proposed_idx,1}, algorithms{proposed_idx,1});
    fprintf(fid_wilcoxon, '\\label{tab:wilcoxon_advanced}\n');
    fprintf(fid_wilcoxon, '\\small\n');
    fprintf(fid_wilcoxon, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_wilcoxon, 'c');
    end
    fprintf(fid_wilcoxon, 'c}\n');
    fprintf(fid_wilcoxon, '\\toprule\n');
    
    fprintf(fid_wilcoxon, '\\textbf{Algorithm}');
    for prob_idx = 1:num_problems
        fprintf(fid_wilcoxon, ' & \\textbf{P%d}', prob_idx);
    end
    fprintf(fid_wilcoxon, ' & \\textbf{W/T/L} \\\\\n');
    fprintf(fid_wilcoxon, '\\midrule\n');
    
    for a = 1:length(other_indices)
        comp_idx = other_indices(a);
        fprintf(fid_wilcoxon, '%s', algorithms{comp_idx,1});
        w = 0; t = 0; l = 0;
        
        for prob_idx = 1:num_problems
            p_adj = P_holm(a, prob_idx);
            d = EffectSize(a, prob_idx);
            
            if isnan(p_adj)
                fprintf(fid_wilcoxon, ' & --');
            elseif p_adj >= alpha
                fprintf(fid_wilcoxon, ' & $=$');
                t = t + 1;
            else
                kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
                comp_mean = mean(squeeze(best_values(comp_idx, prob_idx, :)), 'omitnan');
                d_str = sprintf('(%.2f)', d);
                if kdo_mean < comp_mean
                    fprintf(fid_wilcoxon, ' & $+$ %s', d_str);
                    w = w + 1;
                else
                    fprintf(fid_wilcoxon, ' & $-$ %s', d_str);
                    l = l + 1;
                end
            end
        end
        
        fprintf(fid_wilcoxon, ' & %d/%d/%d \\\\\n', w, t, l);
    end
    
    fprintf(fid_wilcoxon, '\\bottomrule\n');
    fprintf(fid_wilcoxon, '\\end{tabular}\n');
    fprintf(fid_wilcoxon, '\\end{table}\n');
    fclose(fid_wilcoxon);
    fprintf('    Saved: %s\n', latex_file_wilcoxon);
end

% ==================== TABLE 5: Runtime Analysis ====================
fprintf('  Generating Table 5: Runtime Analysis...\n');

latex_file_runtime = fullfile(latex_dir, 'Table5_Runtime_Analysis.tex');
fid_runtime = fopen(latex_file_runtime, 'w');

if fid_runtime ~= -1
    fprintf(fid_runtime, '%%%% Table 5: Runtime Analysis (Mean ± Std in seconds)\n');
    fprintf(fid_runtime, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_runtime, '\\begin{table}[htbp]\n');
    fprintf(fid_runtime, '\\centering\n');
    fprintf(fid_runtime, '\\caption{Execution Time Analysis: Mean $\\pm$ Std (seconds) - HARD VERSION (FAST)}\n');
    fprintf(fid_runtime, '\\label{tab:runtime}\n');
    fprintf(fid_runtime, '\\small\n');
    fprintf(fid_runtime, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_runtime, 'c');
    end
    fprintf(fid_runtime, 'c}\n');
    fprintf(fid_runtime, '\\toprule\n');
    
    fprintf(fid_runtime, '\\textbf{Algorithm}');
    for prob_idx = 1:num_problems
        fprintf(fid_runtime, ' & \\textbf{P%d}', prob_idx);
    end
    fprintf(fid_runtime, ' & \\textbf{Overall} \\\\\n');
    fprintf(fid_runtime, '\\midrule\n');
    
    for algo_idx = 1:num_algorithms
        if algo_idx == proposed_idx
            fprintf(fid_runtime, '\\textbf{%s}', algorithms{algo_idx,1});
        else
            fprintf(fid_runtime, '%s', algorithms{algo_idx,1});
        end
        
        for prob_idx = 1:num_problems
            fprintf(fid_runtime, ' & $%.2f \\pm %.2f$', ...
                mean_runtime(algo_idx, prob_idx), std_runtime(algo_idx, prob_idx));
        end
        
        if algo_idx == proposed_idx
            fprintf(fid_runtime, ' & \\textbf{%.2f} \\\\\n', overall_mean_rt(algo_idx));
        else
            fprintf(fid_runtime, ' & %.2f \\\\\n', overall_mean_rt(algo_idx));
        end
    end
    
    fprintf(fid_runtime, '\\bottomrule\n');
    fprintf(fid_runtime, '\\end{tabular}\n');
    fprintf(fid_runtime, '\\end{table}\n');
    fclose(fid_runtime);
    fprintf('    Saved: %s\n', latex_file_runtime);
end

% ==================== TABLE 6: AUC Convergence Speed ====================
fprintf('  Generating Table 6: AUC Convergence Speed...\n');

latex_file_auc = fullfile(latex_dir, 'Table6_AUC_Convergence.tex');
fid_auc = fopen(latex_file_auc, 'w');

if fid_auc ~= -1
    fprintf(fid_auc, '%%%% Table 6: AUC Convergence Speed Analysis (lower is better)\n');
    fprintf(fid_auc, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_auc, '\\begin{table}[htbp]\n');
    fprintf(fid_auc, '\\centering\n');
    fprintf(fid_auc, '\\caption{AUC Convergence Speed Analysis: Mean $\\pm$ Std (lower = faster convergence) - HARD VERSION (FAST)}\n');
    fprintf(fid_auc, '\\label{tab:auc}\n');
    fprintf(fid_auc, '\\small\n');
    fprintf(fid_auc, '\\begin{tabular}{l');
    for prob_idx = 1:num_problems
        fprintf(fid_auc, 'c');
    end
    fprintf(fid_auc, 'c}\n');
    fprintf(fid_auc, '\\toprule\n');
    
    fprintf(fid_auc, '\\textbf{Algorithm}');
    for prob_idx = 1:num_problems
        fprintf(fid_auc, ' & \\textbf{P%d}', prob_idx);
    end
    fprintf(fid_auc, ' & \\textbf{Overall} \\\\\n');
    fprintf(fid_auc, '\\midrule\n');
    
    overall_auc = mean(mean_auc, 2);
    
    for algo_idx = 1:num_algorithms
        if algo_idx == proposed_idx
            fprintf(fid_auc, '\\textbf{%s}', algorithms{algo_idx,1});
        else
            fprintf(fid_auc, '%s', algorithms{algo_idx,1});
        end
        
        for prob_idx = 1:num_problems
            fprintf(fid_auc, ' & $%.4f \\pm %.4f$', ...
                mean_auc(algo_idx, prob_idx), std_auc(algo_idx, prob_idx));
        end
        
        if algo_idx == proposed_idx
            fprintf(fid_auc, ' & \\textbf{%.4f} \\\\\n', overall_auc(algo_idx));
        else
            fprintf(fid_auc, ' & %.4f \\\\\n', overall_auc(algo_idx));
        end
    end
    
    fprintf(fid_auc, '\\bottomrule\n');
    fprintf(fid_auc, '\\end{tabular}\n');
    fprintf(fid_auc, '\\end{table}\n');
    fclose(fid_auc);
    fprintf('    Saved: %s\n', latex_file_auc);
end

% ==================== TABLE 7: Detailed per Problem Results ====================
fprintf('  Generating Table 7: Detailed per Problem Results...\n');

for prob_idx = 1:num_problems
    latex_file_detailed = fullfile(latex_dir, sprintf('Table7_P%d_Detailed_Results.tex', prob_idx));
    fid_detailed = fopen(latex_file_detailed, 'w');
    
    if fid_detailed ~= -1
        fprintf(fid_detailed, '%%%% Table 7.%d: Detailed Results for %s\n', prob_idx, problem_names{current_problems(prob_idx)});
        fprintf(fid_detailed, '%%%% Generated: %s\n', datestr(now));
        fprintf(fid_detailed, '\\begin{table}[htbp]\n');
        fprintf(fid_detailed, '\\centering\n');
        fprintf(fid_detailed, '\\caption{Detailed Results for %s over %d Runs (HARD VERSION - FAST)}\n', ...
            problem_names{current_problems(prob_idx)}, num_runs);
        fprintf(fid_detailed, '\\label{tab:P%d_detailed}\n', prob_idx);
        fprintf(fid_detailed, '\\small\n');
        fprintf(fid_detailed, '\\begin{tabular}{lccccc}\n');
        fprintf(fid_detailed, '\\toprule\n');
        fprintf(fid_detailed, '\\textbf{Algorithm} & \\textbf{Best} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Median} & \\textbf{Worst} \\\\\n');
        fprintf(fid_detailed, '\\midrule\n');
        
        % Find best mean for this problem
        prob_means = mean_values(:, prob_idx);
        [~, best_algo_prob] = min(prob_means);
        
        for algo_idx = 1:num_algorithms
            if mean_values(algo_idx, prob_idx) > 0
                if algo_idx == best_algo_prob
                    fprintf(fid_detailed, '\\textbf{%s} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} \\\\\n', ...
                        algorithms{algo_idx,1}, min_values(algo_idx, prob_idx), ...
                        mean_values(algo_idx, prob_idx), std_values(algo_idx, prob_idx), ...
                        median_values(algo_idx, prob_idx), max_values(algo_idx, prob_idx));
                else
                    fprintf(fid_detailed, '%s & %.2e & %.2e & %.2e & %.2e & %.2e \\\\\n', ...
                        algorithms{algo_idx,1}, min_values(algo_idx, prob_idx), ...
                        mean_values(algo_idx, prob_idx), std_values(algo_idx, prob_idx), ...
                        median_values(algo_idx, prob_idx), max_values(algo_idx, prob_idx));
                end
            else
                fprintf(fid_detailed, '%s & -- & -- & -- & -- & -- \\\\\n', algorithms{algo_idx,1});
            end
        end
        
        fprintf(fid_detailed, '\\bottomrule\n');
        fprintf(fid_detailed, '\\end{tabular}\n');
        fprintf(fid_detailed, '\\end{table}\n');
        fclose(fid_detailed);
    end
end
fprintf('    Saved: 4 detailed problem tables\n');

% ==================== TABLE 8: Friedman Test Statistics ====================
fprintf('  Generating Table 8: Friedman Test Statistics...\n');

latex_file_friedman_stats = fullfile(latex_dir, 'Table8_Friedman_Statistics.tex');
fid_fstats = fopen(latex_file_friedman_stats, 'w');

if fid_fstats ~= -1
    fprintf(fid_fstats, '%%%% Table 8: Friedman Test Statistics per Problem\n');
    fprintf(fid_fstats, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_fstats, '\\begin{table}[htbp]\n');
    fprintf(fid_fstats, '\\centering\n');
    fprintf(fid_fstats, '\\caption{Friedman Test Statistics per Problem (HARD VERSION - FAST)}\n');
    fprintf(fid_fstats, '\\label{tab:friedman_stats}\n');
    fprintf(fid_fstats, '\\small\n');
    fprintf(fid_fstats, '\\begin{tabular}{lccc}\n');
    fprintf(fid_fstats, '\\toprule\n');
    fprintf(fid_fstats, '\\textbf{Problem} & \\textbf{Friedman Statistic} & \\textbf{p-value} & \\textbf{Significant} \\\\\n');
    fprintf(fid_fstats, '\\midrule\n');
    
    for prob_idx = 1:num_problems
        p_val = friedman_stats(prob_idx).p_value;
        % --- FIX: Check if table exists and is a cell ---
        if isfield(friedman_stats(prob_idx), 'table') && iscell(friedman_stats(prob_idx).table) && ~isempty(friedman_stats(prob_idx).table)
            chi2_val = friedman_stats(prob_idx).table{2, 3};
        else
            chi2_val = NaN;
        end
        sig_str = 'No';
        if p_val < alpha && ~isnan(p_val)
            sig_str = 'Yes';
        end
        
        fprintf(fid_fstats, '%s & %.4f & %.4e & %s \\\\\n', ...
            problem_names{current_problems(prob_idx)}, chi2_val, p_val, sig_str);
    end
    
    fprintf(fid_fstats, '\\bottomrule\n');
    fprintf(fid_fstats, '\\end{tabular}\n');
    fprintf(fid_fstats, '\\end{table}\n');
    fclose(fid_fstats);
    fprintf('    Saved: %s\n', latex_file_friedman_stats);
end

% ==================== TABLE 9: Summary Statistics ====================
fprintf('  Generating Table 9: Summary Statistics...\n');

latex_file_summary = fullfile(latex_dir, 'Table9_Summary_Statistics.tex');
fid_summary = fopen(latex_file_summary, 'w');

if fid_summary ~= -1
    fprintf(fid_summary, '%%%% Table 9: Summary Statistics Over All Problems\n');
    fprintf(fid_summary, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_summary, '\\begin{table}[htbp]\n');
    fprintf(fid_summary, '\\centering\n');
    fprintf(fid_summary, '\\caption{Overall Summary Statistics Across All Problems (HARD VERSION - FAST)}\n');
    fprintf(fid_summary, '\\label{tab:summary}\n');
    fprintf(fid_summary, '\\small\n');
    fprintf(fid_summary, '\\begin{tabular}{lccccc}\n');
    fprintf(fid_summary, '\\toprule\n');
    fprintf(fid_summary, '\\textbf{Algorithm} & \\textbf{Avg Rank} & \\textbf{Friedman Rank} & \\textbf{Mean Time (s)} & \\textbf{Mean AUC} & \\textbf{W/T/L (Holm)} \\\\\n');
    fprintf(fid_summary, '\\midrule\n');
    
    [~, rank_sort_sum] = sort(overall_friedman_ranks);
    
    for si = 1:num_algorithms
        algo_idx = rank_sort_sum(si);
        
        % Calculate W/T/L against proposed using Holm-corrected p-values
        if algo_idx == proposed_idx
            w_str = '-'; t_str = '-'; l_str = '-';
        else
            % Find the correct index in other_indices
            a_idx = find(other_indices == algo_idx);
            w = 0; t = 0; l = 0;
            for prob_idx = 1:num_problems
                p_adj = P_holm(a_idx, prob_idx);
                if isnan(p_adj), continue; end
                if p_adj >= alpha, t = t + 1;
                else
                    kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
                    comp_mean = mean(squeeze(best_values(algo_idx, prob_idx, :)), 'omitnan');
                    if kdo_mean < comp_mean, w = w + 1; else, l = l + 1; end
                end
            end
            w_str = num2str(w); t_str = num2str(t); l_str = num2str(l);
        end
        
        if algo_idx == proposed_idx
            fprintf(fid_summary, '\\textbf{%s} & \\textbf{%.2f} & \\textbf{%.2f} & \\textbf{%.2f} & \\textbf{%.4f} & \\textbf{%s/%s/%s} \\\\\n', ...
                algorithms{algo_idx,1}, avg_ranks(algo_idx), overall_friedman_ranks(algo_idx), ...
                overall_mean_rt(algo_idx), overall_auc(algo_idx), w_str, t_str, l_str);
        else
            fprintf(fid_summary, '%s & %.2f & %.2f & %.2f & %.4f & %s/%s/%s \\\\\n', ...
                algorithms{algo_idx,1}, avg_ranks(algo_idx), overall_friedman_ranks(algo_idx), ...
                overall_mean_rt(algo_idx), overall_auc(algo_idx), w_str, t_str, l_str);
        end
    end
    
    fprintf(fid_summary, '\\bottomrule\n');
    fprintf(fid_summary, '\\end{tabular}\n');
    fprintf(fid_summary, '\\end{table}\n');
    fclose(fid_summary);
    fprintf('    Saved: %s\n', latex_file_summary);
end

% ==================== Create Master LaTeX File ====================
fprintf('  Generating Master LaTeX File...\n');

latex_master = fullfile(latex_dir, 'All_Tables_Master.tex');
fid_master = fopen(latex_master, 'w');

if fid_master ~= -1
    fprintf(fid_master, '%%%% Master LaTeX File for All Tables\n');
    fprintf(fid_master, '%%%% Electronics Optimization Results - HARD VERSION - FAST\n');
    fprintf(fid_master, '%%%% Generated: %s\n', datestr(now));
    fprintf(fid_master, '\\documentclass[12pt,a4paper]{article}\n');
    fprintf(fid_master, '\\usepackage[utf8]{inputenc}\n');
    fprintf(fid_master, '\\usepackage{amsmath}\n');
    fprintf(fid_master, '\\usepackage{amsfonts}\n');
    fprintf(fid_master, '\\usepackage{amssymb}\n');
    fprintf(fid_master, '\\usepackage{booktabs}\n');
    fprintf(fid_master, '\\usepackage{multirow}\n');
    fprintf(fid_master, '\\usepackage{rotating}\n');
    fprintf(fid_master, '\\usepackage{array}\n');
    fprintf(fid_master, '\\usepackage{adjustbox}\n');
    fprintf(fid_master, '\\usepackage[left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}\n');
    fprintf(fid_master, '\\begin{document}\n\n');
    
    fprintf(fid_master, '\\title{Electronics Optimization Results - Hard Version (Fast Track)}\n');
    fprintf(fid_master, '\\author{Optimization Benchmark Study}\n');
    fprintf(fid_master, '\\date{%s}\n', datestr(now));
    fprintf(fid_master, '\\maketitle\n\n');
    
    fprintf(fid_master, '\\section{Main Results}\n\n');
    fprintf(fid_master, '\\input{Table1_Main_Results}\n\n');
    
    fprintf(fid_master, '\\section{Detailed Statistics}\n\n');
    fprintf(fid_master, '\\input{Table2_Best_Median_Worst}\n\n');
    fprintf(fid_master, '\\input{Table7_P1_Detailed_Results}\n\n');
    fprintf(fid_master, '\\input{Table7_P2_Detailed_Results}\n\n');
    fprintf(fid_master, '\\input{Table7_P3_Detailed_Results}\n\n');
    fprintf(fid_master, '\\input{Table7_P4_Detailed_Results}\n\n');
    
    fprintf(fid_master, '\\section{Statistical Analysis}\n\n');
    fprintf(fid_master, '\\input{Table3_Friedman_Ranks}\n\n');
    fprintf(fid_master, '\\input{Table4_Wilcoxon_Advanced}\n\n');
    fprintf(fid_master, '\\input{Table8_Friedman_Statistics}\n\n');
    
    fprintf(fid_master, '\\section{Performance Metrics}\n\n');
    fprintf(fid_master, '\\input{Table5_Runtime_Analysis}\n\n');
    fprintf(fid_master, '\\input{Table6_AUC_Convergence}\n\n');
    
    fprintf(fid_master, '\\section{Overall Summary}\n\n');
    fprintf(fid_master, '\\input{Table9_Summary_Statistics}\n\n');
    
    fprintf(fid_master, '\\end{document}\n');
    fclose(fid_master);
    fprintf('    Saved: %s\n', latex_master);
end

fprintf('\nAll LaTeX tables saved in: %s\n', latex_dir);

%% Save MAT File (UPGRADED)
% =========================================================================
fprintf('\n========================================\n');
fprintf('SAVING MAT FILE\n');
fprintf('========================================\n\n');

if save_mat
    mat_file = fullfile(results_folder, 'Raw_Data', ...
        'Electronics_Results_Hard_Fast.mat');

    save(mat_file, ...
        'best_values', ...
        'convergence_curves', ...
        'execution_times', ...
        'auc_values', ...
        'best_positions', ...            % NEW
        'physical_results', ...          % NEW
        'param_defs', ...                % NEW
        'median_values', ...
        'mean_values', ...
        'std_values', ...
        'min_values', ...
        'max_values', ...
        'ranks', ...
        'friedman_ranks', ...
        'overall_friedman_ranks', ...
        'avg_ranks', ...
        'P_raw', 'P_holm', 'EffectSize', ...
        'friedman_stats', ...
        'perf_matrix', ...
        'rho', ...
        'tau_vals', ...
        'mean_runtime', ...
        'std_runtime', ...
        'mean_auc', ...
        'std_auc', ...
        'algorithms', ...
        'problem_names', ...
        'num_algorithms', ...
        'num_problems', ...
        'num_runs', ...
        'alpha', ...
        'CD', ...
        'proposed_idx');

    fprintf('MAT file saved: %s\n', mat_file);
end

%% Final Summary Display
% =========================================================================
total_time = toc(total_timer);

fprintf('\n');
fprintf('========================================\n');
fprintf('           FINAL SUMMARY (HARD VERSION - FAST)\n');
fprintf('========================================\n\n');

fprintf('Total execution time: %.2f seconds (%.2f minutes)\n', ...
    total_time, total_time/60);
fprintf('Results saved in: %s\n\n', results_folder);
fprintf('Output format: PNG (preview) + %s (journal)\n\n', upper(output_format));

fprintf('%-14s  %-8s  %-14s  %-10s\n', ...
    'Algorithm', 'Avg Rank', 'Friedman Rank', 'Mean Time(s)');
fprintf('%s\n', repmat('-', 1, 55));

[~, disp_sort] = sort(overall_friedman_ranks);

for si = 1:num_algorithms
    algo_idx = disp_sort(si);

    if algo_idx == proposed_idx
        marker = ' <-- PROPOSED';
    else
        marker = '';
    end

    fprintf('%-14s  %-8.3f  %-14.3f  %-10.3f%s\n', ...
        algorithms{algo_idx,1}, ...
        avg_ranks(algo_idx), ...
        overall_friedman_ranks(algo_idx), ...
        overall_mean_rt(algo_idx), ...
        marker);
end

fprintf('\n%s\n', repmat('=', 1, 55));
fprintf('Win/Tie/Loss Summary for %s (HARD VERSION - FAST, Holm-corrected):\n', algorithms{proposed_idx,1});
fprintf('%s\n', repmat('-', 1, 30));

for a = 1:length(other_indices)
    comp_idx = other_indices(a);
    w = 0; t = 0; l = 0;
    for prob_idx = 1:num_problems
        p_adj = P_holm(a, prob_idx);
        if isnan(p_adj), continue; end
        if p_adj >= alpha, t = t + 1;
        else
            kdo_mean = mean(squeeze(best_values(proposed_idx, prob_idx, :)), 'omitnan');
            comp_mean = mean(squeeze(best_values(comp_idx, prob_idx, :)), 'omitnan');
            if kdo_mean < comp_mean, w = w + 1; else, l = l + 1; end
        end
    end
    fprintf('vs %-10s : %d W / %d T / %d L\n', algorithms{comp_idx,1}, w, t, l);
end

fprintf('%s\n\n', repmat('=', 1, 55));

fprintf('Files generated:\n');
fprintf('  Excel   : %s\n', excel_filename);
fprintf('  LaTeX   : %s (and %d additional table files)\n', latex_master, 10);
fprintf('  MAT     : %s\n', mat_file);
fprintf('  Diary   : %s\n', diary_file);
fprintf('\nImages saved in:\n');
fprintf('  PNG (preview) : %s\n', fullfile(results_folder, 'Convergence_Curves'));
fprintf('  PDF (journal) : %s (same folders as PNG)\n', upper(output_format));
fprintf('\nDone.\n');

diary off;

%% ========================================================
%  NEW Helper Functions for Advanced Statistics
%  MUST BE ADDED TO THE END OF THE SCRIPT
%% ========================================================
function d = cliffs_delta(x, y)
    x = x(:); y = y(:);
    if numel(x) == 0 || numel(y) == 0, d = 0; return; end
    try diff_mat = sign(x - y'); catch, diff_mat = sign(bsxfun(@minus, x, y')); end
    d = sum(diff_mat(:)) / (numel(x) * numel(y));
end

function adj_p = holm_bonferroni(pvals)
    pvals = pvals(:);
    m = length(pvals);
    if m == 0, adj_p = []; return; end
    [sorted_p, idx] = sort(pvals, 'ascend');
    raw_adj = min(1, sorted_p .* (m:-1:1)');
    adj_sorted = raw_adj;
    for i = m-1:-1:1
        adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
    end
    adj_p = zeros(m,1);
    adj_p(idx) = adj_sorted;
end

%% ========================================================
%  Helper Function: save_figure_journal_safe (FIXED FOR R2017b)
%  Saves figures in PNG (preview) and PDF/EPS (journal format)
%% ========================================================
function save_figure_journal_safe(fig, filepath_base, output_format, save_png_flag)
    if ~ishandle(fig) || ~isvalid(fig)
        warning('save_figure_journal_safe: Invalid figure handle. Skipping save.');
        return;
    end
    try
        set(fig, 'Renderer', 'painters'); drawnow;
        if save_png_flag
            try, print(fig, [filepath_base '.png'], '-dpng', '-r300'); catch, end
        end
        if strcmp(output_format, 'pdf') || strcmp(output_format, 'both')
            try, set(fig, 'PaperPositionMode', 'auto'); print(fig, [filepath_base '.pdf'], '-dpdf', '-r300', '-bestfit'); catch, end
        end
        if strcmp(output_format, 'eps') || strcmp(output_format, 'both')
            try, print(fig, [filepath_base '.eps'], '-depsc', '-r300'); catch, end
        end
    catch
        warning('save_figure_journal_safe: Error saving figure.');
    end
end

%% ========================================================
%  Helper Function: writecell_safe
%  Compatible with ALL MATLAB versions (no writecell needed)
%% ========================================================
function writecell_safe(data, filename, sheetname)
    if isempty(data), return; end
    [nrows, ncols] = size(data);
    num_matrix  = NaN(nrows, ncols); str_matrix  = repmat({''}, nrows, ncols);
    has_numeric = false(nrows, ncols); has_string  = false(nrows, ncols);
    for r = 1:nrows
        for c = 1:ncols
            val = data{r, c};
            if isnumeric(val) && isscalar(val), num_matrix(r, c) = val; has_numeric(r, c) = true;
            elseif ischar(val) || islogical(val)
                if islogical(val), if val, str_matrix{r,c} = 'TRUE'; else, str_matrix{r,c} = 'FALSE'; end
                else, str_matrix{r, c} = val; end
                has_string(r, c) = true;
            elseif isnumeric(val) && isempty(val), num_matrix(r, c) = NaN; has_numeric(r, c) = true;
            else, str_matrix{r, c} = ''; has_string(r, c) = true;
            end
        end
    end
    try, xlswrite(filename, str_matrix, sheetname, 'A1'); catch, end
    num_rows = find(any(has_numeric, 2)); num_cols = find(any(has_numeric, 1));
    if ~isempty(num_rows) && ~isempty(num_cols)
        r_min = min(num_rows); r_max = max(num_rows); c_min = min(num_cols); c_max = max(num_cols);
        for r = r_min:r_max
            numeric_in_row = find(has_numeric(r, :));
            if ~isempty(numeric_in_row)
                blocks = {}; blk_start = numeric_in_row(1); blk_end = numeric_in_row(1);
                for ci = 2:length(numeric_in_row)
                    if numeric_in_row(ci) == blk_end + 1, blk_end = numeric_in_row(ci);
                    else, blocks{end+1} = [blk_start, blk_end]; blk_start = numeric_in_row(ci); blk_end = numeric_in_row(ci);
                    end
                end
                blocks{end+1} = [blk_start, blk_end];
                for bi = 1:length(blocks)
                    c1 = blocks{bi}(1); c2 = blocks{bi}(2); cell_ref = sprintf('%s%d', col_num_to_letter(c1), r);
                    try, xlswrite(filename, num_matrix(r, c1:c2), sheetname, cell_ref); catch, end
                end
            end
        end
    end
end

function letter = col_num_to_letter(col_num)
    letter = '';
    while col_num > 0
        remainder = mod(col_num - 1, 26);
        letter = [char(65 + remainder), letter];
        col_num = floor((col_num - 1) / 26);
    end
end