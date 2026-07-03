%% Complete Benchmark Script for CEC (2017 or 2022)
% FULL VERSION with ALL Outputs: LaTeX Tables, Plots, Convergence Curves, Boxplots, etc.
% Compatible with MATLAB 2017b
% Generates PNG (preview) and PDF (for journal submission)
%
% VERSION: 8.0 - CEC ONLY (CEC2017 & CEC2022)
%   - 11 algorithms: KDO + 10 state-of-the-art competitors
%   - Correct dimensions: CEC2017 (10,30,50,100), CEC2022 (10,20 only)
%   - F2 removed from CEC2017 (unstable function)
%   - Proper MaxFEs-based stopping criterion
%   - Paired Wilcoxon Signed-Rank Test with Holm-Bonferroni Correction
%   - Effect Size (Cliff's Delta) calculation
%   - Friedman Test with average ranks
%   - Ablation Study support
%   - Parallel execution support
%   - NO Classic benchmark

clear; clc; close all;
warning off all;

%% ==================== CONFIGURATION ====================
benchmark_choice = 'CEC2017';  % 'CEC2017' or 'CEC2022'

% Main parameters
num_runs = 30;                 % Standard for high statistical power in optimization benchmarking

if strcmp(benchmark_choice, 'CEC2022')
    dim = 20;                  % CEC2022 max dimension is 20
elseif strcmp(benchmark_choice, 'CEC2017')
    dim = 30;                  % CEC2017 supports 10, 30, 50, 100
else
    error('Invalid benchmark choice. Use ''CEC2017'' or ''CEC2022''.');
end

enable_ablation = true;        % Enable Ablation Study
enable_parallel = false;       % Use parfor for speed (set to true if Parallel Toolbox available)
output_format = 'pdf';         % PDF is preferred by Nature Portfolio

% Generate outputs flags
generate_convergence_curves = true;
generate_boxplots = true;
generate_performance_profiles = true;
generate_cd_diagrams1 = true;
generate_ablation_results = true;

% CEC bias values
cec2017_bias = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, ...
                1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, ...
                2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000];
cec2022_bias = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200];

%% ==================== ALGORITHM LIST (11 ALGORITHMS) ====================
algorithms = {
    'KDO',      'Karma-Dharma Optimizer (Proposed)';
    'L_SHADE',  'L-SHADE Algorithm';
    'CMA_ES',   'Covariance Matrix Adaptation ES';
    'RUN',      'Runge-Kutta Optimizer';
    'GBO',      'Gradient-Based Optimizer';
    'DE',       'Differential Evolution';
    'HHO',      'Harris Hawks Optimization';
    'GWO',      'Grey Wolf Optimizer';
    'WOA',      'Whale Optimization Algorithm';
    'AVOA',     'African Vultures Optimization Algorithm';
    'COA',      'Crayfish Optimization Algorithm';
};
total_algorithms = size(algorithms, 1);
kdo_idx = find(strcmp(algorithms(:,1), 'KDO'), 1);
if isempty(kdo_idx), kdo_idx = 1; end

%% ==================== ABLATION STUDY ====================
if enable_ablation
    ablation_variants = {
        'KDO_NoNirvana', 'KDO without Nirvana Reset';
        'KDO_NoMemory',  'KDO without Cosmic Memory';
        'KDO_NoDharma',  'KDO without Dharma Phase';
    };
    all_algorithms = [algorithms; ablation_variants];
    total_algorithms_full = size(all_algorithms, 1);
else
    all_algorithms = algorithms;
    total_algorithms_full = total_algorithms;
end
kdo_idx = find(strcmp(all_algorithms(:,1), 'KDO'), 1);
if isempty(kdo_idx), kdo_idx = 1; end

%% ==================== BENCHMARK SETUP ====================
if strcmp(benchmark_choice, 'CEC2017')
    benchmark_func = @cec17_func; benchmark_name = 'CEC2017'; is_cec = true;
    bias_values = cec2017_bias; merged_func_order = [1, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 27, 29];
    functions_to_test = setdiff(1:29, 2);
    function_levels = cell(length(functions_to_test), 1);
    for i = 1:length(functions_to_test)
        f = functions_to_test(i);
        if ismember(f, [6, 9, 22, 27, 29]), function_levels{i} = 'Full';
        elseif ismember(f, [3, 4, 5, 10, 21]), function_levels{i} = 'Partial';
        else, function_levels{i} = 'Numerical'; end
    end
elseif strcmp(benchmark_choice, 'CEC2022')
    benchmark_func = @cec22_test_func; benchmark_name = 'CEC2022'; is_cec = true;
    bias_values = cec2022_bias; merged_func_order = 1:12;
    functions_to_test = 1:12;
    function_levels = cell(12, 1);
    for i = 1:12
        if ismember(i, [1, 2, 4, 5]), function_levels{i} = 'Full';
        elseif ismember(i, [3, 6, 9]), function_levels{i} = 'Partial';
        else, function_levels{i} = 'Numerical'; end
    end
else
    error('Invalid benchmark choice.');
end

% Create results folder
results_folder = [benchmark_name '_Analysis_D' num2str(dim) '_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
if ~exist(results_folder, 'dir'), mkdir(results_folder); end
folders = {'Convergence_Curves', 'Boxplots', 'Tables/LaTeX', 'Summary', 'Performance_Profiles', 'CD_Diagrams', 'Ablation_Results', 'Raw_Data'};
for i = 1:length(folders)
    if ~exist(fullfile(results_folder, folders{i}), 'dir'), mkdir(fullfile(results_folder, folders{i})); end
end

fprintf('\n================================================================================\n');
fprintf('Benchmark: %s (Dim=%d, Runs=%d, Algorithms=%d)\n', benchmark_name, dim, num_runs, total_algorithms_full);
fprintf('Functions to test: %d\n', length(functions_to_test));
fprintf('Output format for journal: %s\n', upper(output_format));
fprintf('Ablation Study: %s | Parallel Execution: %s\n', iff(enable_ablation, 'ENABLED', 'DISABLED'), iff(enable_parallel, 'ENABLED', 'DISABLED'));
if strcmp(benchmark_choice, 'CEC2022'), fprintf('NOTE: CEC2022 max dimension is 20. Current dim=%d is correct.\n', dim); end
fprintf('================================================================================\n');

%% ==================== VALIDATE FUNCTIONS ====================
fprintf('\nValidating functions...\n');
x_test = rand(dim, 1) * 200 - 100;
valid_functions = []; valid_levels = {};
for idx = 1:length(functions_to_test)
    func_num = functions_to_test(idx); level = function_levels{idx};
    try
        f = benchmark_func(x_test, func_num);
        if ~isnan(f) && ~isinf(f) && f < 1e30
            valid_functions = [valid_functions, func_num]; valid_levels{end+1} = level;
            fprintf('  F%d (%s): OK\n', func_num, level);
        end
    catch ME, fprintf('  F%d: Error - %s\n', func_num, ME.message); end
end
functions_to_test = valid_functions; function_levels = valid_levels;
fprintf('\nValid functions: %d\n', length(functions_to_test));

%% ==================== INITIALIZATION ====================
total_functions = length(functions_to_test);
all_best_values = zeros(total_functions, total_algorithms_full, num_runs);
all_runtimes = zeros(total_functions, total_algorithms_full, num_runs);
all_convergence = cell(total_functions, total_algorithms_full, num_runs);
rankings_per_function = zeros(total_functions, total_algorithms_full);
function_optimum = NaN(total_functions, 1);
if is_cec
    for f_idx = 1:length(functions_to_test), function_optimum(f_idx) = bias_values(functions_to_test(f_idx)); end
end

population_size = 50; MaxFEs = 10000 * dim; max_iterations = floor(MaxFEs / population_size);
fprintf('\nDimension: %d | Population: %d | MaxFEs: %d | Iterations: %d\n', dim, population_size, MaxFEs, max_iterations);

%% ==================== MAIN ANALYSIS LOOP ====================
total_time = tic;
if enable_parallel && exist('parfor', 'file')
    fprintf('\nStarting parallel execution...\n');
    all_best_values_temp = zeros(total_functions, total_algorithms_full, num_runs);
    all_runtimes_temp = zeros(total_functions, total_algorithms_full, num_runs);
    all_convergence_temp = cell(total_functions, total_algorithms_full, num_runs);
    parfor f_idx = 1:total_functions
        func_id = functions_to_test(f_idx);
        fobj = @(x) benchmark_func(x(:), func_id); lb = -100 * ones(1, dim); ub = 100 * ones(1, dim);
        f_best = zeros(num_runs, total_algorithms_full); f_time = zeros(num_runs, total_algorithms_full); f_conv = cell(num_runs, total_algorithms_full);
        for algo_idx = 1:total_algorithms_full
            algo_name = all_algorithms{algo_idx, 1};
            for run = 1:num_runs
                rng(run + func_id*1000 + algo_idx*100000, 'twister');
                try
                    tic_local = tic; [best, ~, conv] = feval(algo_name, population_size, max_iterations, lb, ub, dim, fobj, MaxFEs); elapsed = toc(tic_local);
                    f_best(run, algo_idx) = best; f_time(run, algo_idx) = elapsed;
                    if ~isempty(conv) && length(conv) > 1
                        if length(conv) < max_iterations, conv = [conv, best * ones(1, max_iterations - length(conv))]; end
                        f_conv{run, algo_idx} = conv;
                    end
                catch ME, f_best(run, algo_idx) = NaN; f_time(run, algo_idx) = NaN; end
            end
        end
        all_best_values_temp(f_idx, :, :) = f_best'; all_runtimes_temp(f_idx, :, :) = f_time'; all_convergence_temp(f_idx, :, :) = f_conv;
    end
    all_best_values = all_best_values_temp; all_runtimes = all_runtimes_temp; all_convergence = all_convergence_temp;
else
    fprintf('\nStarting sequential execution...\n');
    current_task = 0; total_tasks = total_functions * total_algorithms_full * num_runs;
    for f_idx = 1:total_functions
        func_id = functions_to_test(f_idx); level = function_levels{f_idx};
        fprintf('\n--- Function F%d [%s Analysis] ---\n', func_id, level);
        fobj = @(x) benchmark_func(x(:), func_id); lb = -100 * ones(1, dim); ub = 100 * ones(1, dim);
        func_best = zeros(num_runs, total_algorithms_full); func_time = zeros(num_runs, total_algorithms_full); func_convergence = cell(num_runs, total_algorithms_full);
        for algo_idx = 1:total_algorithms_full
            algo_name = all_algorithms{algo_idx, 1}; fprintf('  %s: [', algo_name);
            for run = 1:num_runs
                current_task = current_task + 1; if mod(run, 5) == 0, fprintf('.'); end
                rng(run + func_id*1000 + algo_idx*100000, 'twister');
                try
                    tic_local = tic; [best, ~, conv] = feval(algo_name, population_size, max_iterations, lb, ub, dim, fobj, MaxFEs); elapsed = toc(tic_local);
                    func_best(run, algo_idx) = best; func_time(run, algo_idx) = elapsed;
                    if ~isempty(conv) && length(conv) > 1
                        if length(conv) < max_iterations, conv = [conv, best * ones(1, max_iterations - length(conv))]; end
                        func_convergence{run, algo_idx} = conv; all_convergence{f_idx, algo_idx, run} = conv;
                    end
                catch ME, fprintf('\nError in %s, F%d, run %d: %s\n', algo_name, func_id, run, ME.message); func_best(run, algo_idx) = NaN; func_time(run, algo_idx) = NaN;
                end
            end
            fprintf('] %d%%\n', floor(100 * current_task / total_tasks));
            all_best_values(f_idx, algo_idx, :) = func_best(:, algo_idx); all_runtimes(f_idx, algo_idx, :) = func_time(:, algo_idx);
            valid = func_best(:, algo_idx); valid = valid(~isnan(valid) & ~isinf(valid));
            if ~isempty(valid), fprintf('     Best: %.4e | Median: %.4e | Mean: %.4e | Std: %.4e\n', min(valid), median(valid), mean(valid), std(valid)); end
        end
        % Calculate rankings
        medians = inf(1, total_algorithms_full);
        for algo_idx = 1:total_algorithms_full
            valid = func_best(:, algo_idx); valid = valid(~isnan(valid) & ~isinf(valid));
            if ~isempty(valid), medians(algo_idx) = median(valid); end
        end
        algo_with_data = find(~isinf(medians));
        if ~isempty(algo_with_data)
            [~, sort_idx] = sort(medians(algo_with_data));
            ranks = zeros(size(algo_with_data)); ranks(sort_idx) = 1:length(algo_with_data);
            [~, ~, ic] = unique(medians(algo_with_data));
            for i = 1:max(ic)
                tied = find(ic == i);
                if length(tied) > 1, avg_rank = mean(ranks(tied)); ranks(tied) = avg_rank; end
            end
            for i = 1:length(algo_with_data), rankings_per_function(f_idx, algo_with_data(i)) = ranks(i); end
        end
    end
end
total_elapsed = toc(total_time);
fprintf('\n\nTotal execution time: %.2f hours\n', total_elapsed/3600);

%% ==================== ADVANCED PAIRED STATISTICAL ANALYSIS ====================
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('ADVANCED PAIRED STATISTICAL ANALYSIS');
fprintf('\n%s\n', repmat('=', 1, 80));

alpha = 0.05; num_competitors = total_algorithms_full - 1;
win_count = 0; tie_count = 0; loss_count = 0;
P_raw = zeros(total_functions, num_competitors);
P_holm = zeros(total_functions, num_competitors);
EffectSize = zeros(total_functions, num_competitors);

fprintf('\nComputing paired tests, effect sizes, and Holm correction...\n');
for f = 1:total_functions
    kdo_results = squeeze(all_best_values(f, kdo_idx, :));
    raw_p_vec = zeros(num_competitors, 1);
    competitor_indices = setdiff(1:total_algorithms_full, kdo_idx);
    for a = 1:num_competitors
        comp_idx = competitor_indices(a);
        comp_results = squeeze(all_best_values(f, comp_idx, :));
        valid_mask = isfinite(kdo_results) & isfinite(comp_results);
        kdo_valid = kdo_results(valid_mask); comp_valid = comp_results(valid_mask);
        if length(kdo_valid) >= 5
            [p, ~] = signrank(kdo_valid, comp_valid, 'alpha', alpha, 'method', 'approximate');
            raw_p_vec(a) = p; P_raw(f, a) = p;
            EffectSize(f, a) = cliffs_delta(kdo_valid, comp_valid);
        else
            raw_p_vec(a) = NaN; P_raw(f, a) = NaN; EffectSize(f, a) = NaN;
        end
    end
    p_holm_vec = holm_bonferroni(raw_p_vec(~isnan(raw_p_vec)));
    P_holm(f, ~isnan(raw_p_vec)) = p_holm_vec;
end

fprintf('\n--- Summary of Paired Wilcoxon Signed-Rank Test vs KDO ---\n');
for a = 1:num_competitors
    comp_idx = competitor_indices(a); w = 0; t = 0; l = 0;
    fprintf('\n  vs %s:\n', all_algorithms{comp_idx,1});
    for f = 1:total_functions
        p_adj = P_holm(f, a); d = EffectSize(f, a);
        if isnan(p_adj)
            fprintf('    F%d: Insufficient data\n', functions_to_test(f));
        elseif p_adj >= alpha
            t = t + 1; fprintf('    F%d: No significant difference (p-adj=%.4f)\n', functions_to_test(f), p_adj);
        else
            kdo_mean = mean(squeeze(all_best_values(f, kdo_idx, :)), 'omitnan');
            comp_mean = mean(squeeze(all_best_values(f, comp_idx, :)), 'omitnan');
            if kdo_mean < comp_mean
                w = w + 1; fprintf('    F%d: KDO significantly better (p-adj=%.4f, delta=%.3f)\n', functions_to_test(f), p_adj, d);
            else
                l = l + 1; fprintf('    F%d: KDO significantly worse (p-adj=%.4f, delta=%.3f)\n', functions_to_test(f), p_adj, d);
            end
        end
    end
    win_count = win_count + w; tie_count = tie_count + t; loss_count = loss_count + l;
    fprintf('    Summary: %d W / %d T / %d L\n', w, t, l);
end
fprintf('\nOverall KDO performance: %d Wins / %d Ties / %d Losses (Holm-corrected)\n', win_count, tie_count, loss_count);

% Friedman Test
fprintf('\n--- Friedman Test (Overall Ranks) ---\n');
friedman_data = zeros(num_runs, total_algorithms_full);
for f_idx = 1:total_functions
    for algo_idx = 1:total_algorithms_full
        data = squeeze(all_best_values(f_idx, algo_idx, :)); friedman_data(:, algo_idx) = data(:);
    end
end
friedman_data(~isfinite(friedman_data)) = NaN;
valid_rows = all(~isnan(friedman_data), 2); friedman_data = friedman_data(valid_rows, :);
if size(friedman_data, 1) >= 2
    [p_friedman, ~, stats_friedman] = friedman(friedman_data, 1, 'off');
    fprintf('  Friedman test p-value: %.6f\n', p_friedman);
    [sorted_ranks, rank_idx] = sort(stats_friedman.meanranks);
    for i = 1:length(rank_idx)
        algo_idx = rank_idx(i);
        if algo_idx == kdo_idx, fprintf('    %-12s: %.3f  <-- PROPOSED\n', all_algorithms{algo_idx,1}, sorted_ranks(i));
        else, fprintf('    %-12s: %.3f\n', all_algorithms{algo_idx,1}, sorted_ranks(i)); end
    end
else
    stats_friedman = struct('meanranks', zeros(1, total_algorithms_full));
    fprintf('  Insufficient valid data for Friedman test\n');
end

%% ==================== GENERATE ALL OUTPUTS ====================
fprintf('\n%s\n', repmat('=', 1, 80)); fprintf('GENERATING ALL OUTPUTS'); fprintf('\n%s\n', repmat('=', 1, 80));
latex_dir = fullfile(results_folder, 'Tables', 'LaTeX');

fprintf('\n1. Generating LaTeX tables...\n');
generate_merged_latex_table(all_best_values, functions_to_test, function_levels, all_algorithms, benchmark_name, latex_dir, function_optimum, merged_func_order);
for f_idx = 1:length(functions_to_test)
    func_id = functions_to_test(f_idx); level = function_levels{f_idx}; is_full = strcmp(level, 'Full');
    generate_detailed_latex_table(all_best_values, all_runtimes, zeros(total_functions, total_algorithms_full), f_idx, func_id, level, is_full, all_algorithms, benchmark_name, latex_dir, function_optimum(f_idx));
end
generate_friedman_latex_table(rankings_per_function, function_levels, all_algorithms, benchmark_name, latex_dir, kdo_idx);

fprintf('\n1a. Generating enhanced statistical tables...\n');
generate_enhanced_latex_table(all_best_values, functions_to_test, all_algorithms, kdo_idx, P_holm, EffectSize, benchmark_name, latex_dir, dim);
generate_summary_latex_table(all_best_values, functions_to_test, all_algorithms, kdo_idx, P_holm, EffectSize, benchmark_name, latex_dir, stats_friedman);

if generate_convergence_curves, fprintf('\n2. Generating convergence curves...\n'); generate_all_convergence_curves(all_convergence, functions_to_test, function_levels, all_algorithms, benchmark_name, results_folder, output_format); end
if generate_boxplots, fprintf('\n3. Generating boxplots...\n'); generate_all_boxplots(all_best_values, functions_to_test, function_levels, all_algorithms, benchmark_name, results_folder, output_format); end
if generate_performance_profiles, fprintf('\n4. Generating performance profiles...\n'); generate_performance_profiles_plot(all_best_values, functions_to_test, all_algorithms, benchmark_name, results_folder, output_format); end
if generate_cd_diagrams1, fprintf('\n5. Generating CD diagrams...\n'); generate_cd_diagrams(rankings_per_function, function_levels, all_algorithms, benchmark_name, results_folder, output_format, kdo_idx); end
if enable_ablation && generate_ablation_results, fprintf('\n6. Generating Ablation Study results...\n'); generate_ablation_results_table(all_best_values, functions_to_test, function_levels, all_algorithms, kdo_idx, benchmark_name, results_folder); end
fprintf('\n7. Generating summary statistics...\n'); generate_summary_statistics(all_best_values, all_runtimes, rankings_per_function, functions_to_test, function_levels, all_algorithms, benchmark_name, results_folder, kdo_idx);

save(fullfile(results_folder, 'Raw_Data', 'raw_data.mat'), ...
     'all_best_values', 'all_runtimes', 'rankings_per_function', ...
     'function_levels', 'all_algorithms', 'benchmark_name', 'dim', 'num_runs', ...
     'win_count', 'tie_count', 'loss_count', 'kdo_idx', ...
     'P_raw', 'P_holm', 'EffectSize');

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE!\nResults saved in: %s\n', results_folder);
fprintf('========================================\n');
warning on all;

%% ==================== HELPER FUNCTIONS ====================
function s = iff(condition, true_val, false_val)
    if condition, s = true_val; else, s = false_val; end
end

function d = cliffs_delta(x, y)
    x = x(:); y = y(:);
    if numel(x) == 0 || numel(y) == 0, d = 0; return; end
    try diff_mat = sign(x - y'); catch, diff_mat = sign(bsxfun(@minus, x, y')); end
    d = sum(diff_mat(:)) / (numel(x) * numel(y));
end

function adj_p = holm_bonferroni(pvals)
    pvals = pvals(:); m = length(pvals);
    if m == 0, adj_p = []; return; end
    [sorted_p, idx] = sort(pvals, 'ascend');
    raw_adj = min(1, sorted_p .* (m:-1:1)');
    adj_sorted = raw_adj;
    for i = m-1:-1:1, adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1)); end
    adj_p = zeros(m,1); adj_p(idx) = adj_sorted;
end

%% ==================== OUTPUT GENERATION FUNCTIONS ====================
function save_figure_journal_format(fig, filepath_base, output_format)
    saveas(fig, [filepath_base '.png']); print(fig, [filepath_base '.png'], '-dpng', '-r300');
    if strcmp(output_format, 'pdf') || strcmp(output_format, 'both')
        set(fig, 'PaperPositionMode', 'auto'); print(fig, [filepath_base '.pdf'], '-dpdf', '-r300', '-bestfit');
    end
    if strcmp(output_format, 'eps') || strcmp(output_format, 'both'), print(fig, [filepath_base '.eps'], '-depsc', '-r300'); end
end

function generate_all_convergence_curves(all_convergence, functions_to_test, function_levels, algorithms, benchmark_name, results_folder, output_format)
    conv_dir = fullfile(results_folder, 'Convergence_Curves'); n_algorithms = size(algorithms, 1); n_functions = length(functions_to_test); colors = lines(n_algorithms);
    for f = 1:n_functions
        func_id = functions_to_test(f); level = function_levels{f}; fig = figure('Position', [100, 100, 900, 600], 'Visible', 'on'); hold on; legends = {};
        for algo = 1:n_algorithms
            mean_curve = []; count = 0;
            for run = 1:size(all_convergence, 3)
                curve = all_convergence{f, algo, run};
                if ~isempty(curve) && isvector(curve)
                    if isempty(mean_curve), mean_curve = zeros(size(curve)); end
                    mean_curve(1:length(curve)) = mean_curve(1:length(curve)) + curve; count = count + 1;
                end
            end
            if count > 0
                mean_curve = mean_curve / count;
                if min(mean_curve) > 0, semilogy(mean_curve, 'Color', colors(algo,:), 'LineWidth', 1.5); else, plot(mean_curve, 'Color', colors(algo,:), 'LineWidth', 1.5); end
                legends{end+1} = algorithms{algo, 1};
            end
        end
        if ~isempty(legends)
            xlabel('Iteration'); ylabel('Best Fitness'); title(sprintf('%s F%d (%s) - Convergence', benchmark_name, func_id, level)); legend(legends, 'Location', 'best'); grid on;
            filepath_base = fullfile(conv_dir, sprintf('%s_F%d_Convergence', benchmark_name, func_id)); save_figure_journal_format(fig, filepath_base, output_format);
        end
        close(fig);
    end
end

function generate_all_boxplots(all_best_values, functions_to_test, function_levels, algorithms, benchmark_name, results_folder, output_format)
    box_dir = fullfile(results_folder, 'Boxplots'); n_functions = length(functions_to_test); n_algorithms = size(algorithms, 1);
    for f = 1:n_functions
        func_id = functions_to_test(f); level = function_levels{f}; fig = figure('Position', [100, 100, 1000, 600], 'Visible', 'on'); data_matrix = []; group_labels = {};
        for algo = 1:n_algorithms
            values = squeeze(all_best_values(f, algo, :)); values = values(~isnan(values) & ~isinf(values));
            if ~isempty(values)
                data_matrix = [data_matrix; log10(abs(values) + 1e-10)]; group_labels = [group_labels, repmat({algorithms{algo,1}}, 1, length(values))];
            end
        end
        if ~isempty(data_matrix)
            boxplot(data_matrix, group_labels, 'Colors', lines(n_algorithms)); ylabel('log10(Fitness)'); title(sprintf('%s F%d (%s) - Boxplot', benchmark_name, func_id, level)); grid on; xtickangle(45);
            filepath_base = fullfile(box_dir, sprintf('%s_F%d_Boxplot', benchmark_name, func_id)); save_figure_journal_format(fig, filepath_base, output_format);
        end
        close(fig);
    end
end

function generate_performance_profiles_plot(all_best_values, functions_to_test, algorithms, benchmark_name, results_folder, output_format)
    profile_dir = fullfile(results_folder, 'Performance_Profiles'); n_algorithms = size(algorithms, 1); n_functions = length(functions_to_test); n_runs = size(all_best_values, 3);
    tau_values = logspace(-3, 3, 100); perf_profiles = zeros(length(tau_values), n_algorithms);
    for algo = 1:n_algorithms
        for t = 1:length(tau_values)
            tau = tau_values(t); count = 0; total = 0;
            for f = 1:n_functions
                for run = 1:n_runs
                    val = all_best_values(f, algo, run);
                    if ~isnan(val) && ~isinf(val)
                        total = total + 1; best_val = inf;
                        for a = 1:n_algorithms
                            v = all_best_values(f, a, run);
                            if ~isnan(v) && ~isinf(v) && v < best_val, best_val = v; end
                        end
                        if val <= best_val * tau, count = count + 1; end
                    end
                end
            end
            if total > 0, perf_profiles(t, algo) = count / total; end
        end
    end
    fig = figure('Position', [100, 100, 800, 600], 'Visible', 'on'); colors = lines(n_algorithms);
    for algo = 1:n_algorithms, semilogx(tau_values, perf_profiles(:, algo), 'Color', colors(algo,:), 'LineWidth', 2); hold on; end
    xlabel('\tau'); ylabel('Proportion Solved'); title(sprintf('%s - Performance Profile', benchmark_name)); legend(algorithms(:,1), 'Location', 'southeast'); grid on;
    filepath_base = fullfile(profile_dir, sprintf('%s_Performance_Profile', benchmark_name)); save_figure_journal_format(fig, filepath_base, output_format); close(fig);
end

function generate_cd_diagrams(rankings_per_function, function_levels, algorithms, benchmark_name, results_folder, output_format, kdo_idx)
    cd_dir = fullfile(results_folder, 'CD_Diagrams'); n_algorithms = size(algorithms, 1); fig = figure('Position', [100, 100, 800, 400], 'Visible', 'on');
    avg_ranks = mean(rankings_per_function, 1, 'omitnan'); [sorted_ranks, idx] = sort(avg_ranks); algo_names = algorithms(idx, 1);
    bar(sorted_ranks); set(gca, 'XTickLabel', algo_names, 'XTick', 1:n_algorithms); xtickangle(45); ylabel('Average Rank'); title(sprintf('%s - Overall CD Diagram', benchmark_name)); grid on;
    if ~isempty(kdo_idx)
        kdo_pos = find(idx == kdo_idx);
        if ~isempty(kdo_pos), hold on; bar(kdo_pos, sorted_ranks(kdo_pos), 'r'); hold off; end
    end
    filepath_base = fullfile(cd_dir, sprintf('%s_Overall_CD_Diagram', benchmark_name)); save_figure_journal_format(fig, filepath_base, output_format); close(fig);
end

function generate_ablation_results_table(all_best_values, functions_to_test, function_levels, algorithms, kdo_idx, benchmark_name, results_folder)
    ablation_dir = fullfile(results_folder, 'Ablation_Results');
    ablation_indices = find(contains(algorithms(:,1), 'No') | contains(algorithms(:,1), 'Fixed'));
    if isempty(ablation_indices), return; end
    n_ablation = length(ablation_indices); n_functions = length(functions_to_test); ablation_results = zeros(n_ablation + 1, n_functions);
    for f = 1:n_functions, vals = squeeze(all_best_values(f, kdo_idx, :)); vals = vals(isfinite(vals)); if ~isempty(vals), ablation_results(1, f) = median(vals); end; end
    for a = 1:n_ablation
        algo_idx = ablation_indices(a);
        for f = 1:n_functions, vals = squeeze(all_best_values(f, algo_idx, :)); vals = vals(isfinite(vals)); if ~isempty(vals), ablation_results(a+1, f) = median(vals); end; end
    end
    fid = fopen(fullfile(ablation_dir, 'Ablation_Study_Results.tex'), 'w');
    fprintf(fid, '\\begin{table}[htbp]\n\\centering\\small\n'); fprintf(fid, '\\caption{Ablation Study Results: Median Fitness Values}\n\\label{tab:ablation}\n');
    fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n\\textbf{Version}', repmat('c', 1, n_functions));
    for f = 1:n_functions, fprintf(fid, ' & \\textbf{F%d}', functions_to_test(f)); end; fprintf(fid, ' \\\\\n\\midrule\n');
    fprintf(fid, '\\textbf{KDO (Full)}'); for f = 1:n_functions, fprintf(fid, ' & %.2e', ablation_results(1, f)); end; fprintf(fid, ' \\\\\n');
    for a = 1:n_ablation
        algo_idx = ablation_indices(a); fprintf(fid, '%s', algorithms{algo_idx,1});
        for f = 1:n_functions, fprintf(fid, ' & %.2e', ablation_results(a+1, f)); end; fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n'); fclose(fid);
    fprintf('  Ablation Study table saved.\n');
end

function generate_summary_statistics(all_best_values, all_runtimes, rankings_per_function, functions_to_test, function_levels, algorithms, benchmark_name, results_folder, kdo_idx)
    summary_dir = fullfile(results_folder, 'Summary'); n_algorithms = size(algorithms, 1);
    overall_mean = zeros(1, n_algorithms); overall_std = zeros(1, n_algorithms); overall_rank = zeros(1, n_algorithms); overall_time = zeros(1, n_algorithms);
    for algo = 1:n_algorithms
        all_values = []; all_times = [];
        for f = 1:length(functions_to_test)
            values = squeeze(all_best_values(f, algo, :)); values = values(isfinite(values)); all_values = [all_values; values];
            times = squeeze(all_runtimes(f, algo, :)); times = times(isfinite(times)); all_times = [all_times; times];
        end
        if ~isempty(all_values), overall_mean(algo) = mean(all_values); overall_std(algo) = std(all_values); end
        overall_rank(algo) = mean(rankings_per_function(:, algo), 'omitnan');
        if ~isempty(all_times), overall_time(algo) = mean(all_times); end
    end
    fid = fopen(fullfile(summary_dir, [benchmark_name '_Summary.txt']), 'w');
    fprintf(fid, '%-12s %-12s %-12s %-12s %-12s\n', 'Algorithm', 'Mean', 'Std', 'Avg Rank', 'Time(s)');
    for algo = 1:n_algorithms
        if algo == kdo_idx, fprintf(fid, '%-12s %-12.4e %-12.4e %-12.4f %-12.4f  <-- PROPOSED\n', algorithms{algo,1}, overall_mean(algo), overall_std(algo), overall_rank(algo), overall_time(algo));
        else, fprintf(fid, '%-12s %-12.4e %-12.4e %-12.4f %-12.4f\n', algorithms{algo,1}, overall_mean(algo), overall_std(algo), overall_rank(algo), overall_time(algo)); end
    end
    fclose(fid);
end

function generate_enhanced_latex_table(all_best_values, functions_to_test, algorithms, kdo_idx, P_holm, EffectSize, benchmark_name, latex_dir, dim)
    competitor_indices = setdiff(1:size(algorithms,1), kdo_idx); num_competitors = length(competitor_indices); num_funcs = length(functions_to_test);
    fid = fopen(fullfile(latex_dir, sprintf('%s_Enhanced_Statistics.tex', benchmark_name)), 'w');
    fprintf(fid, '\\begin{table*}[!t]\n\\centering\\small\n');
    fprintf(fid, '\\caption{Comparative Results on %s (D=%d). Mean$\\pm$Std. ', benchmark_name, dim);
    fprintf(fid, '$^{\\dag}$=Holm-adj. $p<0.05$; $^{\\ddag}$=$p<0.01$. Effect size in brackets: S=Small, M=Medium, L=Large.}\n');
    fprintf(fid, '\\label{tab:enhanced}\n');
    fprintf(fid, '\\begin{tabular}{c|c%s}\n\\hline\n', repmat('c', 1, num_competitors));
    fprintf(fid, 'Func. & \\textbf{KDO}');
    for a = 1:num_competitors, fprintf(fid, ' & %s', algorithms{competitor_indices(a), 1}); end
    fprintf(fid, ' \\\\ \\hline\\hline\n');
    for f = 1:num_funcs
        fprintf(fid, 'F%d ', functions_to_test(f));
        kdo_mean = mean(squeeze(all_best_values(f, kdo_idx, :)), 'omitnan'); kdo_std = std(squeeze(all_best_values(f, kdo_idx, :)), 'omitnan');
        fprintf(fid, ' & %.2e$\\pm$%.2e ', kdo_mean, kdo_std);
        for a = 1:num_competitors
            comp_idx = competitor_indices(a); mean_val = mean(squeeze(all_best_values(f, comp_idx, :)), 'omitnan'); std_val = std(squeeze(all_best_values(f, comp_idx, :)), 'omitnan');
            p_adj = P_holm(f, a); d = EffectSize(f, a);
            if isnan(p_adj), fprintf(fid, ' & -- '); continue; end
            symbol = ''; eff_str = '';
            if p_adj < 0.01, symbol = '$^{\\ddag}$'; elseif p_adj < 0.05, symbol = '$^{\\dag}$'; end
            if p_adj < 0.05
                abs_d = abs(d);
                if abs_d >= 0.474, eff_str = '(L)'; elseif abs_d >= 0.33, eff_str = '(M)'; elseif abs_d >= 0.147, eff_str = '(S)'; else, eff_str = '(N)'; end
            end
            comp_mean = mean(squeeze(all_best_values(f, comp_idx, :)), 'omitnan');
            if kdo_mean < comp_mean && p_adj < 0.05, fprintf(fid, ' & \\textbf{%.2e$\\pm$%.2e}%s %s ', mean_val, std_val, symbol, eff_str);
            else, fprintf(fid, ' & %.2e$\\pm$%.2e%s %s ', mean_val, std_val, symbol, eff_str); end
        end
        fprintf(fid, '\\\\ \\hline\n');
    end
    fprintf(fid, '\\end{tabular}\n\\end{table*}\n'); fclose(fid);
    fprintf('  Enhanced statistical table (with KDO) saved.\n');
end

function generate_summary_latex_table(all_best_values, functions_to_test, algorithms, kdo_idx, P_holm, EffectSize, benchmark_name, latex_dir, stats_friedman)
    competitor_indices = setdiff(1:size(algorithms,1), kdo_idx); num_competitors = length(competitor_indices);
    fid = fopen(fullfile(latex_dir, sprintf('%s_Summary_Statistics.tex', benchmark_name)), 'w');
    fprintf(fid, '\\begin{table}[!t]\n\\centering\\small\n');
    fprintf(fid, '\\caption{Summary of Paired Wilcoxon Signed-Rank Test ($\\alpha=0.05$) with Holm-Bonferroni Correction}\n');
    fprintf(fid, '\\label{tab:summary_stats}\n');
    fprintf(fid, '\\begin{tabular}{l|c|c|c|c|c}\n\\hline\n');
    fprintf(fid, 'Algorithm & Wins & Ties & Losses & Avg $|\\delta|$ & Large Effect (L) \\\\\\hline\\hline\n');
    for a = 1:num_competitors
        comp_idx = competitor_indices(a); w = 0; t = 0; l = 0;
        for f = 1:length(functions_to_test)
            p_adj = P_holm(f, a); if isnan(p_adj), continue; end
            if p_adj >= 0.05, t = t + 1;
            else
                kdo_mean = mean(squeeze(all_best_values(f, kdo_idx, :)), 'omitnan'); comp_mean = mean(squeeze(all_best_values(f, comp_idx, :)), 'omitnan');
                if kdo_mean < comp_mean, w = w + 1; else, l = l + 1; end
            end
        end
        avg_d = mean(abs(EffectSize(:, a)), 'omitnan'); large_eff = sum(abs(EffectSize(:, a)) >= 0.474, 'omitnan');
        fprintf(fid, '%s & %d & %d & %d & %.3f & %d \\\\\\hline\n', algorithms{comp_idx,1}, w, t, l, avg_d, large_eff);
    end
    fprintf(fid, '\\end{tabular}\n\\end{table}\n'); fclose(fid);
    fprintf('  Summary statistical table saved.\n');
end

function generate_merged_latex_table(all_best_values, functions_to_test, function_levels, algorithms, benchmark_name, latex_dir, function_optimum, merged_func_order)
    n_algorithms = size(algorithms, 1); n_funcs = length(functions_to_test); mean_matrix = zeros(n_funcs, n_algorithms); std_matrix = zeros(n_funcs, n_algorithms);
    for f = 1:n_funcs
        for algo = 1:n_algorithms
            values = squeeze(all_best_values(f, algo, :)); values = values(isfinite(values));
            if ~isempty(values), mean_matrix(f, algo) = mean(values); std_matrix(f, algo) = std(values); else, mean_matrix(f, algo) = NaN; end
        end
    end
    fid = fopen(fullfile(latex_dir, sprintf('%s_Merged_Results.tex', benchmark_name)), 'w');
    fprintf(fid, '\\begin{sidewaystable}[p]\n\\centering\n\\caption{%s: Mean $\\pm$ Std results. Best results are highlighted in bold.}\n\\label{tab:%s_merged}\n\\resizebox{\\textheight}{!}{%%\n\\begin{tabular}{l%s}\n\\toprule\n\\textbf{Func}', benchmark_name, lower(benchmark_name), repmat('c', 1, n_algorithms));
    for algo = 1:n_algorithms, fprintf(fid, ' & \\textbf{%s}', algorithms{algo,1}); end; fprintf(fid, ' \\\\\n\\midrule\n');
    for f = 1:n_funcs
        func_id = functions_to_test(f); fprintf(fid, 'F%d', func_id);
        valid_means = mean_matrix(f, :); valid_means(isnan(valid_means)) = inf; [best_mean, best_algo] = min(valid_means);
        for algo = 1:n_algorithms
            if ~isnan(mean_matrix(f, algo))
                if algo == best_algo && best_mean < inf, fprintf(fid, ' & \\textbf{%.2e $\\pm$ %.2e}', mean_matrix(f, algo), std_matrix(f, algo));
                else, fprintf(fid, ' & %.2e $\\pm$ %.2e', mean_matrix(f, algo), std_matrix(f, algo)); end
            else, fprintf(fid, ' & --'); end
        end
        fprintf(fid, ' \\\\\n');
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}%%\n}\n\\end{sidewaystable}\n'); fclose(fid);
end

function generate_detailed_latex_table(all_best_values, all_runtimes, auc_normalized_results, f_idx, func_id, level, is_full, algorithms, benchmark_name, latex_dir, optimum)
    n_algorithms = size(algorithms, 1);
    best_vals = zeros(n_algorithms,1); mean_vals = zeros(n_algorithms,1); std_vals = zeros(n_algorithms,1); median_vals = zeros(n_algorithms,1); worst_vals = zeros(n_algorithms,1); time_vals = zeros(n_algorithms,1);
    for algo = 1:n_algorithms
        values = squeeze(all_best_values(f_idx, algo, :)); values = values(isfinite(values)); times = squeeze(all_runtimes(f_idx, algo, :)); times = times(isfinite(times));
        if ~isempty(values), best_vals(algo) = min(values); mean_vals(algo) = mean(values); std_vals(algo) = std(values); median_vals(algo) = median(values); worst_vals(algo) = max(values); time_vals(algo) = mean(times); else, mean_vals(algo) = NaN; end
    end
    valid_means = mean_vals; valid_means(isnan(valid_means)) = inf; [~, best_algo] = min(valid_means);
    fid = fopen(fullfile(latex_dir, sprintf('%s_F%d_Detailed.tex', benchmark_name, func_id)), 'w');
    fprintf(fid, '\\begin{table}[H]\n\\centering\\small\n\\caption{%s F%d (%s) - Detailed results. Best in bold.}\n\\label{tab:%s_F%d}\n', benchmark_name, func_id, level, lower(benchmark_name), func_id);
    if is_full, fprintf(fid, '\\begin{tabular}{lccccccc}\n\\toprule\n\\textbf{Algorithm} & \\textbf{Best} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Median} & \\textbf{Worst} & \\textbf{Time(s)} & \\textbf{AUC} \\\\\n\\midrule\n');
    else, fprintf(fid, '\\begin{tabular}{lcccccc}\n\\toprule\n\\textbf{Algorithm} & \\textbf{Best} & \\textbf{Mean} & \\textbf{Std} & \\textbf{Median} & \\textbf{Worst} & \\textbf{Time(s)} \\\\\n\\midrule\n'); end
    for algo = 1:n_algorithms
        if ~isnan(mean_vals(algo))
            if algo == best_algo
                if is_full, fprintf(fid, '%s & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.3f} & \\textbf{%.3f} \\\\\n', algorithms{algo,1}, best_vals(algo), mean_vals(algo), std_vals(algo), median_vals(algo), worst_vals(algo), time_vals(algo), 0);
                else, fprintf(fid, '%s & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.2e} & \\textbf{%.3f} \\\\\n', algorithms{algo,1}, best_vals(algo), mean_vals(algo), std_vals(algo), median_vals(algo), worst_vals(algo), time_vals(algo)); end
            else
                if is_full, fprintf(fid, '%s & %.2e & %.2e & %.2e & %.2e & %.2e & %.3f & %.3f \\\\\n', algorithms{algo,1}, best_vals(algo), mean_vals(algo), std_vals(algo), median_vals(algo), worst_vals(algo), time_vals(algo), 0);
                else, fprintf(fid, '%s & %.2e & %.2e & %.2e & %.2e & %.2e & %.3f \\\\\n', algorithms{algo,1}, best_vals(algo), mean_vals(algo), std_vals(algo), median_vals(algo), worst_vals(algo), time_vals(algo)); end
            end
        end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n'); fclose(fid);
end

function generate_friedman_latex_table(rankings_per_function, function_levels, algorithms, benchmark_name, latex_dir, kdo_idx)
    n_algorithms = size(algorithms, 1); overall_ranks = mean(rankings_per_function, 1, 'omitnan'); [~, order] = sort(overall_ranks);
    fid = fopen(fullfile(latex_dir, sprintf('%s_Friedman_Ranks.tex', benchmark_name)), 'w');
    fprintf(fid, '\\begin{table}[H]\n\\centering\\small\n\\caption{Friedman average ranks for %s}\n\\label{tab:%s_friedman}\n\\begin{tabular}{lc}\n\\toprule\n\\textbf{Algorithm} & \\textbf{Overall Rank} \\\\\n\\midrule\n', benchmark_name, lower(benchmark_name));
    for i = 1:length(order)
        algo_idx = order(i);
        if algo_idx == kdo_idx, fprintf(fid, '\\textbf{%s} & \\textbf{%.3f}  <-- PROPOSED \\\\\n', algorithms{algo_idx,1}, overall_ranks(algo_idx));
        else, fprintf(fid, '%s & %.3f \\\\\n', algorithms{algo_idx,1}, overall_ranks(algo_idx)); end
    end
    fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n'); fclose(fid);
end