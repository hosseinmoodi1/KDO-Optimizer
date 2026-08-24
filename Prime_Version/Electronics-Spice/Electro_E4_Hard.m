function cost = Electro_E4_Hard(x)
% Electro_E4_Hard: 6T FinFET SRAM (22nm) Cost Function
% SILENT MODE - No debug output.

verbose = false;

    cfg = config();
    ngspice_exe = cfg.ngspice_exe;
    templates_dir = cfg.templates_dir;
    models_dir = cfg.models_dir;

    x = x(:)';

    temp_base = fullfile(cfg.scripts_dir, 'Outputs', 'temp');
    if ~exist(temp_base, 'dir'), mkdir(temp_base); end
    t = getCurrentTask();
    worker_id = 0; if ~isempty(t), worker_id = t.ID; end
    temp_dir = fullfile(temp_base, sprintf('worker_%d', worker_id));
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end

    param_file = fullfile(temp_dir, 'params_E4.inc');
    netlist_file = fullfile(temp_dir, 'run_E4.sp');
    log_file = fullfile(temp_dir, 'run_E4.log');

    lb = [1, 1, 1, 1, 20e-9, 20e-9, 20e-9, 0.6, 0.8e-9];
    ub = [8, 8, 8, 8, 40e-9, 40e-9, 40e-9, 1.0, 2.0e-9];
    center_point = (ub + lb) / 2;

    Wfin_pd = round(min(max(abs(x(1)), lb(1)), ub(1)));
    Wfin_pu = round(min(max(abs(x(2)), lb(2)), ub(2)));
    Wfin_ax = round(min(max(abs(x(3)), lb(3)), ub(3)));
    Wfin_ps = round(min(max(abs(x(4)), lb(4)), ub(4)));
    Lfin_pd = min(max(abs(x(5)), lb(5)), ub(5));
    Lfin_pu = min(max(abs(x(6)), lb(6)), ub(6));
    Lfin_ax = min(max(abs(x(7)), lb(7)), ub(7));
    VDD = min(max(x(8), lb(8)), ub(8));
    Tox = min(max(abs(x(9)), lb(9)), ub(9));

    x_clamped = [Wfin_pd, Wfin_pu, Wfin_ax, Wfin_ps, Lfin_pd, Lfin_pu, Lfin_ax, VDD, Tox];

    mu_n = 270e-4; mu_p = 80e-4;
    epsilon0 = 8.854e-12;
    epsilon_ox = 3.9;
    Cox = epsilon_ox * epsilon0 / Tox;
    W_single_fin = 40e-9;
    Wpd_eff = Wfin_pd * W_single_fin;
    Wpu_eff = Wfin_pu * W_single_fin;
    Wax_eff = Wfin_ax * W_single_fin;

    beta_pd = mu_n * Cox * (Wpd_eff / Lfin_pd);
    beta_pu = mu_p * Cox * (Wpu_eff / Lfin_pu);
    beta_ax = mu_n * Cox * (Wax_eff / Lfin_ax);

    CR = beta_pd / max(beta_ax, 1e-20);
    PR = beta_pu / max(beta_ax, 1e-20);

    SNM = (VDD / 2) * (1 - 1/(1+CR)) * (1 - PR/(1+PR));
    SNM = min(max(SNM, 0), VDD/2);

    fid = fopen(param_file, 'w');
    if fid == -1, cost = 1e10; return; end
    fprintf(fid, '* Parameters for E4 SIL run\n');
    fprintf(fid, '.param Wfin_pd=%.6e\n', Wfin_pd);
    fprintf(fid, '.param Wfin_pu=%.6e\n', Wfin_pu);
    fprintf(fid, '.param Wfin_ax=%.6e\n', Wfin_ax);
    fprintf(fid, '.param Wfin_ps=%.6e\n', Wfin_ps);
    fprintf(fid, '.param Lfin_pd=%.6e\n', Lfin_pd);
    fprintf(fid, '.param Lfin_pu=%.6e\n', Lfin_pu);
    fprintf(fid, '.param Lfin_ax=%.6e\n', Lfin_ax);
    fprintf(fid, '.param VDD=%.6e\n', VDD);
    fprintf(fid, '.param Tox=%.6e\n', Tox);
    fclose(fid);

    model_file = fullfile(models_dir, 'ptm22hp.lib');
    if ~isfile(model_file), cost = 1e10; return; end
    models_dir_norm = strrep(models_dir, '\', '/');
    param_file_norm = strrep(param_file, '\', '/');

    template = fileread(fullfile(templates_dir, 'template_E4.sp'));
    template = strrep(template, '<PARAM_FILE>', param_file_norm);
    template = strrep(template, '<MODELS_DIR>', models_dir_norm);

    fid = fopen(netlist_file, 'w');
    if fid == -1, cost = 1e10; return; end
    fprintf(fid, '%s', template);
    fclose(fid);

    if exist(log_file, 'file'), delete(log_file); end

    cmd = sprintf('"%s" -b -o "%s" "%s"', ngspice_exe, log_file, netlist_file);
    [status, ~] = system(cmd);

    if status ~= 0 || ~exist(log_file, 'file')
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist;
        cleanup_temp_files(param_file, netlist_file, log_file);
        return;
    end

    txt = fileread(log_file);
    cleanup_temp_files(param_file, netlist_file, log_file);

    function v = parse_scalar(text, pattern, defaultVal)
        tok = regexp(text, pattern, 'tokens', 'once');
        if isempty(tok)
            v = defaultVal;
            return;
        end
        v = str2double(tok{1});
        if ~isfinite(v)
            v = defaultVal;
        end
    end

    has_sentinel = ~isempty(regexpi(txt, 'E4_SIMULATION_DONE', 'once'));

    numPattern = '([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)';
    I_leak = parse_scalar(txt, ['(?m)^\s*I_leak_A\s*=\s*' numPattern '\s*$'], NaN);

    valid_sim = status == 0 && has_sentinel && isfinite(I_leak) && I_leak >= 0;

    if ~valid_sim
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; return;
    end

    yield = (SNM >= 0.15);
    Area_m2 = W_single_fin * (Wfin_pd * Lfin_pd + Wfin_pu * Lfin_pu + 2 * Wfin_ax * Lfin_ax);
    Area = Area_m2 * 1e12;

    snm_target = 0.15;
    cr_target = 1.2;
    pr_target = 0.8;
    leak_target = 10e-9;
    area_target = 100;

    snm_penalty = max(0, (snm_target - SNM) / snm_target);
    cr_penalty = max(0, (cr_target - CR) / cr_target);
    pr_penalty = max(0, (PR - pr_target) / pr_target);
    leak_penalty = I_leak / leak_target;
    area_penalty = Area / area_target;

    w_snm = 0.40;
    w_yield = 0.30;
    w_leak = 0.20;
    w_area = 0.10;

    yield_penalty = 1 - yield;

    cost = w_snm * snm_penalty + ...
           w_yield * yield_penalty + ...
           w_leak * leak_penalty + ...
           w_area * area_penalty;

    constraint_penalty = 1e3 * (snm_penalty^2 + cr_penalty^2 + pr_penalty^2 + leak_penalty^2);
    cost = cost + constraint_penalty;

    if ~isfinite(cost) || cost > 1e10
        cost = 1e10;
    end
end

function cleanup_temp_files(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end