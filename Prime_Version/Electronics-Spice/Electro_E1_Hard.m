function cost = Electro_E1_Hard(x)
% Electro_E1_Hard: Two-Stage Miller-Compensated Op-Amp (180nm)
% FINAL VERSION: SILENT MODE - No debug output.

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

    param_file = fullfile(temp_dir, 'params_E1.inc');
    netlist_file = fullfile(temp_dir, 'run_E1.sp');
    log_file = fullfile(temp_dir, 'run_E1.log');

    ub = [100e-6 * ones(1,8), 1e-6 * ones(1,8), 500e-6, 10e-12];
    lb = [1e-6 * ones(1,8), 200e-9 * ones(1,8), 1e-6, 0.1e-12];
    center_point = (ub + lb) / 2;

    W = min(max(abs(x(1:8)), lb(1:8)), ub(1:8));
    L = min(max(abs(x(9:16)), lb(9:16)), ub(9:16));
    Ibias = min(max(abs(x(17)), lb(17)), ub(17));
    Cc = min(max(abs(x(18)), lb(18)), ub(18));

    VDD = 1.8; VTHN = 0.45; VTHP = -0.45;
    UN_Cox = 180e-6; UP_Cox = 60e-6;
    lambda_n = 0.05; lambda_p = 0.05;
    I_tail = Ibias; I_d = I_tail / 2; I_stage2 = I_tail / 2;
    gm1 = sqrt(2 * UN_Cox * (W(1)/L(1)) * I_d);
    gm6 = sqrt(2 * UN_Cox * (W(6)/L(6)) * I_stage2);
    ro1 = 1 / (lambda_n * I_d + eps); ro3 = 1 / (lambda_p * I_d + eps);
    ro6 = 1 / (lambda_n * I_stage2 + eps); ro7 = 1 / (lambda_p * I_stage2 + eps);
    Rout1 = ro1 * ro3 / (ro1 + ro3 + eps); Rout2 = ro6 * ro7 / (ro6 + ro7 + eps);
    Gain_est = gm1 * Rout1 * gm6 * Rout2;
    GBW_est = gm1 / (2 * pi * Cc + eps);
    PM_est = 90 - atand(GBW_est / (gm6 / (2 * pi * Cc) + eps));

    if Gain_est < 1 || GBW_est < 1e3 || PM_est < 1
        dist = sum(((x - center_point) ./ (ub - lb)).^2);
        cost = 1e9 + 1e6 * dist; cost = cost(1); return;
    end

    fid = fopen(param_file, 'w');
    if fid == -1, cost = 1e9; return; end
    fprintf(fid, '* Parameters for E1 SIL run\n');
    for i = 1:8
        fprintf(fid, '.param W%d=%.6e\n', i, W(i));
        fprintf(fid, '.param L%d=%.6e\n', i, L(i));
    end
    fprintf(fid, '.param Ibias=%.6e\n', Ibias);
    fprintf(fid, '.param Cc=%.6e\n', Cc);
    fclose(fid);

    model_file = fullfile(models_dir, 'ptm180.lib');
    if ~isfile(model_file), cost = 1e9; return; end
    models_dir_spice = strrep(models_dir, '\', '/');
    param_file_spice = strrep(param_file, '\', '/');

    template = fileread(fullfile(templates_dir, 'template_E1.sp'));
    template = strrep(template, '<PARAM_FILE>', param_file_spice);
    template = strrep(template, '<MODELS_DIR>', models_dir_spice);

    fid = fopen(netlist_file, 'w');
    if fid == -1, cost = 1e9; return; end
    fprintf(fid, '%s', template);
    fclose(fid);

    if exist(log_file, 'file'), delete(log_file); end

    cmd = sprintf('"%s" -b -o "%s" "%s"', ngspice_exe, log_file, netlist_file);
    [status, ~] = system(cmd);

    if status ~= 0 || ~exist(log_file, 'file')
        dist = sum(((x - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; cost = cost(1);
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

    Gain_db = parse_scalar(txt, '(?i)Gain_low_db\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    if isnan(Gain_db)
        Gain_db = parse_scalar(txt, '(?i)Gain_max_db\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    end

    GBW = parse_scalar(txt, '(?i)GBW\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    PM  = parse_scalar(txt, '(?i)PM\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    Power = parse_scalar(txt, '(?i)Power_W\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    SR = I_tail / Cc;

    Vout_dc = parse_scalar(txt, '(?i)Vout_dc\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    Vn2_dc = parse_scalar(txt, '(?i)Vn2_dc\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);

    has_ac_data = ~isempty(regexp(txt, 'No\. of Data Rows\s*:\s*\d+', 'once'));
    has_measure_gain = ~isnan(Gain_db);
    simulation_fatal = ~isempty(regexpi(txt, 'fatal|singular matrix|no convergence', 'once'));

    if simulation_fatal || status ~= 0 || ~has_ac_data || ~has_measure_gain
        dist = sum(((x - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; cost = cost(1); return;
    end

    gain_target = 74;
    gbw_target = 5e6;
    pm_target = 60;
    sr_target = 10e6;
    power_target = 5e-3;

    bias_penalty = 0;
    rail_margin = 0.1;
    
    if ~isnan(Vout_dc)
        if Vout_dc < rail_margin || Vout_dc > (1.8 - rail_margin)
            dist_rail = max(rail_margin - Vout_dc, Vout_dc - (1.8 - rail_margin));
            bias_penalty = bias_penalty + 1e4 * (1 + max(0, dist_rail) * 50);
        end
    end
    
    if ~isnan(Vn2_dc)
        if Vn2_dc < rail_margin || Vn2_dc > (1.8 - rail_margin)
            dist_rail = max(rail_margin - Vn2_dc, Vn2_dc - (1.8 - rail_margin));
            bias_penalty = bias_penalty + 1e4 * (1 + max(0, dist_rail) * 50);
        end
    end

    sym_penalty = 0;
    ratio1 = (W(1)/L(1)) / (W(2)/L(2) + eps);
    ratio3 = (W(3)/L(3)) / (W(4)/L(4) + eps);
    ratio7_8 = (W(7)/L(7)) / (W(8)/L(8) + eps);
    sym_penalty = sym_penalty + 1e3 * abs(log(ratio1))^2;
    sym_penalty = sym_penalty + 1e3 * abs(log(ratio3))^2;
    sym_penalty = sym_penalty + 1e2 * abs(log(ratio7_8))^2;

    gbw_failed = (isnan(GBW) || GBW <= 0);
    pm_failed = (isnan(PM) || PM <= -180);

    if gbw_failed || pm_failed
        gain_deficit = max(0, (gain_target - Gain_db) / gain_target);
        if Gain_db > 0, gain_deficit = gain_deficit * 0.5; end
        power_penalty = max(0, (Power - power_target) / power_target);
        sr_penalty = max(0, (sr_target - SR) / sr_target);
        cost = 1e6 + 1e4 * gain_deficit^2 + 1e3 * power_penalty + 1e3 * sr_penalty ...
               + bias_penalty + sym_penalty;
        cost = cost(1);
        return;
    end

    gain_penalty = max(0, (gain_target - Gain_db) / gain_target);
    gbw_penalty  = max(0, (gbw_target - GBW) / gbw_target);
    pm_penalty   = max(0, (pm_target - PM) / pm_target);
    sr_penalty   = max(0, (sr_target - SR) / sr_target);
    power_penalty = max(0, (Power - power_target) / power_target);

    w_gain = 0.25; w_gbw = 0.25; w_pm = 0.20; w_sr = 0.10; w_power = 0.10;
    base_cost = w_gain * (gain_target / max(Gain_db, 1)) + ...
                w_gbw  * (gbw_target / max(GBW, 1e3)) + ...
                w_pm   * (pm_target / max(PM, 1)) + ...
                w_sr   * (sr_target / max(SR, 1e3)) + ...
                w_power * (Power / power_target);

    total_penalty = gain_penalty + gbw_penalty + pm_penalty + sr_penalty + power_penalty;
    cost = base_cost + 1e4 * total_penalty + 0.05 * bias_penalty + 0.05 * sym_penalty;

    if ~isfinite(cost) || cost > 1e10, cost = 1e10; end
end

function cleanup_temp_files(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end