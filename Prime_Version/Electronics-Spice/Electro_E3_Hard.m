function cost = Electro_E3_Hard(x)
% Electro_E3_Hard: LC-VCO (65nm) Cost Function
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

    tag = sprintf('E3_w%d_%d', worker_id, randi(1e9));
    param_file = fullfile(temp_dir, [tag '_params.inc']);
    netlist_file = fullfile(temp_dir, [tag '.sp']);
    log_file = fullfile(temp_dir, [tag '.log']);

    ub = [20e-9, 2e-12, 500e-15, 15e-6, 25e-6, 10e-3, 10e3, 200e-15, 1.2, 1.2, 10e9];
    lb = [0.5e-9, 50e-15, 10e-15, 1e-6, 1e-6, 0.5e-3, 100, 5e-15, 0, 0.1, 50e6];
    center_point = (ub + lb) / 2;

    L = min(max(x(1), lb(1)), ub(1));
    Cvar = min(max(x(2), lb(2)), ub(2));
    Cfix = min(max(x(3), lb(3)), ub(3));
    Wn = min(max(x(4), lb(4)), ub(4));
    Wp = min(max(x(5), lb(5)), ub(5));
    Ibias = min(max(x(6), lb(6)), ub(6));
    Rbias = min(max(x(7), lb(7)), ub(7));
    Cbuf = min(max(x(8), lb(8)), ub(8));

    min_ctrl_span = 0.2;
    Vctrl_min = min(max(x(9), lb(9)), ub(10) - min_ctrl_span);
    Vctrl_max = min(max(x(10), Vctrl_min + min_ctrl_span), ub(10));
    K_vco = min(max(x(11), lb(11)), ub(11));

    x_clamped = [L, Cvar, Cfix, Wn, Wp, Ibias, Rbias, Cbuf, Vctrl_min, Vctrl_max, K_vco];

    Ceq_approx = Cvar + 0.5 * (Cfix + Cbuf);
    f_approx = 1 / (2 * pi * sqrt(L * max(Ceq_approx, 1e-24)));
    
    if f_approx < 0.5e9 || f_approx > 6e9
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; return;
    end
    
    if Ibias < 0.8e-3 || Wn < 1.5e-6 || Wp < 1.5e-6
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; return;
    end

    Q_est = 5;
    Rtank_est = Q_est * 2 * pi * f_approx * L;
    Vamp_est = Ibias * Rtank_est / 2;
    max_safe_amplitude = 1.0;
    
    if Vamp_est > max_safe_amplitude
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; return;
    end

    fid = fopen(param_file, 'w');
    if fid == -1, cost = 1e12; return; end
    fprintf(fid, '* Parameters for E3\n');
    fprintf(fid, '.param L=%.6e\n', L);
    fprintf(fid, '.param Cvar=%.6e\n', Cvar);
    fprintf(fid, '.param Cfix=%.6e\n', Cfix);
    fprintf(fid, '.param Wn=%.6e\n', Wn);
    fprintf(fid, '.param Wp=%.6e\n', Wp);
    fprintf(fid, '.param Ibias=%.6e\n', Ibias);
    fprintf(fid, '.param Rbias=%.6e\n', Rbias);
    fprintf(fid, '.param Cbuf=%.6e\n', Cbuf);
    fprintf(fid, '.param Vctrl_min=%.6e\n', Vctrl_min);
    fprintf(fid, '.param Vctrl_max=%.6e\n', Vctrl_max);
    fprintf(fid, '.param K_vco=%.6e\n', K_vco);
    fclose(fid);

    model_file = fullfile(models_dir, 'ptm65.lib');
    if ~isfile(model_file), cost = 1e12; return; end
    models_dir_norm = strrep(models_dir, '\', '/');
    param_file_norm = strrep(param_file, '\', '/');

    template_file = fullfile(templates_dir, 'template_E3.sp');
    if ~isfile(template_file), cost = 1e12; return; end

    template = fileread(template_file);
    template = strrep(template, '<PARAM_FILE>', param_file_norm);
    template = strrep(template, '<MODELS_DIR>', models_dir_norm);

    fid = fopen(netlist_file, 'w');
    if fid == -1, cost = 1e12; return; end
    fprintf(fid, '%s', template);
    fclose(fid);

    if exist(log_file, 'file'), delete(log_file); end

    cmd = sprintf('"%s" -b -o "%s" "%s" > nul 2>&1', ngspice_exe, log_file, netlist_file);
    [status, ~] = system(cmd);

    if status ~= 0 || ~exist(log_file, 'file')
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist;
        cleanup_temp_files(param_file, netlist_file, log_file);
        return;
    end

    txt = fileread(log_file);
    cleanup_temp_files(param_file, netlist_file, log_file);

    numPattern = '([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)';
    
    tok = regexpi(txt, ['tperiod\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok)
        tperiod = str2double(tok{1});
        if isfinite(tperiod) && tperiod > 0
            f_osc_Hz = 1.0 / tperiod;
        else
            f_osc_Hz = NaN;
        end
    else
        f_osc_Hz = NaN;
    end

    tok = regexpi(txt, ['Power_W\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok)
        Power = str2double(tok{1});
        if ~isfinite(Power) || Power < 0, Power = NaN; end
    else
        Power = NaN;
    end

    tok = regexpi(txt, ['vpp_out\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok)
        Vpp_out = str2double(tok{1});
        if ~isfinite(Vpp_out) || Vpp_out < 0, Vpp_out = NaN; end
    else
        Vpp_out = NaN;
    end

    tok = regexpi(txt, ['voutp_max\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok), voutp_max = str2double(tok{1}); else, voutp_max = NaN; end
    
    tok = regexpi(txt, ['voutp_min\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok), voutp_min = str2double(tok{1}); else, voutp_min = NaN; end
    
    tok = regexpi(txt, ['voutn_max\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok), voutn_max = str2double(tok{1}); else, voutn_max = NaN; end
    
    tok = regexpi(txt, ['voutn_min\s*=\s*' numPattern], 'tokens', 'once');
    if ~isempty(tok), voutn_min = str2double(tok{1}); else, voutn_min = NaN; end

    has_sentinel = ~isempty(regexpi(txt, 'E3_SIMULATION_DONE', 'once'));
    has_tperiod = isfinite(f_osc_Hz) && f_osc_Hz > 0;
    fatal = ~isempty(regexpi(txt, 'fatal|singular matrix|no convergence', 'once'));
    rails_measured = all(isfinite([voutp_max, voutp_min, voutn_max, voutn_min]));

    if has_tperiod && isfinite(Power) && Power > 0
        PN_dBc = -120 - 10*log10(Power/1e-3) - 20*log10(f_osc_Hz/1e9);
    else
        PN_dBc = NaN;
    end

    f_min = 1 / (2 * pi * sqrt(L * (Cvar + Cfix)));
    f_max = 1 / (2 * pi * sqrt(L * (0.3*Cvar + Cfix)));
    TR = (f_max - f_min) / mean([f_max, f_min]) * 100;

    VDD = 1.2;
    margin = 0.03 * VDD;
    rail_violation = 0;
    if rails_measured
        upper_violation = max(0, voutp_max - (VDD + margin)) + ...
                          max(0, voutn_max - (VDD + margin));
        lower_violation = max(0, -margin - voutp_min) + ...
                          max(0, -margin - voutn_min);
        rail_violation = upper_violation + lower_violation;
    end
    severe_rail_violation = rail_violation > 0.1;

    valid_sim = status == 0 && has_sentinel && has_tperiod && ~fatal && ...
                isfinite(Power) && Power > 0 && ...
                isfinite(Vpp_out) && Vpp_out > 1e-3 && ...
                rails_measured && ~severe_rail_violation;

    if ~valid_sim
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; return;
    end

    f_target = 2.5e9;
    pn_target = -120;
    power_target = 8e-3;
    tr_target = 10;

    f_penalty = ((f_osc_Hz - f_target) / f_target)^2;
    freq_ok = f_osc_Hz >= 2.4e9 && f_osc_Hz <= 2.6e9;
    
    if ~freq_ok
        freq_range_penalty = 10 * (max(0, (2.4e9 - f_osc_Hz) / 1e9) + ...
                                   max(0, (f_osc_Hz - 2.6e9) / 1e9))^2;
    else
        freq_range_penalty = 0;
    end

    pn_penalty = max(0, (PN_dBc - pn_target) / 10);
    power_score = Power / power_target;
    power_violation = max(0, (Power - power_target) / power_target);
    power_penalty = 0.5 * power_score + 0.5 * power_violation;

    tr_penalty = max(0, (tr_target - TR) / tr_target);
    tr_score = tr_target / max(TR, 1e-6);
    tr_penalty_combined = 0.1 * tr_score + 0.9 * tr_penalty;

    w_f = 0.30; w_pn = 0.25; w_power = 0.20; w_tr = 0.20; w_freq_range = 0.05;

    cost = w_f * f_penalty + w_pn * pn_penalty + w_power * power_penalty + ...
           w_tr * tr_penalty_combined + w_freq_range * freq_range_penalty;

    rail_penalty = 1e3 * rail_violation^2;
    constraint_penalty = 1e3 * (pn_penalty^2 + power_violation^2 + tr_penalty^2);
    cost = cost + rail_penalty + constraint_penalty;

    if ~isfinite(cost) || cost > 1e10, cost = 1e10; end
end

function cleanup_temp_files(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end