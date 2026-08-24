function cost = Electro_E2_Hard(x)
% Electro_E2_Hard: RF Cascode LNA (65nm)
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

    param_file = fullfile(temp_dir, 'params_E2.inc');
    netlist_file = fullfile(temp_dir, 'run_E2.sp');
    log_file = fullfile(temp_dir, 'run_E2.log');

    ub = [50e-6*ones(1,4), 20e-9*ones(1,3), 500e-15*ones(1,3), 5e-3, 1.8];
    lb = [0.5e-6*ones(1,4), 0.5e-9*ones(1,3), 10e-15*ones(1,3), 0.1e-3, 0.5];
    center_point = (ub + lb) / 2;

    Wm1 = min(max(abs(x(1)), lb(1)), ub(1));
    Wm2 = min(max(abs(x(2)), lb(2)), ub(2));
    Wm3 = min(max(abs(x(3)), lb(3)), ub(3));
    Wm4 = min(max(abs(x(4)), lb(4)), ub(4));
    Lg1 = min(max(abs(x(5)), lb(5)), ub(5));
    Lg2 = min(max(abs(x(6)), lb(6)), ub(6));
    Ld  = min(max(abs(x(7)), lb(7)), ub(7));
    Cgd1 = min(max(abs(x(8)), lb(8)), ub(8));
    Cgd2 = min(max(abs(x(9)), lb(9)), ub(9));
    Cgd3 = min(max(abs(x(10)), lb(10)), ub(10));
    Id = min(max(abs(x(11)), lb(11)), ub(11));
    Vdd = min(max(x(12), lb(12)), ub(12));

    x_clamped = [Wm1, Wm2, Wm3, Wm4, Lg1, Lg2, Ld, Cgd1, Cgd2, Cgd3, Id, Vdd];

    if Id < 0.5e-3 || Vdd < 0.5 || Wm1 < 0.5e-6 || Wm2 < 0.5e-6
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e9 + 1e4 * dist; cost = cost(1); return;
    end

    fid = fopen(param_file, 'w');
    if fid == -1, cost = 1e9; return; end
    fprintf(fid, '* Parameters for E2 SIL run\n');
    fprintf(fid, '.param Wm1=%.6e\n', Wm1);
    fprintf(fid, '.param Wm2=%.6e\n', Wm2);
    fprintf(fid, '.param Wm3=%.6e\n', Wm3);
    fprintf(fid, '.param Wm4=%.6e\n', Wm4);
    fprintf(fid, '.param Lg1=%.6e\n', Lg1);
    fprintf(fid, '.param Lg2=%.6e\n', Lg2);
    fprintf(fid, '.param Ld=%.6e\n', Ld);
    fprintf(fid, '.param Cgd1=%.6e\n', Cgd1);
    fprintf(fid, '.param Cgd2=%.6e\n', Cgd2);
    fprintf(fid, '.param Cgd3=%.6e\n', Cgd3);
    fprintf(fid, '.param Id=%.6e\n', Id);
    fprintf(fid, '.param Vdd=%.6e\n', Vdd);
    fclose(fid);

    model_file = fullfile(models_dir, 'ptm65.lib');
    if ~isfile(model_file), cost = 1e9; return; end
    model_file_spice = strrep(model_file, '\', '/');
    param_file_spice = strrep(param_file, '\', '/');

    template = fileread(fullfile(templates_dir, 'template_E2.sp'));
    template = strrep(template, '<PARAM_FILE>', param_file_spice);
    template = strrep(template, '<MODEL_FILE>', model_file_spice);

    fid = fopen(netlist_file, 'w');
    if fid == -1, cost = 1e9; return; end
    fprintf(fid, '%s', template);
    fclose(fid);

    if exist(log_file, 'file'), delete(log_file); end

    cmd = sprintf('"%s" -b -o "%s" "%s"', ngspice_exe, log_file, netlist_file);
    [status, ~] = system(cmd);

    if status ~= 0 || ~exist(log_file, 'file')
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e8 + 1e4 * dist; cost = cost(1);
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

    Gain_db = parse_scalar(txt, '(?i)Gain_db\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);
    Power = parse_scalar(txt, '(?i)Power_W\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', NaN);

    has_ac_data = ~isempty(regexp(txt, 'No\. of Data Rows\s*:\s*\d+', 'once'));
    fatal = ~isempty(regexpi(txt, 'fatal|singular matrix|no convergence', 'once'));

    if ~isfinite(Gain_db) || ~isfinite(Power) || Gain_db < -100 || Power < 0
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e10 + 1e6 * dist; cost = cost(1); return;
    end

    gain_target = 15;
    power_target = 20e-3;

    gain_penalty = max(0, (gain_target - Gain_db) / gain_target);
    power_penalty = max(0, (Power - power_target) / power_target);

    w_gain = 0.40;
    w_power = 0.30;
    base_cost = w_gain * (gain_target / max(Gain_db, 0.1)) + ...
                w_power * (Power / power_target);

    if ~has_ac_data || fatal || status ~= 0
        dist = sum(((x_clamped - center_point) ./ (ub - lb)).^2);
        cost = 1e12 + 1e6 * dist; cost = cost(1); return;
    end

    cost = base_cost + 1e4 * (gain_penalty + power_penalty);
    if ~isfinite(cost) || cost > 1e10, cost = 1e10; end
end

function cleanup_temp_files(varargin)
    for i = 1:nargin
        if exist(varargin{i}, 'file')
            delete(varargin{i});
        end
    end
end