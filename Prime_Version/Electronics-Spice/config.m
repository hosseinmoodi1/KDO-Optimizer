function cfg = config()
% config.m - Global Configuration File
% UPDATED: Uses persistent paths from paths_local().

    paths = paths_local();

    cfg.simulator   = 'ngspice';
    cfg.ngspice_exe = paths.ngspice_exe;
    
    cfg.project_root     = paths.project_root;
    cfg.scripts_dir      = paths.scripts_dir;
    cfg.models_dir       = paths.models_dir;
    cfg.templates_dir    = paths.templates_dir;
    cfg.params_dir       = paths.params_dir;
    cfg.runs_dir         = paths.runs_dir;
    cfg.results_dir      = paths.specific_results_dir;
    cfg.raw_data_dir     = paths.raw_data_dir;

    cfg.results_mat_file = fullfile(cfg.raw_data_dir, 'Electronics_Results_Hard.mat');
    cfg.algorithms_to_validate = {'KDO', 'L_SHADE', 'DE', 'GBO'};

    cfg.problem_names = {'E1', 'E2', 'E3', 'E4'};

    % SPICE metrics
    cfg.spice_metrics.E1 = {'Gain_db', 'GBW', 'PM', 'Power_W'};
    cfg.spice_metrics.E2 = {'Gain_db', 'Power_W'};
    cfg.spice_metrics.E3 = {'f_osc_Hz', 'Power_W'};
    cfg.spice_metrics.E4 = {'SNM_V', 'I_leak_A'};

    % MATLAB metrics
    cfg.matlab_metrics.E1 = {'Gain_db', 'GBW', 'PM', 'Power_W'};
    cfg.matlab_metrics.E2 = {'Gain_db', 'Power_W'};
    cfg.matlab_metrics.E3 = {'f_osc_Hz', 'Power_W'};
    cfg.matlab_metrics.E4 = {'SNM_V', 'I_leak_A', 'CR', 'PR'};

    % Unit scaling (all 1 = SI units)
    cfg.unit_scaling.E1.Gain_db = 1;
    cfg.unit_scaling.E1.GBW     = 1;
    cfg.unit_scaling.E1.PM      = 1;
    cfg.unit_scaling.E1.Power_W = 1;

    cfg.unit_scaling.E2.Gain_db = 1;
    cfg.unit_scaling.E2.Power_W = 1;

    cfg.unit_scaling.E3.f_osc_Hz = 1;
    cfg.unit_scaling.E3.Power_W  = 1;

    cfg.unit_scaling.E4.SNM_V    = 1;
    cfg.unit_scaling.E4.I_leak_A = 1;
    cfg.unit_scaling.E4.CR       = 1;
    cfg.unit_scaling.E4.PR       = 1;

    cfg.epsilon = 1e-12;
    cfg.results_dir = paths.specific_results_dir;
    cfg.runs_dir = paths.runs_dir;
    cfg.debug_mode = false;

end