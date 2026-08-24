%% ========================================================================
%  paths_local.m - Local Path Configuration (Portable for GitHub)
%  ========================================================================
%  AUTHOR: Hossein Moodi
%  DATE:   16-Aug-2026 (FINAL)
%  DESCRIPTION:
%    Defines all project directory paths. Uses persistent variable to
%    ensure the timestamp folder is created only once per MATLAB session.
%  ========================================================================

function paths = paths_local()

    % Persistent variable to store results folder path across calls
    persistent saved_paths;

    % If already computed, return the saved paths
    if ~isempty(saved_paths) && isstruct(saved_paths)
        paths = saved_paths;
        return;
    end

    % Get the directory where this file resides (project root)
    project_root = fileparts(mfilename('fullpath'));
    
    paths.project_root = project_root;
    paths.scripts_dir  = project_root;
    paths.templates_dir = fullfile(project_root, 'Templates');
    paths.models_dir    = fullfile(project_root, 'Models');
    paths.params_dir    = fullfile(project_root, 'Params');
    paths.runs_dir      = fullfile(project_root, 'Outputs', 'Runs');
    paths.results_base_dir = fullfile(project_root, 'Results_Electronics_Hard');
    
    % Create required directories if they don't exist
    if ~exist(paths.templates_dir, 'dir'), mkdir(paths.templates_dir); end
    if ~exist(paths.models_dir, 'dir'), mkdir(paths.models_dir); end
    if ~exist(paths.params_dir, 'dir'), mkdir(paths.params_dir); end
    if ~exist(paths.runs_dir, 'dir'), mkdir(paths.runs_dir); end
    if ~exist(paths.results_base_dir, 'dir'), mkdir(paths.results_base_dir); end
    
    % --------------------------------------------------------------------
    % Timestamp folder (created only once per MATLAB session)
    % --------------------------------------------------------------------
    % Check if there's an existing timestamp folder from this session
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    paths.timestamp_run = timestamp;
    paths.specific_results_dir = fullfile(paths.results_base_dir, timestamp);
    
    % Create the specific results folder if it doesn't exist
    if ~exist(paths.specific_results_dir, 'dir')
        mkdir(paths.specific_results_dir);
        fprintf('[PATHS] Created results folder: %s\n', paths.specific_results_dir);
    else
        fprintf('[PATHS] Using existing results folder: %s\n', paths.specific_results_dir);
    end
    
    % Raw data directory inside specific results
    paths.raw_data_dir = fullfile(paths.specific_results_dir, 'Raw_Data');
    if ~exist(paths.raw_data_dir, 'dir')
        mkdir(paths.raw_data_dir);
    end
    
    % --------------------------------------------------------------------
    % NGSPICE executable path
    % --------------------------------------------------------------------
    if ispc
        paths.ngspice_exe = 'F:\Software\ngspice-46_64\Spice64\bin\ngspice_con.exe';
    else
        paths.ngspice_exe = 'ngspice';
    end

    % Save the computed paths for future calls
    saved_paths = paths;

end