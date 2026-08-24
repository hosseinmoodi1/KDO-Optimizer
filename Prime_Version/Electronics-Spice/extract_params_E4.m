function params = extract_params_E4(x)
% extract_params_E4: Extract physical parameters from E4 (SRAM) solution
% This function is used for reporting only, not during optimization.
%
% SNM is computed analytically from CR, PR, and VDD.

    params = struct();
    params.SNM_mean_V = NaN;
    params.SNM_std_V = NaN;
    params.SNM_3sigma_V = NaN;
    params.Yield = NaN;
    params.I_leak_A = NaN;
    params.I_leak_std = NaN;
    params.CR = NaN;
    params.PR = NaN;
    params.Area_um2 = NaN;

    % Physical constants
    mu_n = 270e-4; mu_p = 80e-4;
    Vth_n0 = 0.28; Vth_p0 = -0.28;
    Vt = 0.026;
    n = 1.5;
    epsilon0 = 8.854e-12;
    epsilon_ox = 3.9;
    W_single_fin = 40e-9;

    % Clamp variables (same as cost function)
    Wfin_pd = round(min(max(abs(x(1)), 1), 8));
    Wfin_pu = round(min(max(abs(x(2)), 1), 8));
    Wfin_ax = round(min(max(abs(x(3)), 1), 8));
    Wfin_ps = round(min(max(abs(x(4)), 1), 8));
    Lfin_pd = min(max(abs(x(5)), 20e-9), 40e-9);
    Lfin_pu = min(max(abs(x(6)), 20e-9), 40e-9);
    Lfin_ax = min(max(abs(x(7)), 20e-9), 40e-9);
    VDD = min(max(x(8), 0.6), 1.0);
    Tox = max(abs(x(9)), 0.8e-9);

    % Compute Cox from Tox
    Cox = epsilon_ox * epsilon0 / Tox;

    sigma_Vth = 0.025; sigma_L = 0.02; sigma_W = 0.02; Nmc = 200;
    SNM_samples = []; I_leak_samples = []; CR_samples = []; PR_samples = [];

    for mc = 1:Nmc
        dVth_n = sigma_Vth * randn; dVth_p = sigma_Vth * randn;
        Lpd_eff = max(Lfin_pd + Lfin_pd*sigma_L*randn, 10e-9);
        Lpu_eff = max(Lfin_pu + Lfin_pu*sigma_L*randn, 10e-9);
        Lax_eff = max(Lfin_ax + Lfin_ax*sigma_L*randn, 10e-9);
        Wpd_eff = max(Wfin_pd + Wfin_pd*sigma_W*randn, 0.5) * W_single_fin;
        Wpu_eff = max(Wfin_pu + Wfin_pu*sigma_W*randn, 0.5) * W_single_fin;
        Wax_eff = max(Wfin_ax + Wfin_ax*sigma_W*randn, 0.5) * W_single_fin;

        mu_n_Cox = mu_n * Cox; mu_p_Cox = mu_p * Cox;
        beta_pd = mu_n_Cox * (Wpd_eff / Lpd_eff);
        beta_pu = mu_p_Cox * (Wpu_eff / Lpu_eff);
        beta_ax = mu_n_Cox * (Wax_eff / Lax_eff);

        CR = beta_pd / max(beta_ax, 1e-20);
        PR = beta_pu / max(beta_ax, 1e-20);

        % Analytical SNM
        SNM = (VDD / 2) * (1 - 1/(1+CR)) * (1 - PR/(1+PR));
        Vth_n_eff = Vth_n0 + dVth_n; Vth_p_eff = Vth_p0 + dVth_p;
        SNM = SNM * (1 - 0.5 * abs(Vth_n_eff - Vth_n0) / Vth_n0);
        SNM = max(min(SNM, VDD/2), 0);

        % Leakage current
        I_leak_n = beta_pd * Vt^2 * exp(-max(Vth_n_eff, 0.1) / (n * Vt)) * (1 - exp(-VDD/Vt));
        I_leak_p = beta_pu * Vt^2 * exp(-max(abs(Vth_p_eff), 0.1) / (n * Vt)) * (1 - exp(-VDD/Vt));
        I_leak = I_leak_n + I_leak_p;

        SNM_samples = [SNM_samples; SNM];
        I_leak_samples = [I_leak_samples; I_leak];
        CR_samples = [CR_samples; CR];
        PR_samples = [PR_samples; PR];
    end

    SNM_mean = mean(SNM_samples); SNM_std = std(SNM_samples);
    SNM_3sigma = SNM_mean - 3*SNM_std;
    I_leak_mean = mean(I_leak_samples); I_leak_std = std(I_leak_samples);
    CR_mean = mean(CR_samples); PR_mean = mean(PR_samples);
    yield = sum(SNM_samples >= 0.15) / Nmc;
    Area_m2 = W_single_fin * (Wfin_pd * Lfin_pd + Wfin_pu * Lfin_pu + 2 * Wfin_ax * Lfin_ax);
    Area = Area_m2 * 1e12;

    params.SNM_mean_V = SNM_mean;
    params.SNM_std_V = SNM_std;
    params.SNM_3sigma_V = SNM_3sigma;
    params.Yield = yield;
    params.I_leak_A = I_leak_mean;
    params.I_leak_std = I_leak_std;
    params.CR = CR_mean;
    params.PR = PR_mean;
    params.Area_um2 = Area;
end