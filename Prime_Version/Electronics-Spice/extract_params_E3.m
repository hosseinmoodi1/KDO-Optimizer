function params = extract_params_E3(x)
% EXTRACT_PARAMS_E3 - Extract physical parameters from E3 solution
% OPTIMIZED: Reduced Nmc from 50 to 10 for speed

    params = struct();
    params.f_osc_Hz = NaN;
    params.f_osc_std = NaN;
    params.PN_dBc   = NaN;
    params.PN_std   = NaN;
    params.Power_W  = NaN;
    params.TR_percent = NaN;
    params.FOM_dB   = NaN;

    Vdd = 1.2;
    gamma = 1.5;
    k_B = 1.38e-23;
    T = 300;
    QL = 10;

    L = min(max(abs(x(1)), 1e-9), 20e-9);
    Cvar = min(max(abs(x(2)), 100e-15), 2e-12);
    Cfix = min(max(abs(x(3)), 50e-15), 500e-15);
    Wn = min(max(abs(x(4)), 2e-6), 100e-6);
    Wp = min(max(abs(x(5)), 2e-6), 200e-6);
    Ibias = min(max(abs(x(6)), 1e-3), 10e-3);
    Rbias = min(max(abs(x(7)), 200), 10e3);
    Cbuf = min(max(abs(x(8)), 50e-15), 1e-9);
    Vctrl_min = min(max(x(9), 0.1), 1.2);
    Vctrl_max = min(max(x(10), Vctrl_min + 0.2), 1.2);
    K_vco = min(max(abs(x(11)), 100e6), 10e9);

    % ----- REDUCED Monte Carlo (10 instead of 50) -----
    sigma_L = 0.05;
    sigma_C = 0.03;
    sigma_Vth = 0.02;
    Nmc = 10;  % Reduced from 50 for speed

    f_osc_samples = zeros(Nmc, 1);
    PN_samples = zeros(Nmc, 1);
    Power_samples = zeros(Nmc, 1);

    for mc = 1:Nmc
        L_eff = L * (1 + sigma_L * randn);
        Cvar_eff = Cvar * (1 + sigma_C * randn);
        Cfix_eff = Cfix * (1 + sigma_C * randn);
        Vth_shift = sigma_Vth * randn;

        Ctank = Cvar_eff + Cfix_eff + 0.2 * Cbuf;
        omega0 = 1 / sqrt(max(L_eff * Ctank, 1e-24));
        f_osc = omega0 / (2 * pi);
        f_osc = min(max(f_osc, 0.5e9), 6e9);
        f_osc_samples(mc) = f_osc;

        Id_half = Ibias / 2;
        gm_n = sqrt(2 * 300e-6 * (Wn / 1e-6) * Id_half);
        gm_p = sqrt(2 * 100e-6 * (Wp / 1e-6) * Id_half);
        gm_total = gm_n + gm_p;

        Rtank = QL * omega0 * L_eff;
        V_sig = Id_half * Rtank;
        P_sig = 0.5 * V_sig^2 / max(Rtank, 1);

        F = 1 + 4 * gamma * gm_total * Rtank / (1 + (Vdd / 0.5)^2);
        f_offset = 1e6;
        PN_lin = (F * k_B * T / max(P_sig, 1e-20)) * (f_osc / (2 * QL * f_offset))^2;
        PN_dB = 10 * log10(max(PN_lin, 1e-30));
        PN_samples(mc) = min(max(PN_dB, -150), -80);

        Power_samples(mc) = Vdd * Ibias;
    end

    f_osc_mean = mean(f_osc_samples);
    f_osc_std = std(f_osc_samples);
    PN_mean = mean(PN_samples);
    PN_std = std(PN_samples);
    Power_mean = mean(Power_samples);

    f_min = 1 / (2 * pi * sqrt(L * (Cvar + Cfix + 0.5 * Cbuf)));
    f_max = 1 / (2 * pi * sqrt(L * (Cvar*0.3 + Cfix + 0.5 * Cbuf)));
    TR = (f_max - f_min) / mean([f_max, f_min]) * 100;
    FOM = PN_mean + 20 * log10(f_osc_mean / 1e6) - 10 * log10(Power_mean / 1e-3);

    params.f_osc_Hz = f_osc_mean;
    params.f_osc_std = f_osc_std;
    params.PN_dBc   = PN_mean;
    params.PN_std   = PN_std;
    params.Power_W  = Power_mean;
    params.TR_percent = TR;
    params.FOM_dB   = FOM;
end