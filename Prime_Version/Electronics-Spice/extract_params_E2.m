function params = extract_params_E2(x)
% extract_params_E2: Analytical estimate of physical parameters for E2 (LNA)
% Dimension: 12 variables (same as Electro_E2_Hard.m and template_E2.sp)
% UPDATED: Inductor bounds [0.5nH, 20nH], Id [0.1mA, 5mA].

    params = struct();
    params.Gain_dB = NaN;
    params.Gain_std = NaN;
    params.NF_dB   = NaN;
    params.S11_dB  = NaN;
    params.S22_dB  = NaN;
    params.Power_W = NaN;
    params.Yield   = NaN;

    % Technology parameters (65nm CMOS)
    Vth0 = 0.35;
    mu_n_Cox = 300e-6;
    Cox = 15e-3;
    f0 = 2.4e9;
    omega0 = 2 * pi * f0;
    Z0 = 50;

    % Extract variables with bounds (aligned with Electro_E2_Hard.m)
    Wm1 = min(max(abs(x(1)), 0.5e-6), 50e-6);
    Wm2 = min(max(abs(x(2)), 0.5e-6), 50e-6);
    Wm3 = min(max(abs(x(3)), 0.5e-6), 50e-6);
    Wm4 = min(max(abs(x(4)), 0.5e-6), 50e-6);
    Lg1 = min(max(abs(x(5)), 0.5e-9), 20e-9);   % UPDATED: max 20nH
    Lg2 = min(max(abs(x(6)), 0.5e-9), 20e-9);   % UPDATED
    Ld  = min(max(abs(x(7)), 0.5e-9), 20e-9);   % UPDATED
    Cgd1 = min(max(abs(x(8)), 10e-15), 500e-15);
    Cgd2 = min(max(abs(x(9)), 10e-15), 500e-15);
    Cgd3 = min(max(abs(x(10)), 10e-15), 500e-15);
    Id = min(max(abs(x(11)), 0.1e-3), 5e-3);    % UPDATED: min 0.1mA
    Vdd = min(max(x(12), 0.5), 1.8);

    % Process variation parameters
    sigma_Vth = 0.035;
    sigma_W = 0.03;
    sigma_L = 0.03;
    Nmc = 50;

    Gain_samples = zeros(Nmc, 1);
    NF_samples = zeros(Nmc, 1);
    S11_samples = zeros(Nmc, 1);
    S22_samples = zeros(Nmc, 1);

    for mc = 1:Nmc
        Vth_shift = sigma_Vth * randn;
        Wm1_eff = Wm1 * (1 + sigma_W * randn);
        Wm2_eff = Wm2 * (1 + sigma_W * randn);
        Lg1_eff = Lg1 * (1 + sigma_L * randn);
        Lg2_eff = Lg2 * (1 + sigma_L * randn);
        Ld_eff  = Ld * (1 + sigma_L * randn);

        Vth_eff = Vth0 + Vth_shift;
        Vov = sqrt(2 * Id / (mu_n_Cox * (Wm1_eff / 1e-6)));
        gm = 2 * Id / max(Vov, 0.05);
        ro = 10 / max(Id, 1e-9);

        Cgs = (2/3) * Wm1_eff * 1e-6 * Cox;
        Zin = 1/(1j*omega0*Cgd1) + 1j*omega0*Lg1_eff + 1/(gm + 1j*omega0*Cgs);

        Gain = gm * ro / (1 + gm * ro * (Lg2_eff * omega0 / Z0));
        Gain_dB = 20 * log10(min(max(Gain, 1), 100));
        Gain_samples(mc) = Gain_dB;

        gamma = 2/3;
        NF = 1 + gamma * gm * Z0 / (1 + gm * Z0)^2;
        NF_dB = 10 * log10(max(NF, 1));
        NF_samples(mc) = min(max(NF_dB, 1), 10);

        S11 = 20 * log10(abs((Zin - Z0) / (Zin + Z0)));
        S11_samples(mc) = max(min(S11, 0), -40);

        Zout = 1j * omega0 * Ld_eff + ro;
        S22 = 20 * log10(abs((Zout - Z0) / (Zout + Z0)));
        S22_samples(mc) = max(min(S22, 0), -40);
    end

    % Statistical metrics
    Gain_mean = mean(Gain_samples);
    Gain_std = std(Gain_samples);
    NF_mean = mean(NF_samples);
    S11_mean = mean(S11_samples);
    S22_mean = mean(S22_samples);
    Power = Vdd * Id;
    yield = sum(Gain_samples > 15 & NF_samples < 2.5 & S11_samples < -10) / Nmc;

    params.Gain_dB = Gain_mean;
    params.Gain_std = Gain_std;
    params.NF_dB   = NF_mean;
    params.S11_dB  = S11_mean;
    params.S22_dB  = S22_mean;
    params.Power_W = Power;
    params.Yield   = yield;
end