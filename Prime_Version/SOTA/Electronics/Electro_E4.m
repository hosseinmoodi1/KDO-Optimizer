function cost = Electro_E4(x)
% =========================================================
% Problem E4: Robust FinFET 6T-SRAM Design under Process Variation
% =========================================================
%
% Objective (minimization):
%   - Maximize mean Static Noise Margin (SNM) and manufacturing yield
%   - Minimize leakage power and cell area
%   All under process-induced parameter variations (Monte Carlo)
%
% Design Variables (8D):
%   x(1) = Wfin_pd : Fin width/count, pull-down NMOS  [1  .. 10]  (fins)
%   x(2) = Wfin_pu : Fin width/count, pull-up  PMOS   [1  .. 10]  (fins)
%   x(3) = Wfin_ax : Fin width/count, access   NMOS   [1  .. 10]  (fins)
%   x(4) = Lfin_pd : Channel length, pull-down         [10e-9 .. 40e-9] m
%   x(5) = Lfin_pu : Channel length, pull-up           [10e-9 .. 40e-9] m
%   x(6) = Lfin_ax : Channel length, access            [10e-9 .. 40e-9] m
%   x(7) = VDD     : Supply voltage                    [0.6   .. 1.0]   V
%   x(8) = Tox     : Equivalent gate oxide thickness   [0.5e-9.. 1.5e-9] m
%
% Soft Constraints:
%   Mean SNM  >= 0.15 V
%   Yield (fraction of MC samples with SNM >= 0.15 V) >= 95%
%   Max leakage current <= 1e-5 A
%
% Suggested bounds for optimizer:
%   lb = [1,  1,  1,  10e-9, 10e-9, 10e-9, 0.6, 0.5e-9]
%   ub = [10, 10, 10, 40e-9, 40e-9, 40e-9, 1.0, 1.5e-9]
% =========================================================

    % -------------------------------------------------------
    % Nominal FinFET Technology Parameters
    % -------------------------------------------------------
    mu_n   = 200e-4;    % Effective electron mobility (m^2/V/s)
    mu_p   = 80e-4;     % Effective hole mobility (m^2/V/s)
    Vth_n0 = 0.25;      % Nominal NMOS threshold voltage (V)
    Vth_p0 = -0.25;     % Nominal PMOS threshold voltage (V)
    Vt     = 0.026;     % Thermal voltage kT/q at 300 K (V)

    % -------------------------------------------------------
    % Process Variation Parameters (1-sigma Gaussian)
    % -------------------------------------------------------
    sigma_Vth = 20e-3;  % Threshold voltage mismatch (V), typical FinFET
    sigma_L   = 1e-9;   % Channel length variation (m), LER contribution
    sigma_W   = 0.2;    % Fin count variation (fins), discrete granularity

    % -------------------------------------------------------
    % Monte Carlo Setup
    % -------------------------------------------------------
    Nmc = 50;           % Number of MC samples (reduced for runtime speed)

    % -------------------------------------------------------
    % Extract & Clip Design Variables
    % -------------------------------------------------------
    Wfin_pd = max(abs(x(1)), 1);        % pull-down fins (dimensionless)
    Wfin_pu = max(abs(x(2)), 1);        % pull-up fins
    Wfin_ax = max(abs(x(3)), 1);        % access fins
    Lfin_pd = max(abs(x(4)), 10e-9);    % pull-down channel length (m)
    Lfin_pu = max(abs(x(5)), 10e-9);    % pull-up channel length (m)
    Lfin_ax = max(abs(x(6)), 10e-9);    % access channel length (m)
    VDD     = min(max(x(7), 0.6), 1.0); % supply voltage (V), hard-clipped
    Tox     = max(abs(x(8)), 0.5e-9);   % gate oxide thickness (m)

    % -------------------------------------------------------
    % Oxide Capacitance Effect on Threshold (Tox influence)
    % -------------------------------------------------------
    % Thinner Tox ? stronger gate control ? lower effective Vth
    % Simple linear model: delta_Vth ~ -0.05 * (Tox_ref - Tox) / Tox_ref
    Tox_ref   = 1e-9;                   % reference oxide thickness (m)
    dVth_tox  = -0.05 * (Tox_ref - Tox) / Tox_ref;  % Tox-induced Vth shift (V)
    Vth_n_nom = Vth_n0 + dVth_tox;     % effective nominal NMOS Vth (V)
    Vth_p_nom = Vth_p0 - dVth_tox;     % effective nominal PMOS Vth (V)

    % -------------------------------------------------------
    % Monte Carlo: SNM and Leakage under Process Variation
    % -------------------------------------------------------
    SNM_samples    = zeros(Nmc, 1);
    I_leak_samples = zeros(Nmc, 1);

    for mc = 1:Nmc
        % --- Random process perturbations (independent per device) ---
        dVth_n = sigma_Vth * randn;     % NMOS threshold variation (V)
        dVth_p = sigma_Vth * randn;     % PMOS threshold variation (V)
        dLpd   = sigma_L   * randn;     % pull-down length variation (m)
        dLpu   = sigma_L   * randn;     % pull-up length variation (m)
        dLax   = sigma_L   * randn;     % access length variation (m)
        dWpd   = sigma_W   * randn;     % pull-down fin variation
        dWpu   = sigma_W   * randn;     % pull-up fin variation
        dWax   = sigma_W   * randn;     % access fin variation

        % --- Perturbed device strengths (beta = mu * W/L) ---
        % Clamp to physical minimums to avoid negative/zero values
        beta_pd = mu_n * max(Wfin_pd + dWpd, 0.1) / max(Lfin_pd + dLpd, 1e-9);
        beta_pu = mu_p * max(Wfin_pu + dWpu, 0.1) / max(Lfin_pu + dLpu, 1e-9);
        beta_ax = mu_n * max(Wfin_ax + dWax, 0.1) / max(Lfin_ax + dLax, 1e-9);

        % --- Perturbed threshold voltages ---
        Vth_n_var = Vth_n_nom + dVth_n;    % varied NMOS Vth (V)
        Vth_p_var = Vth_p_nom + dVth_p;    % varied PMOS Vth (V)

        % --- Static Noise Margin (SNM) approximation ---
        % Cell ratio CR = beta_pd / beta_ax (read stability)
        % Pull-up ratio PR = beta_pu / beta_ax (write-ability)
        CR = beta_pd / max(beta_ax, 1e-20);
        PR = beta_pu / max(beta_ax, 1e-20);

        % Simplified SNM model based on CR and PR
        % Higher CR ? better read SNM; lower PR ? better write margin
        SNM_k = (VDD / 2) * (1 - 1/(1 + CR)) * (1 - PR/(1 + PR));
        SNM_k = max(SNM_k, 1e-9);          % clamp to small positive value
        SNM_samples(mc) = SNM_k;

        % --- Leakage current approximation (subthreshold) ---
        % I_leak ~ beta * exp(-Vth / Vt) for each device
        % PMOS leakage uses |Vth_p_var| (Vth_p_var is negative)
        I_leak_n = (beta_pd + beta_ax) * exp(-max(Vth_n_var, 0.05) / Vt) * 1e-7;
        I_leak_p = beta_pu             * exp(-max(abs(Vth_p_var), 0.05) / Vt) * 1e-7;
        I_leak_samples(mc) = I_leak_n + I_leak_p;
    end

    % -------------------------------------------------------
    % Statistical Metrics
    % -------------------------------------------------------
    SNM_mean    = mean(SNM_samples);        % mean SNM across MC (V)
    SNM_min     = min(SNM_samples);         % worst-case SNM (V)
    I_leak_mean = mean(I_leak_samples);     % mean leakage (A)
    I_leak_max  = max(I_leak_samples);      % worst-case leakage (A)

    % Yield: fraction of MC samples meeting SNM requirement
    SNM_req = 0.15;                         % minimum acceptable SNM (V)
    yield   = sum(SNM_samples >= SNM_req) / Nmc;  % yield [0, 1]

    % -------------------------------------------------------
    % Cell Area Approximation
    % -------------------------------------------------------
    % Area ~ sum of (fin_count * channel_length) for all 3 device types
    % Factor of 2 accounts for complementary pairs (NMOS + PMOS per node)
    Area = 2 * (Wfin_pd*Lfin_pd + Wfin_pu*Lfin_pu + Wfin_ax*Lfin_ax); % (m)

    % -------------------------------------------------------
    % Soft Constraint Penalties
    % -------------------------------------------------------
    penalty = 0;

    % Constraint 1: Mean SNM >= 0.15 V
    if SNM_mean < SNM_req
        penalty = penalty + ((SNM_req - SNM_mean) / SNM_req)^2;
    end

    % Constraint 2: Worst-case SNM >= 0.10 V (additional robustness)
    SNM_min_req = 0.10;
    if SNM_min < SNM_min_req
        penalty = penalty + ((SNM_min_req - SNM_min) / SNM_min_req)^2;
    end

    % Constraint 3: Yield >= 95%
    yield_req = 0.95;
    if yield < yield_req
        penalty = penalty + ((yield_req - yield) / yield_req)^2;
    end

    % Constraint 4: Worst-case leakage <= 1e-5 A
    I_leak_max_req = 1e-5;
    if I_leak_max > I_leak_max_req
        penalty = penalty + ((I_leak_max - I_leak_max_req) / I_leak_max_req)^2;
    end

    % -------------------------------------------------------
    % Multi-Objective Cost (weighted sum, minimization)
    % -------------------------------------------------------
    % Weights: SNM (0.4), yield (0.3), leakage (0.2), area (0.1)
    w_snm   = 0.4;
    w_yield = 0.3;
    w_leak  = 0.2;
    w_area  = 0.1;

    % Area normalization: reference = 10 fins * 40 nm length
    Area_ref = 10 * 40e-9;

    cost = w_snm   * (1 / max(SNM_mean, 1e-9))          + ...  % maximize SNM
           w_yield * (1 / max(yield, 1e-3))              + ...  % maximize yield
           w_leak  * (I_leak_mean / I_leak_max_req)      + ...  % minimize leakage
           w_area  * (Area / Area_ref)                   + ...  % minimize area
           1e4     * penalty;                                    % constraint violations

    % Guard: return large finite value instead of NaN/Inf
    if ~isfinite(cost)
        cost = 1e10;
    end

end
