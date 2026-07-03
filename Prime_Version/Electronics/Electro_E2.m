function cost = Electro_E2(x)
% =========================================================
% E2: 6T CMOS SRAM Cell Optimization
% =========================================================
% Objective : Maximize Static Noise Margin (SNM),
%             minimize leakage power
%             (formulated as minimization problem)
%
% Design Variables (6D):
%   x(1) = Wpd : Pull-Down  NMOS width  [100 nm .. 2 µm]
%   x(2) = Wpu : Pull-Up    PMOS width  [100 nm .. 2 µm]
%   x(3) = Wax : Access     NMOS width  [100 nm .. 2 µm]
%   x(4) = Lpd : Pull-Down  length      [45 nm  .. 250 nm]
%   x(5) = Lpu : Pull-Up    length      [45 nm  .. 250 nm]
%   x(6) = Lax : Access     length      [45 nm  .. 250 nm]
%
% Stability Constraints:
%   CR = beta_pd / beta_ax >= 1  (read stability)
%   PR = beta_pu / beta_ax <= 1  (write ability)
%   SNM >= 150 mV
%
% Technology: 45 nm CMOS, VDD = 1.0 V
%
% Reference: Seevinck et al., IEEE JSSC, 1987
% =========================================================

    % --- Technology parameters (45 nm node) ---
    VDD  = 1.0;      % Supply voltage (V)
    mu_n = 270e-4;   % NMOS carrier mobility (m^2/V·s)
    mu_p = 100e-4;   % PMOS carrier mobility (m^2/V·s)
    Vth  = 0.35;     % Threshold voltage (V)

    % --- Extract and sanitize design variables ---
    Wpd = max(abs(x(1)), 1e-10);   % pull-down width
    Wpu = max(abs(x(2)), 1e-10);   % pull-up width
    Wax = max(abs(x(3)), 1e-10);   % access transistor width
    Lpd = max(abs(x(4)), 1e-11);   % pull-down length
    Lpu = max(abs(x(5)), 1e-11);   % pull-up length
    Lax = max(abs(x(6)), 1e-11);   % access transistor length

    % --- Drive strengths (proportional to mu*W/L, Cox absorbed) ---
    beta_pd = mu_n * (Wpd / Lpd);   % pull-down strength
    beta_pu = mu_p * (Wpu / Lpu);   % pull-up strength
    beta_ax = mu_n * (Wax / Lax);   % access transistor strength

    % --- Cell Ratio (CR) and Pull-up Ratio (PR) ---
    CR = beta_pd / max(beta_ax, 1e-20);   % CR > 1 ensures read stability
    PR = beta_pu / max(beta_ax, 1e-20);   % PR < 1 ensures write ability

    % --- Static Noise Margin (analytical approximation, Seevinck model) ---
    SNM = (VDD/2) * (1 - 1/(1+CR)) * (1 - PR/(1+PR));
    SNM = max(SNM, 1e-9);   % clip to avoid 1/SNM singularity

    % --- Subthreshold leakage power ---
    % I_leak proportional to beta * exp(-Vth/Vt), Vt = 26 mV at 300 K
    I_leak = (beta_pd + beta_pu + beta_ax) * exp(-Vth/0.026) * 1e-6;
    Power  = VDD * I_leak;

    % --- Soft penalty for constraint violations ---
    penalty = 0;
    if beta_ax > beta_pd                          % read stability
        penalty = penalty + ((beta_ax - beta_pd) / beta_pd)^2;
    end
    if beta_pu > beta_ax                          % write ability
        penalty = penalty + ((beta_pu - beta_ax) / beta_ax)^2;
    end
    if SNM < 0.15                                 % SNM >= 150 mV
        penalty = penalty + ((0.15 - SNM) / 0.15)^2;
    end

    % --- Cost function (minimization) ---
    % w1 penalizes low SNM, w2 penalizes high leakage power
    w1 = 0.6;   % SNM importance weight
    w2 = 0.4;   % power importance weight
    cost = w1*(1/SNM) + w2*Power*1e9 + 1e4*penalty;

end
