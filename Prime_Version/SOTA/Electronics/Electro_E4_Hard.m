function cost = Electro_E4_Hard(x)
% Electro_E4_Hard: Robust 6T FinFET SRAM Cell Optimization
% Dimension: 9 variables
%   x(1) = Wfin_pd (pull-down fin width)
%   x(2) = Wfin_pu (pull-up fin width)
%   x(3) = Wfin_ax (access transistor fin width)
%   x(4) = Wfin_ps (precharge fin width)
%   x(5) = Lfin_pd (pull-down fin length)
%   x(6) = Lfin_pu (pull-up fin length)
%   x(7) = Lfin_ax (access fin length)
%   x(8) = VDD     (supply voltage)
%   x(9) = Tox     (oxide thickness)
%
% References:
%   - FinFET 22nm technology parameters
%   - IEEE T-VLSI, "SRAM Design for FinFET Technology"

    % FinFET technology parameters (22nm node)
    mu_n = 270e-4;          % NMOS mobility (m^2/V*s)
    mu_p = 80e-4;           % PMOS mobility (m^2/V*s)
    Vth_n0 = 0.28;          % NMOS threshold (V)
    Vth_p0 = -0.28;         % PMOS threshold (V)
    Vt = 0.026;             % Thermal voltage (V)
    
    % Extract variables with bounds
    Wfin_pd = max(abs(x(1)), 0.5);
    Wfin_pu = max(abs(x(2)), 0.5);
    Wfin_ax = max(abs(x(3)), 0.5);
    Wfin_ps = max(abs(x(4)), 0.5);
    Lfin_pd = max(abs(x(5)), 14e-9);
    Lfin_pu = max(abs(x(6)), 14e-9);
    Lfin_ax = max(abs(x(7)), 14e-9);
    VDD = min(max(x(8), 0.6), 1.0);
    Tox = max(abs(x(9)), 0.8e-9);
    
    % Process variation parameters
    sigma_Vth = 0.025;      % 25mV threshold variation
    sigma_L = 0.02;         % 2% length variation
    sigma_W = 0.02;         % 2% width variation
    Nmc = 200;              % Monte Carlo samples
    
    % Monte Carlo simulation
    SNM_samples = zeros(Nmc, 1);
    I_leak_samples = zeros(Nmc, 1);
    CR_samples = zeros(Nmc, 1);
    PR_samples = zeros(Nmc, 1);
    
    for mc = 1:Nmc
        % Apply process variations
        dVth_n = sigma_Vth * randn;
        dVth_p = sigma_Vth * randn;
        dLpd = Lfin_pd * sigma_L * randn;
        dLpu = Lfin_pu * sigma_L * randn;
        dLax = Lfin_ax * sigma_L * randn;
        dWpd = Wfin_pd * sigma_W * randn;
        dWpu = Wfin_pu * sigma_W * randn;
        dWax = Wfin_ax * sigma_W * randn;
        
        % Effective dimensions
        Lpd_eff = max(Lfin_pd + dLpd, 10e-9);
        Lpu_eff = max(Lfin_pu + dLpu, 10e-9);
        Lax_eff = max(Lfin_ax + dLax, 10e-9);
        Wpd_eff = max(Wfin_pd + dWpd, 0.1);
        Wpu_eff = max(Wfin_pu + dWpu, 0.1);
        Wax_eff = max(Wfin_ax + dWax, 0.1);
        
        % Beta ratios (drive strengths)
        beta_pd = mu_n * (Wpd_eff / Lpd_eff);
        beta_pu = mu_p * (Wpu_eff / Lpu_eff);
        beta_ax = mu_n * (Wax_eff / Lax_eff);
        
        % Cell Ratio (CR) and Pull-up Ratio (PR)
        CR = beta_pd / max(beta_ax, 1e-20);
        PR = beta_pu / max(beta_ax, 1e-20);
        CR_samples(mc) = min(max(CR, 0.5), 3);
        PR_samples(mc) = min(max(PR, 0.2), 1.5);
        
        % Static Noise Margin (Seevinck model)
        SNM = (VDD / 2) * (1 - 1/(1+CR)) * (1 - PR/(1+PR));
        
        % Add Vth variation effect
        Vth_n_eff = Vth_n0 + dVth_n;
        Vth_p_eff = Vth_p0 + dVth_p;
        SNM = SNM * (1 - 0.1 * abs(Vth_n_eff - Vth_n0)/Vth_n0);
        SNM_samples(mc) = max(min(SNM, VDD/2), 0.01);
        
        % Leakage current
        I_leak_n = beta_pd * exp(-max(Vth_n_eff, 0.1) / Vt);
        I_leak_p = beta_pu * exp(-max(abs(Vth_p_eff), 0.1) / Vt);
        I_leak = (I_leak_n + I_leak_p) * VDD * 1e-6;
        I_leak_samples(mc) = min(max(I_leak, 1e-12), 1e-4);
    end
    
    % Statistical metrics
    SNM_mean = mean(SNM_samples);
    SNM_std = std(SNM_samples);
    SNM_3sigma = SNM_mean - 3 * SNM_std;
    I_leak_mean = mean(I_leak_samples);
    I_leak_std = std(I_leak_samples);
    CR_mean = mean(CR_samples);
    PR_mean = mean(PR_samples);
    
    % Yield: SNM >= 150mV
    SNM_req = 0.15;
    yield = sum(SNM_samples >= SNM_req) / Nmc;
    
    % Area estimation
    Area = (Wfin_pd * Lfin_pd + Wfin_pu * Lfin_pu + 2 * Wfin_ax * Lfin_ax) * 1e9;
    
    % ========== CONSTRAINTS ==========
    penalty = 0;
    
    % Constraint 1: Mean SNM >= 150mV
    if SNM_mean < SNM_req
        penalty = penalty + 30 * ((SNM_req - SNM_mean) / 0.05)^2;
    end
    
    % Constraint 2: 3-sigma SNM >= 100mV (robustness)
    if SNM_3sigma < 0.10
        penalty = penalty + 20 * ((0.10 - SNM_3sigma) / 0.05)^2;
    end
    
    % Constraint 3: Yield >= 99%
    if yield < 0.99
        penalty = penalty + 50 * ((0.99 - yield) / 0.01)^2;
    end
    
    % Constraint 4: Cell Ratio >= 1.2
    if CR_mean < 1.2
        penalty = penalty + 10 * ((1.2 - CR_mean) / 0.5)^2;
    end
    
    % Constraint 5: Pull-up Ratio <= 0.8
    if PR_mean > 0.8
        penalty = penalty + 10 * ((PR_mean - 0.8) / 0.5)^2;
    end
    
    % Constraint 6: Leakage current <= 10nA
    if I_leak_mean > 10e-9
        penalty = penalty + 15 * ((I_leak_mean - 10e-9) / 10e-9)^2;
    end
    
    % Constraint 7: Leakage variation <= 30%
    if I_leak_std / I_leak_mean > 0.3
        penalty = penalty + 10 * ((I_leak_std/I_leak_mean - 0.3) / 0.3)^2;
    end
    
    % ========== COST FUNCTION ==========
    % Weights
    w_snm = 0.40;
    w_yield = 0.30;
    w_leak = 0.20;
    w_area = 0.10;
    
    cost = w_snm * (0.15 / max(SNM_mean, 0.01)) + ...
           w_yield * (1 / max(yield, 0.001)) + ...
           w_leak * (I_leak_mean / 10e-9) + ...
           w_area * (Area / 100) + ...
           1e5 * penalty;
    
    if ~isfinite(cost) || cost > 1e10
        cost = 1e10;
    end
end