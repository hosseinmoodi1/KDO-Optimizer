function cost = Electro_E3_Hard(x)
% Electro_E3_Hard: Low Phase-Noise LC-VCO with PVT Variation
% Dimension: 11 variables
%   x(1)  = L      (inductor value)
%   x(2)  = Cvar   (varactor capacitance)
%   x(3)  = Cfix   (fixed capacitance)
%   x(4)  = Wn     (NMOS width)
%   x(5)  = Wp     (PMOS width)
%   x(6)  = Ibias  (bias current)
%   x(7)  = Rbias  (bias resistor)
%   x(8)  = Cbuf   (buffer capacitance)
%   x(9)  = Vctrl_min (minimum control voltage)
%   x(10) = Vctrl_max (maximum control voltage)
%   x(11) = K_vco  (VCO gain)
%
% References:
%   - Hajimiri & Lee, "Phase Noise in LC Oscillators", JSSC 1999
%   - CMOS 65nm technology

    % Technology parameters
    Vdd = 1.2;              % Supply voltage (V)
    gamma = 1.5;            % MOSFET excess noise factor
    k_B = 1.38e-23;         % Boltzmann constant
    T = 300;                % Temperature (K)
    QL = 10;                % Inductor quality factor
    
    % Extract variables with bounds
    L = max(abs(x(1)), 0.5e-9);
    Cvar = max(abs(x(2)), 50e-15);
    Cfix = max(abs(x(3)), 10e-15);
    Wn = max(abs(x(4)), 1e-6);
    Wp = max(abs(x(5)), 1e-6);
    Ibias = max(abs(x(6)), 1e-6);
    Rbias = max(abs(x(7)), 100);
    Cbuf = max(abs(x(8)), 10e-15);
    Vctrl_min = min(max(x(9), 0), 1.2);
    Vctrl_max = min(max(x(10), Vctrl_min + 0.1), 1.2);
    K_vco = max(abs(x(11)), 50e6);
    
    % Limit to realistic values
    L = min(L, 20e-9);
    Cvar = min(Cvar, 2e-12);
    Cfix = min(Cfix, 500e-15);
    Wn = min(Wn, 100e-6);
    Wp = min(Wp, 200e-6);
    Ibias = min(Ibias, 10e-3);
    
    % Process variation parameters
    sigma_L = 0.05;         % Inductor variation
    sigma_C = 0.03;         % Capacitor variation
    sigma_Vth = 0.02;       % Threshold variation
    Nmc = 50;
    
    % Monte Carlo simulation
    f_osc_samples = zeros(Nmc, 1);
    PN_samples = zeros(Nmc, 1);
    Power_samples = zeros(Nmc, 1);
    
    for mc = 1:Nmc
        % Apply variations
        L_eff = L * (1 + sigma_L * randn);
        Cvar_eff = Cvar * (1 + sigma_C * randn);
        Cfix_eff = Cfix * (1 + sigma_C * randn);
        Vth_shift = sigma_Vth * randn;
        
        % Total capacitance
        Ctank = Cvar_eff + Cfix_eff + 0.2 * Cbuf;
        
        % Oscillation frequency
        omega0 = 1 / sqrt(max(L_eff * Ctank, 1e-24));
        f_osc = omega0 / (2 * pi);
        f_osc_samples(mc) = min(max(f_osc, 0.5e9), 6e9);
        
        % Effective transconductance
        Id_half = Ibias / 2;
        gm_n = sqrt(2 * 300e-6 * (Wn / 1e-6) * Id_half);
        gm_p = sqrt(2 * 100e-6 * (Wp / 1e-6) * Id_half);
        gm_total = gm_n + gm_p;
        
        % Tank resistance
        Rtank = QL * omega0 * L_eff;
        
        % Signal power
        P_sig = 0.5 * (Id_half * Rtank)^2 / max(Rtank, 1);
        
        % Phase noise (Leeson's model with Hajimiri correction)
        f_offset = 1e6;
        F = (gamma * 4 * k_B * T) / max(P_sig, 1e-20);
        PN_lin = F * (f_osc / (2 * QL * f_offset))^2;
        PN_dB = 10 * log10(max(PN_lin, 1e-30));
        PN_samples(mc) = min(max(PN_dB, -150), -80);
        
        % Power consumption
        Power_samples(mc) = Vdd * Ibias;
    end
    
    % Statistical metrics
    f_osc_mean = mean(f_osc_samples);
    f_osc_std = std(f_osc_samples);
    PN_mean = mean(PN_samples);
    PN_std = std(PN_samples);
    Power_mean = mean(Power_samples);
    
    % Tuning range
    f_min = 1 / (2 * pi * sqrt(L * (Cvar + Cfix + 0.5 * Cbuf)));
    f_max = 1 / (2 * pi * sqrt(L * (Cvar*0.3 + Cfix + 0.5 * Cbuf)));
    TR = (f_max - f_min) / mean([f_max, f_min]) * 100;
    
    % Figure of Merit (FOM)
    FOM = PN_mean + 20 * log10(f_osc_mean / 1e6) - 10 * log10(Power_mean / 1e-3);
    
    % ========== CONSTRAINTS ==========
    penalty = 0;
    
    % Constraint 1: Oscillation frequency in [2.4, 2.6] GHz
    if f_osc_mean < 2.4e9
        penalty = penalty + 20 * ((2.4e9 - f_osc_mean) / 1e9)^2;
    elseif f_osc_mean > 2.6e9
        penalty = penalty + 20 * ((f_osc_mean - 2.6e9) / 1e9)^2;
    end
    
    % Constraint 2: Frequency variation <= 5%
    if f_osc_std / f_osc_mean > 0.05
        penalty = penalty + 10 * ((f_osc_std/f_osc_mean - 0.05) / 0.05)^2;
    end
    
    % Constraint 3: Phase Noise <= -115 dBc/Hz @ 1MHz
    if PN_mean > -115
        penalty = penalty + 30 * ((PN_mean + 115) / 15)^2;
    end
    
    % Constraint 4: PN variation <= 2 dB
    if PN_std > 2
        penalty = penalty + 10 * ((PN_std - 2) / 2)^2;
    end
    
    % Constraint 5: Power <= 8 mW
    if Power_mean > 8e-3
        penalty = penalty + 15 * ((Power_mean - 8e-3) / 4e-3)^2;
    end
    
    % Constraint 6: Tuning range >= 10%
    if TR < 10
        penalty = penalty + 10 * ((10 - TR) / 10)^2;
    end
    
    % Constraint 7: FOM <= -180 dBc/Hz
    if FOM > -180
        penalty = penalty + 5 * ((FOM + 180) / 20)^2;
    end
    
    % ========== COST FUNCTION ==========
    % Weights
    w_f = 0.20;
    w_pn = 0.40;
    w_power = 0.25;
    w_fom = 0.15;
    
    cost = w_f * (abs(f_osc_mean - 2.5e9) / 1e9) + ...
           w_pn * ((-115 - PN_mean) / 20) + ...
           w_power * (Power_mean / 8e-3) + ...
           w_fom * ((-180 - FOM) / 20) + ...
           1e5 * penalty;
    
    if ~isfinite(cost) || cost > 1e10
        cost = 1e10;
    end
end