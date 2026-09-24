function cost = Electro_E2_Hard(x)
% Electro_E2_Hard: 2.4 GHz CMOS Cascode LNA with Process Variation
% Dimension: 12 variables
%   x(1)  = Wm1   (input transistor width)
%   x(2)  = Wm2   (cascode transistor width)
%   x(3)  = Wm3   (load transistor width)
%   x(4)  = Wm4   (current source width)
%   x(5)  = Lg1   (gate inductance)
%   x(6)  = Lg2   (source inductance)
%   x(7)  = Ld    (load inductance)
%   x(8)  = Cgd1  (gate-drain cap M1)
%   x(9)  = Cgd2  (gate-drain cap M2)
%   x(10) = Cgd3  (gate-drain cap M3)
%   x(11) = Id    (bias current)
%   x(12) = Vdd   (supply voltage)
%
% References:
%   - Razavi, "RF Microelectronics", 2nd Edition
%   - CMOS 65nm technology parameters

    % Technology parameters (65nm CMOS)
    Vth0 = 0.35;            % Threshold voltage (V)
    mu_n_Cox = 300e-6;      % NMOS transconductance (A/V^2)
    Cox = 15e-3;            % Gate oxide capacitance (F/m^2)
    f0 = 2.4e9;             % Operating frequency (Hz)
    omega0 = 2 * pi * f0;
    Z0 = 50;                % Characteristic impedance (ohms)
    
    % Extract variables with bounds
    Wm1 = max(abs(x(1)), 0.5e-6);
    Wm2 = max(abs(x(2)), 0.5e-6);
    Wm3 = max(abs(x(3)), 0.5e-6);
    Wm4 = max(abs(x(4)), 0.5e-6);
    Lg1 = max(abs(x(5)), 0.1e-9);
    Lg2 = max(abs(x(6)), 0.1e-9);
    Ld = max(abs(x(7)), 0.1e-9);
    Cgd1 = max(abs(x(8)), 10e-15);
    Cgd2 = max(abs(x(9)), 10e-15);
    Cgd3 = max(abs(x(10)), 10e-15);
    Id = max(abs(x(11)), 0.1e-3);
    Vdd = min(max(x(12), 0.5), 1.8);
    
    % Process variation parameters
    sigma_Vth = 0.035;      % Threshold voltage variation
    sigma_W = 0.03;         % Width variation
    sigma_L = 0.03;         % Length variation
    Nmc = 50;               % Monte Carlo samples
    
    % Monte Carlo simulation
    Gain_samples = zeros(Nmc, 1);
    NF_samples = zeros(Nmc, 1);
    S11_samples = zeros(Nmc, 1);
    S22_samples = zeros(Nmc, 1);
    
    for mc = 1:Nmc
        % Apply process variations
        Vth_shift = sigma_Vth * randn;
        Wm1_eff = Wm1 * (1 + sigma_W * randn);
        Wm2_eff = Wm2 * (1 + sigma_W * randn);
        
        % Effective threshold
        Vth_eff = Vth0 + Vth_shift;
        
        % Overdrive voltage
        Vov = sqrt(2 * Id / (mu_n_Cox * (Wm1_eff / 1e-6)));
        
        % Transconductance
        gm = 2 * Id / max(Vov, 0.05);
        
        % Output resistance
        ro = 10 / max(Id, 1e-9);
        
        % Gain calculation
        Gain = gm * ro / (1 + gm * ro * (Lg2 * omega0 / Z0));
        Gain = min(max(Gain, 1), 100);
        Gain_samples(mc) = 20 * log10(Gain);
        
        % Noise Figure (Frias model)
        gamma = 2 / 3;
        NF = 1 + gamma * gm * Z0 / (1 + gm * Z0)^2;
        NF_dB = 10 * log10(NF);
        NF_samples(mc) = min(max(NF_dB, 1), 10);
        
        % Input matching S11
        Zin = 1/(j*omega0*Cgd1) + j*omega0*Lg1 + 1/(gm);
        S11 = 20 * log10(abs((Zin - Z0) / (Zin + Z0)));
        S11_samples(mc) = max(min(S11, 0), -40);
        
        % Output matching S22
        Zout = j*omega0*Ld + ro;
        S22 = 20 * log10(abs((Zout - Z0) / (Zout + Z0)));
        S22_samples(mc) = max(min(S22, 0), -40);
    end
    
    % Statistical metrics
    Gain_mean = mean(Gain_samples);
    Gain_std = std(Gain_samples);
    NF_mean = mean(NF_samples);
    S11_mean = mean(S11_samples);
    S22_mean = mean(S22_samples);
    
    % Yield calculation
    yield = sum(Gain_samples > 15 & NF_samples < 2.5 & S11_samples < -10) / Nmc;
    
    % Power consumption
    Power = Vdd * Id;
    
    % ========== CONSTRAINTS ==========
    penalty = 0;
    
    % Constraint 1: Gain >= 15 dB
    if Gain_mean < 15
        penalty = penalty + 20 * ((15 - Gain_mean) / 10)^2;
    end
    
    % Constraint 2: Gain variation <= 3 dB
    if Gain_std > 3
        penalty = penalty + 10 * ((Gain_std - 3) / 3)^2;
    end
    
    % Constraint 3: Noise Figure <= 2.5 dB
    if NF_mean > 2.5
        penalty = penalty + 15 * ((NF_mean - 2.5) / 1)^2;
    end
    
    % Constraint 4: Input matching S11 <= -10 dB
    if S11_mean > -10
        penalty = penalty + 10 * ((S11_mean + 10) / 10)^2;
    end
    
    % Constraint 5: Output matching S22 <= -10 dB
    if S22_mean > -10
        penalty = penalty + 5 * ((S22_mean + 10) / 10)^2;
    end
    
    % Constraint 6: Yield >= 99%
    if yield < 0.99
        penalty = penalty + 50 * ((0.99 - yield) / 0.01)^2;
    end
    
    % Constraint 7: Power <= 20 mW
    if Power > 20e-3
        penalty = penalty + 10 * ((Power - 20e-3) / 10e-3)^2;
    end
    
    % ========== COST FUNCTION ==========
    % Weights
    w_gain = 0.30;
    w_nf = 0.30;
    w_s11 = 0.20;
    w_power = 0.20;
    
    cost = w_gain * (20 / max(Gain_mean, 1)) + ...
           w_nf * (NF_mean / 2) + ...
           w_s11 * (abs(S11_mean) / 20) + ...
           w_power * (Power / 20e-3) + ...
           1e4 * penalty;
    
    if ~isfinite(cost) || cost > 1e10
        cost = 1e10;
    end
end