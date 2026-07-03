function cost = Electro_E1_Hard(x)
% Electro_E1_Hard: Two-Stage Folded-Cascode Operational Amplifier
% Dimension: 18 variables
%   x(1:8)   = W1..W8 (transistor widths)
%   x(9:16)  = L1..L8 (channel lengths)
%   x(17)    = Ibias (bias current)
%   x(18)    = Cc (compensation capacitor)
%
% References:
%   - Mahdavi et al., "Metaheuristic optimization of analog IC design", 2020
%   - CMOS technology parameters: 180nm process

    % Technology parameters
    VDD = 1.8;          % Supply voltage (V)
    VTHN = 0.45;        % NMOS threshold voltage (V)
    VTHP = -0.45;       % PMOS threshold voltage (V)
    UN_Cox = 180e-6;    % NMOS transconductance parameter (A/V^2)
    UP_Cox = 60e-6;     % PMOS transconductance parameter (A/V^2)
    lambda_n = 0.05;    % Channel length modulation (1/V)
    lambda_p = 0.05;
    
    % Extract variables with bounds
    W = max(abs(x(1:8)), 1e-6);
    L = max(abs(x(9:16)), 45e-9);
    Ibias = max(abs(x(17)), 1e-6);
    Cc = max(abs(x(18)), 0.1e-12);
    
    % Limit Cc to reasonable range
    Cc = min(Cc, 10e-12);
    
    % Current mirrors
    I_tail = Ibias;
    I_d = I_tail / 2;
    
    % Transconductance calculations
    gm1 = sqrt(2 * UN_Cox * (W(1)/L(1)) * I_d);
    gm2 = sqrt(2 * UN_Cox * (W(2)/L(2)) * I_d);
    gm3 = sqrt(2 * UP_Cox * (W(3)/L(3)) * I_d);
    gm4 = sqrt(2 * UP_Cox * (W(4)/L(4)) * I_d);
    gm5 = sqrt(2 * UN_Cox * (W(5)/L(5)) * I_tail);
    gm6 = sqrt(2 * UN_Cox * (W(6)/L(6)) * I_d);
    gm7 = sqrt(2 * UP_Cox * (W(7)/L(7)) * I_d);
    gm8 = sqrt(2 * UP_Cox * (W(8)/L(8)) * I_d);
    
    % Output resistances
    ro1 = 1 / (lambda_n * I_d);
    ro2 = 1 / (lambda_p * I_d);
    ro3 = 1 / (lambda_p * I_d);
    ro4 = 1 / (lambda_n * I_d);
    ro5 = 1 / (lambda_n * I_tail);
    ro6 = 1 / (lambda_n * I_d);
    ro7 = 1 / (lambda_p * I_d);
    ro8 = 1 / (lambda_p * I_tail);
    
    % Gain stages
    A1 = gm1 * (ro1 || ro3);
    A2 = gm6 * (ro6 || ro7);
    Gain = A1 * A2;
    
    % Gain-Bandwidth Product
    GBW = gm1 / (2 * pi * Cc);
    
    % Phase Margin (simplified)
    p1 = 1 / (2 * pi * (ro6 || ro7) * Cc);
    p2 = gm6 / (2 * pi * Cc);
    PM = 90 - atand(GBW / p2);
    PM = max(PM, 10);
    
    % Slew Rate
    SR = I_tail / Cc;
    
    % Power dissipation
    Power = VDD * (I_tail + 2 * I_d);
    
    % Input Common Mode Range
    ICMR_min = VDD - abs(VTHP) - sqrt(I_d / (UP_Cox * (W(3)/L(3)))) - sqrt(I_d / (UN_Cox * (W(1)/L(1))));
    ICMR_max = VDD - abs(VTHP) - sqrt(I_d / (UP_Cox * (W(3)/L(3))));
    
    % Output Swing
    Vout_max = VDD - abs(VTHP) - sqrt(I_d / (UP_Cox * (W(7)/L(7))));
    Vout_min = VTHN + sqrt(I_d / (UN_Cox * (W(6)/L(6))));
    Output_swing = Vout_max - Vout_min;
    
    % ========== CONSTRAINTS ==========
    penalty = 0;
    
    % Constraint 1: DC Gain >= 5000 (74dB)
    if Gain < 5000
        penalty = penalty + 10 * ((5000 - Gain) / 1000)^2;
    end
    
    % Constraint 2: GBW >= 5 MHz
    if GBW < 5e6
        penalty = penalty + 10 * ((5e6 - GBW) / 1e6)^2;
    end
    
    % Constraint 3: Phase Margin >= 60 degrees
    if PM < 60
        penalty = penalty + 10 * ((60 - PM) / 20)^2;
    end
    
    % Constraint 4: Slew Rate >= 10 V/us
    if SR < 10e6
        penalty = penalty + 5 * ((10e6 - SR) / 5e6)^2;
    end
    
    % Constraint 5: Power <= 5 mW
    if Power > 5e-3
        penalty = penalty + 5 * ((Power - 5e-3) / 2e-3)^2;
    end
    
    % Constraint 6: Input Common Mode Range >= 1V
    if (ICMR_max - ICMR_min) < 1.0
        penalty = penalty + 5 * ((1.0 - (ICMR_max - ICMR_min)) / 0.5)^2;
    end
    
    % Constraint 7: Output Swing >= 1.2V
    if Output_swing < 1.2
        penalty = penalty + 5 * ((1.2 - Output_swing) / 0.5)^2;
    end
    
    % ========== COST FUNCTION ==========
    % Weights
    w_gain = 0.30;
    w_gbw = 0.25;
    w_pm = 0.20;
    w_power = 0.15;
    w_sr = 0.10;
    
    % Normalized objectives (minimization)
    cost = w_gain * (10000 / max(Gain, 1)) + ...
           w_gbw * (10e6 / max(GBW, 1e6)) + ...
           w_pm * (90 / max(PM, 1)) + ...
           w_sr * (20e6 / max(SR, 1e6)) + ...
           w_power * (Power / 5e-3) + ...
           1e6 * penalty;
    
    if ~isfinite(cost) || cost > 1e10
        cost = 1e10;
    end
end