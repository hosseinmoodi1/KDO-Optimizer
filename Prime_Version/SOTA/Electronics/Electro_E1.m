function cost = Electro_E1(x)
% =========================================================
% E1: CMOS Two-Stage Operational Amplifier Optimization
% =========================================================
% Objective : Minimize power while satisfying gain and GBW specs
%             (formulated as minimization problem)
%
% Design Variables (10D):
%   x(1..8) = W1..W8 : transistor widths  [1e-6 .. 100e-6] m
%             W1,W2  = differential input pair (NMOS)
%             W3,W4  = active load (PMOS)
%             W5,W6  = second stage
%             W7     = tail current mirror
%             W8     = output load transistor
%   x(9)    = Ibias  : tail bias current   [1e-6 .. 500e-6] A
%   x(10)   = Cc     : Miller compensation capacitor [0.1p .. 10p] F
%
% Performance Specs (constraints):
%   Gain  >= 1000  (60 dB)
%   GBW   >= 1 MHz
%   Power <= 5 mW
%
% Technology: 180 nm CMOS, VDD = 1.8 V
%
% Reference: Hershenson et al., IEEE TCAD, 2001
% =========================================================

    % --- Technology parameters ---
    VDD  = 1.8;       % Supply voltage (V)
    mu_n = 400e-4;    % NMOS carrier mobility (m^2/V·s)
    mu_p = 200e-4;    % PMOS carrier mobility (m^2/V·s)
    Cox  = 10e-3;     % Gate oxide capacitance per unit area (F/m^2)
    L    = 180e-9;    % Minimum channel length, fixed at 180 nm

    % --- Extract and sanitize design variables ---
    W     = x(1:8);
    Ibias = max(abs(x(9)),  1e-9);   % bias current (A), clipped to avoid /0
    Cc    = max(abs(x(10)), 1e-15);  % compensation capacitor (F)

    % --- Input pair transconductance (M1, M2 — NMOS differential pair) ---
    gm1 = sqrt(2 * mu_n * Cox * (W(1)/L) * Ibias);

    % --- Output resistance (Early voltage model) ---
    VA  = 10;                          % Early voltage approximation (V)
    ro2 = VA / max(Ibias, 1e-9);       % output resistance of M2
    ro4 = VA / max(Ibias, 1e-9);       % output resistance of M4 (load)
    Rout = (ro2 * ro4) / (ro2 + ro4);  % parallel combination

    % --- DC Gain (open-loop) ---
    Gain = max(gm1 * Rout, 1e-6);      % V/V, clipped to avoid log(0)

    % --- Gain-Bandwidth Product ---
    GBW = max(gm1 / (2 * pi * Cc), 1e-6);  % Hz

    % --- Power dissipation ---
    Power = VDD * Ibias;   % W (static power, dominant term)

    % --- Soft penalty for constraint violations ---
    penalty = 0;
    if Gain < 1000                          % Gain >= 60 dB
        penalty = penalty + ((1000  - Gain)  / 1000 )^2;
    end
    if GBW < 1e6                            % GBW >= 1 MHz
        penalty = penalty + ((1e6   - GBW)   / 1e6  )^2;
    end
    if Power > 5e-3                         % Power <= 5 mW
        penalty = penalty + ((Power - 5e-3)  / 5e-3 )^2;
    end

    % --- Cost function (minimization) ---
    % Terms: minimize 1/Gain, 1/GBW (scaled), Power, and constraint violations
    cost = (1/Gain) + (1/GBW)*1e-6 + Power*200 + 1e4*penalty;

end
