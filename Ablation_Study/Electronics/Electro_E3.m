function cost = Electro_E3(x)
% =========================================================
% Problem E3: CMOS LC Voltage-Controlled Oscillator (LC-VCO)
% =========================================================
% Objective (minimization):
%   - Maximize oscillation frequency (f_osc) and tuning range (TR)
%   - Minimize phase noise (PN) at 1 MHz offset
%   - Minimize DC power consumption
%
% Design Variables (8D):
%   x(1) = L      : Tank inductance          [0.5e-9 .. 10e-9]  H
%   x(2) = Cvar   : Varactor capacitance     [50e-15 .. 2e-12]  F
%   x(3) = Cfix   : Fixed tank capacitance   [10e-15 .. 1e-12]  F
%   x(4) = Wn     : NMOS cross-coupled width [1e-6   .. 200e-6] m
%   x(5) = Wp     : PMOS cross-coupled width [1e-6   .. 200e-6] m
%   x(6) = Ibias  : Tail bias current        [10e-6  .. 10e-3]  A
%   x(7) = Rbias  : Bias resistor            [100    .. 50e3]   Ohm
%   x(8) = Cbuf   : Output buffer capacitance[10e-15 .. 1e-12]  F
%
% Soft Constraints:
%   f_osc in [1 GHz, 5 GHz]
%   PN(1 MHz offset) <= -100 dBc/Hz
%   Power <= 10 mW
%
% Suggested bounds for optimizer:
%   lb = [0.5e-9, 50e-15, 10e-15, 1e-6,   1e-6,   10e-6, 100,  10e-15]
%   ub = [10e-9,  2e-12,  1e-12,  200e-6, 200e-6, 10e-3, 50e3, 1e-12 ]
% =========================================================

    % -------------------------------------------------------
    % Technology & Physical Constants
    % -------------------------------------------------------
    VDD      = 1.2;       % Supply voltage (V)
    gamma    = 1.5;       % MOSFET channel noise factor (short-channel)
    k_B      = 1.38e-23;  % Boltzmann constant (J/K)
    T        = 300;       % Ambient temperature (K)
    QL       = 10;        % Inductor loaded quality factor (dimensionless)
    Cox_n    = 8e-3;      % NMOS gate oxide capacitance per area (F/m^2)
    Cox_p    = 8e-3;      % PMOS gate oxide capacitance per area (F/m^2)
    Lch      = 180e-9;    % Channel length (m), 180 nm technology node
    mu_n_Cox = 270e-4;    % NMOS process transconductance (A/V^2), mu_n*Cox
    mu_p_Cox = 90e-4;     % PMOS process transconductance (A/V^2), mu_p*Cox

    % -------------------------------------------------------
    % Extract & Clip Design Variables (prevent degenerate values)
    % -------------------------------------------------------
    L     = max(abs(x(1)), 0.5e-9);   % Tank inductance (H)
    Cvar  = max(abs(x(2)), 50e-15);   % Varactor capacitance (F)
    Cfix  = max(abs(x(3)), 10e-15);   % Fixed capacitance (F)
    Wn    = max(abs(x(4)), 1e-6);     % NMOS width (m)
    Wp    = max(abs(x(5)), 1e-6);     % PMOS width (m)
    Ibias = max(abs(x(6)), 10e-6);    % Bias current (A)
    Rbias = max(abs(x(7)), 100);      % Bias resistor (Ohm)
    Cbuf  = max(abs(x(8)), 10e-15);   % Buffer load capacitance (F)

    % -------------------------------------------------------
    % LC Tank: Total Capacitance & Oscillation Frequency
    % -------------------------------------------------------
    % Buffer capacitance adds parasitic load to the tank
    Ctank   = Cvar + Cfix + 0.1*Cbuf;              % effective tank cap (F)
    omega_0 = 1 / sqrt(max(L * Ctank, 1e-24));     % resonant frequency (rad/s)
    f_osc   = omega_0 / (2*pi);                     % oscillation frequency (Hz)

    % -------------------------------------------------------
    % Transconductance of Cross-Coupled Pair
    % -------------------------------------------------------
    % Each half-circuit carries Ibias/2; gm = sqrt(2 * mu*Cox * W/L * Id)
    Id_half = Ibias / 2;                            % drain current per transistor (A)
    gm_n = sqrt(max(2 * mu_n_Cox * (Wn/Lch) * Id_half, 0)); % NMOS gm (A/V)
    gm_p = sqrt(max(2 * mu_p_Cox * (Wp/Lch) * Id_half, 0)); % PMOS gm (A/V)
    gm_total = gm_n + gm_p;                         % total effective gm (A/V)

    % -------------------------------------------------------
    % Phase Noise (Leeson's Model, extended with gm)
    % -------------------------------------------------------
    % Offset frequency for PN specification
    f_offset = 1e6;                                 % 1 MHz offset (Hz)

    % Tank signal power: P_sig ~ (Ibias * Rtank)^2 / 2
    % Rtank at resonance ~ QL * omega_0 * L
    Rtank  = QL * omega_0 * L;                      % tank parallel resistance (Ohm)
    P_sig  = 0.5 * (Ibias * Rtank)^2 / max(Rtank, 1); % signal power proxy (W)
    P_sig  = max(P_sig, 1e-20);                     % guard against zero

    % Leeson's single-sideband phase noise (linear scale, dBc/Hz)
    % L(f_offset) = (gamma * 4kT / P_sig) * (f_osc / (2*QL*f_offset))^2
    PN_lin = (gamma * 4 * k_B * T / P_sig) * (f_osc / (2 * QL * f_offset))^2;
    PN_dBc = 10 * log10(max(PN_lin, 1e-30));        % phase noise (dBc/Hz)

    % -------------------------------------------------------
    % DC Power Consumption
    % -------------------------------------------------------
    Power = VDD * Ibias;                            % total DC power (W)

    % -------------------------------------------------------
    % Tuning Range (varactor swept ±30% around nominal)
    % -------------------------------------------------------
    Cvar_min = 0.7 * Cvar;
    Cvar_max = 1.3 * Cvar;
    f_lo = 1 / (2*pi * sqrt(max(L*(Cvar_max + Cfix), 1e-24))); % lowest freq (Hz)
    f_hi = 1 / (2*pi * sqrt(max(L*(Cvar_min + Cfix), 1e-24))); % highest freq (Hz)
    TR   = f_hi - f_lo;                             % absolute tuning range (Hz)

    % -------------------------------------------------------
    % Oscillation Start-Up Condition Check
    % -------------------------------------------------------
    % VCO starts if gm_total > 1/Rtank (negative resistance condition)
    % Soft penalty if start-up margin is insufficient
    startup_margin = gm_total * Rtank;              % should be > 1 for oscillation
    if startup_margin < 1.0
        startup_penalty = (1.0 - startup_margin)^2;
    else
        startup_penalty = 0;
    end

    % -------------------------------------------------------
    % Soft Constraint Penalties
    % -------------------------------------------------------
    penalty = startup_penalty;

    % Constraint 1: f_osc in [1 GHz, 5 GHz]
    f_min_req = 1e9;
    f_max_req = 5e9;
    if f_osc < f_min_req
        penalty = penalty + ((f_min_req - f_osc) / f_min_req)^2;
    elseif f_osc > f_max_req
        penalty = penalty + ((f_osc - f_max_req) / f_max_req)^2;
    end

    % Constraint 2: Phase noise <= -100 dBc/Hz at 1 MHz offset
    PN_req = -100;
    if PN_dBc > PN_req
        penalty = penalty + ((PN_dBc - PN_req) / abs(PN_req))^2;
    end

    % Constraint 3: Power <= 10 mW
    Pmax = 10e-3;
    if Power > Pmax
        penalty = penalty + ((Power - Pmax) / Pmax)^2;
    end

    % -------------------------------------------------------
    % Multi-Objective Cost (weighted sum, minimization)
    % -------------------------------------------------------
    % Weights: frequency (0.3), tuning range (0.2),
    %          phase noise (0.3), power (0.2)
    w_f  = 0.3;
    w_TR = 0.2;
    w_PN = 0.3;
    w_P  = 0.2;

    % Normalize each term to be dimensionless and O(1) at typical values
    cost = w_f  * (1 / max(f_osc, 1e6))                    + ...  % maximize f_osc
           w_TR * (1 / max(TR,    1e6))                    + ...  % maximize TR
           w_PN * ((PN_dBc - PN_req + 200) / 200)          + ...  % minimize PN (scaled)
           w_P  * (Power / Pmax)                           + ...  % minimize power
           1e4  * penalty;                                         % constraint violations

    % Guard: return large finite value instead of NaN/Inf
    if ~isfinite(cost)
        cost = 1e10;
    end

end
