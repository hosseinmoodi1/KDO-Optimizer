function params = extract_params_E1(x)
% extract_params_E1: Analytical estimate of physical parameters for E1
% Dimension: 18 variables (same as Electro_E1_Hard.m and template_E1.sp)
% UPDATED: L bounds aligned with cost function: [200e-9, 1e-6].

    params = struct();
    params.Gain_dB = NaN;
    params.GBW_Hz  = NaN;
    params.PM_deg  = NaN;
    params.SR_Vus  = NaN;
    params.Power_W = NaN;
    params.ICMR_V  = NaN;
    params.Output_swing_V = NaN;

    VDD = 1.8; VTHN=0.45; VTHP=-0.45;
    UN_Cox=180e-6; UP_Cox=60e-6;
    lambda_n=0.05; lambda_p=0.05;

    % Clamp variables (aligned with Electro_E1_Hard.m)
    W = min(max(abs(x(1:8)), 1e-6), 100e-6);
    L = min(max(abs(x(9:16)), 200e-9), 1e-6);   % FIXED: min 200e-9 (was 500e-9)
    Ibias = min(max(abs(x(17)), 1e-6), 500e-6);
    Cc = min(max(abs(x(18)), 0.1e-12), 10e-12);

    W1=W(1); L1=L(1); W2=W(2); L2=L(2); W3=W(3); L3=L(3);
    W4=W(4); L4=L(4); W5=W(5); L5=L(5); W6=W(6); L6=L(6);
    W7=W(7); L7=L(7); W8=W(8); L8=L(8);

    I_tail = Ibias;
    I_d = I_tail/2;
    I_stage2 = I_tail/2;
    I_bias_ref = I_tail/2;

    Vov1 = sqrt(2*I_d/(UN_Cox*(W1/L1)+eps));
    Vov3 = sqrt(2*I_d/(UP_Cox*(W3/L3)+eps));
    gm1 = 2*I_d/max(Vov1,eps);
    gm3 = 2*I_d/max(Vov3,eps);
    ro1 = 1/(lambda_n*I_d+eps);
    ro3 = 1/(lambda_p*I_d+eps);
    Rout1 = ro1*ro3/(ro1+ro3+eps);
    A1 = gm1*Rout1;

    Vov6 = sqrt(2*I_stage2/(UN_Cox*(W6/L6)+eps));
    gm6 = 2*I_stage2/max(Vov6,eps);
    ro6 = 1/(lambda_n*I_stage2+eps);
    ro7 = 1/(lambda_p*I_stage2+eps);
    Rout2 = ro6*ro7/(ro6+ro7+eps);
    A2 = gm6*Rout2;
    Gain = A1*A2;

    GBW = gm1/(2*pi*Cc+eps);
    p2 = gm6/(2*pi*Cc+eps);
    PM = 90 - atand(GBW / max(p2, eps));

    SR = I_tail/(Cc+eps);
    Power = VDD * (I_tail + I_d + I_stage2 + I_bias_ref);

    Vov1_t = sqrt(I_d/(UN_Cox*(W1/L1)+eps));
    Vov3_t = sqrt(I_d/(UP_Cox*(W3/L3)+eps));
    ICMR_min = VDD - abs(VTHP) - Vov3_t - Vov1_t;
    ICMR_max = VDD - abs(VTHP) - Vov3_t;
    ICMR = ICMR_max - ICMR_min;

    Vov6_t = sqrt(I_stage2/(UN_Cox*(W6/L6)+eps));
    Vov7_t = sqrt(I_stage2/(UP_Cox*(W7/L7)+eps));
    Vout_max = VDD - abs(VTHP) - Vov7_t;
    Vout_min = VTHN + Vov6_t;
    Output_swing = Vout_max - Vout_min;

    params.Gain_dB = 20*log10(max(Gain,1));
    params.GBW_Hz = GBW;
    params.PM_deg = PM;
    params.SR_Vus = SR/1e6;
    params.Power_W = Power;
    params.ICMR_V = ICMR;
    params.Output_swing_V = Output_swing;
end