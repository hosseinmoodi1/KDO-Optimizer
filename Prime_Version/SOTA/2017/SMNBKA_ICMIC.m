function [BestScore, BestPos, ConvergenceCurve] = SMNBKA_ICMIC(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% SMNBKA-ICMIC - Binary Evolution Enhanced Black-kite Algorithm
% -------------------------------------------------------------------------
% Reference:
%   Sun, H., Tang, N., Li, Z., & Chen, H. (2026). Based on binary evolution
%   operator-enhanced black-kite algorithm with natural replacement for
%   engineering numerical optimization problems. Scientific Reports, 16,
%   6881. https://doi.org/10.1038/s41598-026-35846-2
%
% Description:
%   Enhanced Black-winged Kite Algorithm (BKA) with four key improvements:
%     1. ICMIC chaotic-map initialization for improved population diversity
%     2. Simulated Binary Crossover (SBX) in the attack phase
%     3. Globalized leader/follower migration with Cauchy mutation
%     4. Natural replacement mechanism (Darwinian culling of worst 20%)
%
% Inputs / Outputs: see SFO.m for the standard interface convention.
% -------------------------------------------------------------------------

    lb = lb(:)';  ub = ub(:)';
    FEs = 0;

    % ---- Parameters (Table 4 of the paper) ----
    a_icmic        = 0.7;    % ICMIC map parameter a
    b_icmic        = 1.5;    % ICMIC map parameter b
    p_attack       = 0.9;    % Attack probability (from BKA)
    eta_c          = 1.8;    % SBX distribution index
    gamma_replace  = 0.20;   % Fraction eliminated per generation
    natural_radius = 0.10;   % Neighborhood radius around best

    % ===== ICMIC chaotic-map initialization (Eq. 9-10) =====
    X = zeros(nPop, dim);
    for i = 1:nPop
        c = rand(1, dim);
        for k = 1:100
            c = a_icmic * sin(pi * c) + b_icmic * (c - abs(c));
        end
        c = abs(c - floor(c));           % Map to [0, 1)
        X(i, :) = lb + c .* (ub - lb);
    end
    X = min(max(X, lb), ub);

    Fit = inf(nPop, 1);
    for i = 1:nPop
        if FEs >= MaxFEs, break; end
        Fit(i) = fobj(X(i, :));  FEs = FEs + 1;
    end
    [BestScore, idx] = min(Fit);
    BestPos = X(idx, :);

    ConvergenceCurve = zeros(1, MaxIter);

    % ===== Main loop =====
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end
        n_factor = 0.05 * exp(-2 * (t / MaxIter)^2);   % Dynamic factor (Eq. 5)

        % ===== Attack phase with SBX operator =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % SBX evolution factor (Eq. 11)
            u = rand();
            if u <= 0.5
                beta_sbx = (2 * u)^(1 / (eta_c + 1));
            else
                beta_sbx = (1 / (2 * (1 - u)))^(1 / (eta_c + 1));
            end
            EF = 0.5 * ((1 + beta_sbx) * BestPos + (1 - beta_sbx) * X(i, :)) - X(i, :);

            r_val = rand();
            if p_attack < r_val
                Xnew = X(i, :) + n_factor * (1 + sin(r_val)) * X(i, :) + EF;
            else
                Xnew = X(i, :) + n_factor * (2 * r_val - 1)    * X(i, :) + EF;
            end
            Xnew   = min(max(Xnew, lb), ub);
            NewFit = fobj(Xnew);  FEs = FEs + 1;

            if NewFit < Fit(i)
                X(i, :) = Xnew;  Fit(i) = NewFit;
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = Xnew;
                end
            end
        end

        % ===== Migration phase with Cauchy mutation =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end
            r_idx = randi(nPop);
            % Cauchy mutation C(0, 1) via inverse-CDF
            C_cauchy = tan(pi * (rand(1, dim) - 0.5));

            if Fit(i) < Fit(r_idx)
                Xnew = X(i, :) + C_cauchy .* (X(i, :) - BestPos);
            else
                m_factor = 2 * sin(rand() + pi/2);     % Eq. 7
                Xnew = X(i, :) + C_cauchy .* (BestPos - m_factor * X(i, :));
            end
            Xnew   = min(max(Xnew, lb), ub);
            NewFit = fobj(Xnew);  FEs = FEs + 1;

            if NewFit < Fit(i)
                X(i, :) = Xnew;  Fit(i) = NewFit;
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = Xnew;
                end
            end
        end

        % ===== Natural replacement mechanism (Eq. 15-16) =====
        [~, sidx] = sort(Fit);
        n_elim = max(1, floor(gamma_replace * nPop));
        for k = 1:n_elim
            i = sidx(end - k + 1);
            if FEs >= MaxFEs, break; end
            r_scalar = rand();                                      % Scalar radius factor
            X(i, :)  = BestPos + natural_radius * r_scalar * (ub - lb);
            X(i, :)  = min(max(X(i, :), lb), ub);
            Fit(i)   = fobj(X(i, :));  FEs = FEs + 1;
            if Fit(i) < BestScore
                BestScore = Fit(i);  BestPos = X(i, :);
            end
        end

        ConvergenceCurve(t) = BestScore;
    end
end