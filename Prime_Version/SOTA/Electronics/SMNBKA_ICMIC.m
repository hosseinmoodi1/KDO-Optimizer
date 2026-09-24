function [BestScore, BestPos, ConvergenceCurve] = SMNBKA_ICMIC(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% SMNBKA-ICMIC - Binary Evolution Enhanced Black-kite Algorithm
% -------------------------------------------------------------------------
% Reference:
%   Sun, H., Tang, N., Li, Z., & Chen, H. (2026). Based on binary evolution
%   operator-enhanced black-kite algorithm with natural replacement for
%   engineering numerical optimization problems. Scientific Reports, 16,
%   6881. https://doi.org/10.1038/s41598-026-35846-2
%
% Inputs / Outputs: standard SFO.m / INFO.m interface convention.
% -------------------------------------------------------------------------

    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';  ub = ub(:)';
    end

    FEs = 0;

    % ---- Parameters (Table 4 of the paper) ----
    a_icmic        = 0.7;
    b_icmic        = 1.5;
    p_attack       = 0.9;
    eta_c          = 1.8;
    gamma_replace  = 0.20;
    natural_radius = 0.10;

    % ===== ICMIC chaotic-map initialization =====
    X = zeros(nPop, dim);
    for i = 1:nPop
        c = rand(1, dim);
        for k = 1:100
            c = a_icmic * sin(pi * c) + b_icmic * (c - abs(c));
        end
        c = abs(c - floor(c));
        X(i, :) = lb + c .* (ub - lb);
    end
    X = min(max(X, lb), ub);

    Fit = inf(nPop, 1);
    for i = 1:nPop
        if FEs >= MaxFEs, break; end
        Fit(i) = fobj(X(i, :));  FEs = FEs + 1;
    end

    if all(~isfinite(Fit))
        BestScore = inf;  BestPos = X(1, :);
        ConvergenceCurve = [];
        return;
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
        n_factor = 0.05 * exp(-2 * (t / MaxIter)^2);

        % ===== Attack phase with SBX =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

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
            C_cauchy = tan(pi * (rand(1, dim) - 0.5));

            if Fit(i) < Fit(r_idx)
                Xnew = X(i, :) + C_cauchy .* (X(i, :) - BestPos);
            else
                m_factor = 2 * sin(rand() + pi/2);
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

        % ===== Natural replacement mechanism =====
        [~, sidx] = sort(Fit);
        n_elim = max(1, floor(gamma_replace * nPop));
        for k = 1:n_elim
            i = sidx(end - k + 1);
            if FEs >= MaxFEs, break; end
            r_scalar = rand();
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