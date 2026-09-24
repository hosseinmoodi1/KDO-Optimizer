function [BestScore, BestPos, ConvergenceCurve] = CAO(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% CAO - Cell Autophagy Optimization Algorithm
% -------------------------------------------------------------------------
% Reference:
%   Al Tawil, A., Edinat, A., Chelloug, S. A., & Fathi, H. (2026).
%   Cell autophagy optimization algorithm based on cellular self-degradation
%   and recycling for global optimization. Scientific Reports.
%   https://doi.org/10.1038/s41598-026-65895-6
%
% Description:
%   A metaheuristic inspired by the cellular autophagy mechanism. Seven
%   computational phases:
%     1. Chaotic tent-map population seeding
%     2. Opposition-based learning (OBL)
%     3. Stress detection and autophagy triggering
%     4. Autophagosome formation
%     5. Lysosomal degradation
%     6. Lévy-flight recycling from the elite archive
%     7. Mitophagy-inspired escape from stagnation
%
% Inputs / Outputs: see SFO.m for the standard interface convention.
% -------------------------------------------------------------------------

    lb = lb(:)';  ub = ub(:)';
    FEs = 0;

    % ---- Parameters (Table 7 of the paper) ----
    sigma_max   = 0.9;    % Initial stress threshold
    kappa       = 3.5;    % Threshold decay rate
    pdeg0       = 0.92;   % Initial degradation probability
    pdegT       = 0.25;   % Final degradation probability
    beta_levy   = 1.5;    % Lévy exponent
    n_elite     = 5;      % Elite archive size
    tau_OBL     = 20;     % OBL refresh interval (iterations)
    rho_OBL     = 0.30;   % OBL refresh fraction (worst agents)
    kappa_stag  = 10;     % Stagnation threshold for mitophagy
    beta_osc    = 0.5;    % Oscillation amplitude
    gamma_osc   = 1.0;    % Oscillation frequency

    % ===== Phase 1: Chaotic tent-map seeding (Eq. 1) =====
    X = zeros(nPop, dim);
    for i = 1:nPop
        c = rand(1, dim);
        for k = 1:50
            c = tent_map(c);
        end
        X(i, :) = lb + c .* (ub - lb);
    end

    % ===== Phase 2: Full opposition-based learning (Eq. 2) =====
    X_opp  = lb + ub - X;
    AllX   = [X; X_opp];
    AllFit = inf(2*nPop, 1);
    for i = 1:2*nPop
        if FEs >= MaxFEs, break; end
        AllFit(i) = fobj(AllX(i, :));  FEs = FEs + 1;
    end
    [~, sortidx] = sort(AllFit);
    X = AllX(sortidx(1:nPop), :);

    Fit = inf(nPop, 1);
    for i = 1:nPop
        if FEs >= MaxFEs, break; end
        Fit(i) = fobj(X(i, :));  FEs = FEs + 1;
    end
    [BestScore, idx] = min(Fit);
    BestPos = X(idx, :);

    EliteX           = X;
    StallCounter     = 0;
    ConvergenceCurve = zeros(1, MaxIter);
    range_norm       = norm(ub - lb) + eps;

    % ===== Main loop =====
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end
        p = t / MaxIter;

        % Adaptive stress threshold (Eq. 3)
        sigma_th = sigma_max * exp(-kappa * p);

        % Time-varying degradation probability (Eq. 5)
        pdeg = pdeg0 - (pdeg0 - pdegT) * p;

        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % ===== Phase 3: Stress detection (spatial distance to best) =====
            sigma_i = norm(X(i, :) - BestPos) / range_norm;
            sigma_i = min(max(sigma_i, 0), 1);

            if sigma_i > sigma_th
                % ----- Phases 4-6: Degrade and rebuild via Lévy recycling -----
                r = rand(1, dim);
                D = X(i, :) .* (pdeg * r);                    % Eq. 5 (degrade)

                elite_idx = randi(size(EliteX, 1));
                x_elite   = EliteX(elite_idx, :);
                alpha_step = rand();
                L_scale = 0.01 * (ub - lb) * (1 - p)^2;       % Eq. 6 (scale)
                L = levy_step(dim, beta_levy) .* L_scale;
                Xnew = D + alpha_step * (x_elite - D) + L;    % Eq. 6 (rebuild)
            else
                % ----- Normal metabolism (lightweight exploitation) -----
                L_scale = 0.01 * (ub - lb) * (1 - p);
                L = levy_step(dim, beta_levy) .* L_scale;
                Xnew = X(i, :) + rand(1, dim) .* (BestPos - X(i, :)) + L;
            end

            Xnew   = min(max(Xnew, lb), ub);
            NewFit = fobj(Xnew);  FEs = FEs + 1;

            if NewFit < Fit(i)
                X(i, :) = Xnew;  Fit(i) = NewFit;
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = Xnew;
                    StallCounter = 0;
                else
                    StallCounter = StallCounter + 1;
                end
            else
                StallCounter = StallCounter + 1;
            end
        end

        % ---- Update elite archive ----
        [~, sidx] = sort(Fit);
        EliteX = X(sidx(1:min(n_elite, nPop)), :);

        % ===== Dynamic OBL refresh every tau_OBL iterations =====
        if mod(t, tau_OBL) == 0 && FEs < MaxFEs
            [~, sidx_desc] = sort(Fit, 'descend');
            n_worst = max(1, floor(rho_OBL * nPop));
            for k = 1:n_worst
                i = sidx_desc(k);
                if FEs >= MaxFEs, break; end
                X_opp_i = min(max(lb + ub - X(i, :), lb), ub);
                F_opp   = fobj(X_opp_i);  FEs = FEs + 1;
                if F_opp < Fit(i)
                    X(i, :) = X_opp_i;  Fit(i) = F_opp;
                    if F_opp < BestScore
                        BestScore = F_opp;  BestPos = X_opp_i;
                    end
                end
            end
        end

        % ===== Phase 7: Mitophagy-inspired escape (Eq. 7) =====
        if StallCounter >= kappa_stag && FEs < MaxFEs
            for i = 1:nPop
                if FEs >= MaxFEs, break; end
                L_scale = 0.05 * (ub - lb) * (1 - p);
                L = levy_step(dim, beta_levy) .* L_scale;
                X(i, :) = X(i, :) + beta_osc * sin(gamma_osc * t) + L;
                X(i, :) = min(max(X(i, :), lb), ub);
                Fit(i)  = fobj(X(i, :));  FEs = FEs + 1;
                if Fit(i) < BestScore
                    BestScore = Fit(i);  BestPos = X(i, :);
                end
            end
            StallCounter = 0;
        end

        ConvergenceCurve(t) = BestScore;
    end
end

% -------------------------------------------------------------------------
% Helper: Tent map (Eq. 1)
% -------------------------------------------------------------------------
function c = tent_map(c)
    mu1 = 1 / 0.7;
    mu2 = 1 / 0.3;
    idx = c < 0.7;
    c(idx)  = mu1 * c(idx);
    c(~idx) = mu2 * (1 - c(~idx));
end

% -------------------------------------------------------------------------
% Helper: Mantegna's Lévy step generator
% -------------------------------------------------------------------------
function L = levy_step(d, beta)
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
             (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    L = u ./ (abs(v).^(1/beta));
end