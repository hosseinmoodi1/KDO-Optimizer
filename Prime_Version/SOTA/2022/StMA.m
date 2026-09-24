function [BestScore, BestPos, ConvergenceCurve] = StMA(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% StMA - Sterna Migration Algorithm
% -------------------------------------------------------------------------
% Reference:
%   Bai, H., Tong, W., Wei, B., Duan, C., & Liu, Q. (2025). A Sterna
%   migration algorithm-based efficient bionic engineering optimization
%   algorithm. Scientific Reports, 15, 38261.
%   https://doi.org/10.1038/s41598-025-22038-7
%
% Description:
%   A bio-inspired optimizer based on the collective migration behaviour of
%   the Oriental pratincole (Sterna). Key components:
%     - Multi-cluster sectoral diffusion initialization
%     - Leader-follower dynamic migration with role switching
%     - Three-phase search: Exploratory Launch / Decision Convergence /
%       Spatial Expansion
%     - Adaptive perturbation with spatial compression
%     - Multi-phase termination (max iter / stagnation / diversity)
%
% Inputs / Outputs: see SFO.m for the standard interface convention.
% -------------------------------------------------------------------------

    lb = lb(:)';  ub = ub(:)';
    FEs = 0;

    % ---- Parameters (Table 1 of the paper) ----
    NL_ratio     = 0.20;   % Fraction of leaders
    eps_stall    = 1e-8;   % Stagnation threshold
    tau_stag     = 20;     % Stagnation window
    theta_div    = 1e-4;   % Diversity threshold
    eta0         = 0.20;   % Initial perturbation scale
    lambda_eta   = 0.01;   % Perturbation decay rate
    k_sigmoid    = 10;     % Sigmoid slope for spatial compression
    p_perturb    = 0.20;   % Fraction of worst agents perturbed
    alpha_couple = 0.20;   % Perturbation-diversity coupling weight

    % ---- Update coefficients ----
    w_iner   = 0.5;   mu_lead  = 0.5;   sigma_p = 0.05;
    alpha_f  = 0.5;   beta_f   = 0.3;   delta_f = 0.05;

    NL = max(2, floor(NL_ratio * nPop));
    M  = 3;   % Number of initial clusters

    % ---- Multi-cluster sectoral diffusion initialization (Eq. 5-9) ----
    K = floor(nPop / M);
    X = zeros(nPop, dim);
    for j = 1:M
        Sj = lb + (0.2 + 0.6 * rand(1, dim)) .* (ub - lb);
        for k = 1:K
            idx = (j-1)*K + k;
            if idx > nPop, break; end
            rho = rand(1, dim) * 0.1;
            X(idx, :) = Sj + rho .* (ub - lb) .* (2*rand(1, dim) - 1);
        end
    end
    for idx = K*M+1 : nPop
        X(idx, :) = rand(1, dim) .* (ub - lb) + lb;
    end
    X = min(max(X, lb), ub);

    Fit = inf(nPop, 1);
    for i = 1:nPop
        if FEs >= MaxFEs, break; end
        Fit(i) = fobj(X(i, :));  FEs = FEs + 1;
    end
    [BestScore, idx] = min(Fit);
    BestPos = X(idx, :);

    X_prev           = X;
    StallCounter     = zeros(nPop, 1);
    BestHist         = BestScore * ones(1, tau_stag);
    ConvergenceCurve = zeros(1, MaxIter);
    sigma0           = std(X(:)) + eps;

    % ---- Main loop ----
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % Three-phase scheduling
        if t < MaxIter/3
            phase = 1;    % Exploratory Launch
        elseif t < 2*MaxIter/3
            phase = 2;    % Decision Convergence
        else
            phase = 3;    % Spatial Expansion
        end

        % Spatial compression factor (Eq. 25)
        beta_sig = 1 / (1 + exp(-k_sigmoid * (t - MaxIter/2)));

        % Population diversity (Eq. 29)
        W_mean  = mean(X, 1);
        sigma_t = sqrt(mean(sum((X - W_mean).^2, 2)));

        % Adaptive perturbation scale (Eq. 33)
        eta_min_p = 0.01;  eta_max_p = 0.30;
        eta_t = eta_min_p + (eta_max_p - eta_min_p) * (1 - sigma_t / (sigma0 + eps));

        % Identify leaders
        [~, sidx]  = sort(Fit);
        Leaders    = sidx(1:NL);
        C_center   = mean(X, 1);
        X_new      = X;

        % ---- Leader / Follower update ----
        for i = 1:nPop
            if ismember(i, Leaders)
                % Leader update (Eq. 10)
                X_new(i, :) = X(i, :) + w_iner * (X(i, :) - X_prev(i, :)) ...
                            + mu_lead * (BestPos - X(i, :)) ...
                            + sigma_p * randn(1, dim);
            else
                % Follower update (Eq. 12)
                dists = sum((X(Leaders, :) - X(i, :)).^2, 2);
                [~, nearest] = min(dists);
                Lk = X(Leaders(nearest), :);

                if phase == 1
                    % Broad exploration
                    X_new(i, :) = X(i, :) ...
                                + alpha_f * (Lk - X(i, :)) ...
                                + beta_f  * (C_center - X(i, :)) ...
                                + delta_f * (2*rand(1, dim) - 1);
                elseif phase == 2
                    % Decision convergence with spatial compression (Eq. 26)
                    R1 = 2*rand(1, dim) - 1;
                    X_new(i, :) = X(i, :) ...
                                + (1 - beta_sig) * R1 ...
                                + beta_sig * (BestPos - X(i, :)) ...
                                + eta_t * randn(1, dim);
                else
                    % Spatial re-expansion
                    X_new(i, :) = X(i, :) ...
                                + alpha_f * (Lk - X(i, :)) ...
                                + beta_f  * (C_center - X(i, :)) ...
                                + eta_t   * randn(1, dim);
                end
            end
            X_new(i, :) = min(max(X_new(i, :), lb), ub);
        end

        % ---- Evaluation with greedy selection ----
        for i = 1:nPop
            if FEs >= MaxFEs, break; end
            NewFit = fobj(X_new(i, :));  FEs = FEs + 1;
            if NewFit < Fit(i)
                X(i, :) = X_new(i, :);  Fit(i) = NewFit;
                StallCounter(i) = 0;
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = X_new(i, :);
                end
            else
                StallCounter(i) = StallCounter(i) + 1;
            end
        end

        % ---- Local perturbation of worst p% (Eq. 27-30) ----
        [~, sidx_desc] = sort(Fit, 'descend');
        n_perturb = max(1, floor(p_perturb * nPop));
        for k = 1:n_perturb
            i = sidx_desc(k);
            if FEs >= MaxFEs, break; end
            eta_p = eta0 * exp(-lambda_eta * t) ...
                  * (1 + alpha_couple * sigma_t / (sigma0 + eps));
            X(i, :) = X(i, :) + eta_p * randn(1, dim);
            X(i, :) = min(max(X(i, :), lb), ub);
            Fit(i)  = fobj(X(i, :));  FEs = FEs + 1;
            if Fit(i) < BestScore
                BestScore = Fit(i);  BestPos = X(i, :);
            end
        end

        X_prev = X;
        ConvergenceCurve(t) = BestScore;

        % ---- Multi-phase termination check (Eq. 16) ----
        BestHist = [BestHist(2:end), BestScore];
        if t > tau_stag
            stag_diff = max(BestHist) - min(BestHist);
            if stag_diff < eps_stall || sigma_t < theta_div
                ConvergenceCurve = ConvergenceCurve(1:t);
                break;
            end
        end
    end
end