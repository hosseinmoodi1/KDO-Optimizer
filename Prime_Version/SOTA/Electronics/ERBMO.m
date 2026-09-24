function [BestScore, BestPos, ConvergenceCurve] = ERBMO(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% ERBMO - Enhanced Red-billed Blue Magpie Optimizer
% -------------------------------------------------------------------------
% Reference:
%   Wang, H., Xin, Z., Qi, X., Zhou, Q., Liu, J., & Li, Z. (2026). Enhanced
%   Red-billed Blue Magpie Optimizer for engineering optimization problems.
%   Scientific Reports, 16, 10619.
%   https://doi.org/10.1038/s41598-026-44507-3
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

    % ---- RBMO group-size parameters (safe for small nPop) ----
    p_min = min(2, nPop);
    p_max = min(5, nPop);
    if p_max < p_min, p_max = p_min; end
    q_min = min(10, nPop);
    q_max = nPop;
    if q_max < q_min, q_min = q_max; end
    CF_init = 0.8;

    % ---- ERBMO adaptation weights ----
    alpha_cf = 0.95;
    beta_cf  = 0.05;
    gamma_cf = 0.01;

    % ---- Periodic pattern search parameters ----
    delta_step_init = 0.5;
    b_acc           = 1.2;
    a_step          = 0.5;
    alpha_min       = 1e-3;

    % ---- Evolutionary combinatorial mutation ----
    p_mut_init  = 0.1;
    dp_mut      = 0.001;
    p_mut_max   = 0.5;
    p_gauss     = 0.6;
    sigma_gauss = 0.1;

    % ---- Initialization ----
    X   = rand(nPop, dim) .* (ub - lb) + lb;
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

    CF               = CF_init;
    PrevBest         = BestScore;
    delta_local      = delta_step_init;
    ConvergenceCurve = zeros(1, MaxIter);

    % ===== Main loop =====
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % ===== Diversity-adaptive CF update (robust to Inf Fit) =====
        improvement = abs(PrevBest - BestScore);
        finiteFit   = Fit(isfinite(Fit));
        if isempty(finiteFit)
            diversity = 0;
        else
            diversity = std(finiteFit);
        end
        CF = alpha_cf * CF ...
           + beta_cf  * (1 / (1 + improvement^(1 + rand()))) ...
           + gamma_cf * diversity;
        CF = min(max(CF, 0.01), 1.5);

        do_pattern = (mod(t, 100) == 0);
        p_current  = min(p_mut_init + dp_mut * (t - 1), p_mut_max);

        % ===== RBMO core =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            Xnew = X(i, :);
            if rand() < 0.5
                if rand() < 0.5
                    p_size = randi([p_min, p_max]);
                    idxs   = randperm(nPop, p_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xrs    = X(randi(nPop), :);
                    Xnew   = X(i, :) + (Xm - Xrs) * rand();
                else
                    q_size = randi([q_min, q_max]);
                    idxs   = randperm(nPop, q_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xrs    = X(randi(nPop), :);
                    Xnew   = X(i, :) + (Xm - Xrs) * rand();
                end
            else
                if rand() < 0.5
                    p_size = randi([p_min, p_max]);
                    idxs   = randperm(nPop, p_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xnew   = BestPos + CF * (Xm - X(i, :)) * randn();
                else
                    q_size = randi([q_min, q_max]);
                    idxs   = randperm(nPop, q_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xnew   = BestPos + CF * (Xm - X(i, :)) * randn();
                end
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

        % ===== Evolutionary combinatorial mutation =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end
            if rand() < p_current
                if rand() < p_gauss
                    X_mut = X(i, :) + sigma_gauss * randn(1, dim);
                else
                    X_mut = lb + (ub - lb) .* rand(1, dim);
                end
                X_mut = min(max(X_mut, lb), ub);
                f_mut = fobj(X_mut);  FEs = FEs + 1;
                if f_mut < Fit(i)
                    X(i, :) = X_mut;  Fit(i) = f_mut;
                    if f_mut < BestScore
                        BestScore = f_mut;  BestPos = X_mut;
                    end
                end
            end
        end

        % ===== Periodic pattern search =====
        if do_pattern && FEs < MaxFEs
            x1     = BestPos;
            x_curr = BestPos;
            improved_explore = false;

            for j = 1:dim
                if FEs >= MaxFEs, break; end
                Ij = zeros(1, dim);  Ij(j) = 1;

                y_plus = min(max(x_curr + delta_local * Ij, lb), ub);
                f_plus = fobj(y_plus);  FEs = FEs + 1;
                if f_plus < BestScore
                    x_curr = y_plus;  BestPos = y_plus;  BestScore = f_plus;
                    improved_explore = true;
                else
                    if FEs >= MaxFEs, break; end
                    y_minus = min(max(x_curr - delta_local * Ij, lb), ub);
                    f_minus = fobj(y_minus);  FEs = FEs + 1;
                    if f_minus < BestScore
                        x_curr = y_minus;  BestPos = y_minus;  BestScore = f_minus;
                        improved_explore = true;
                    end
                end
            end

            if improved_explore && FEs < MaxFEs
                x = x_curr;
                y_pattern = x + b_acc * (x - x1);
                y_pattern = min(max(y_pattern, lb), ub);
                f_pattern = fobj(y_pattern);  FEs = FEs + 1;
                if f_pattern < BestScore
                    BestPos = y_pattern;  BestScore = f_pattern;
                end
            end

            if delta_local < alpha_min
                delta_local = delta_local * a_step;
            end
        end

        PrevBest = BestScore;
        ConvergenceCurve(t) = BestScore;
    end
end