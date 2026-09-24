function [BestScore, BestPos, ConvergenceCurve] = ERBMO(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% ERBMO - Enhanced Red-billed Blue Magpie Optimizer
% -------------------------------------------------------------------------
% Reference:
%   Wang, H., Xin, Z., Qi, X., Zhou, Q., Liu, J., & Li, Z. (2026). Enhanced
%   Red-billed Blue Magpie Optimizer for engineering optimization problems.
%   Scientific Reports, 16, 10619.
%   https://doi.org/10.1038/s41598-026-44507-3
%
% Description:
%   Enhanced variant of the Red-billed Blue Magpie Optimizer (RBMO). The
%   base RBMO core comprises four sequential phases: food searching (small
%   group + cluster), prey attacking (small group + cluster), and food
%   caching (greedy selection). ERBMO adds three synergistic mechanisms:
%     - Diversity-adaptive weight updating for the cooperation factor CF
%     - Periodic pattern search (every 100 iterations)
%     - Evolutionary probabilistic combinatorial mutation
%
% Inputs / Outputs: see SFO.m for the standard interface convention.
% -------------------------------------------------------------------------

    lb = lb(:)';  ub = ub(:)';
    FEs = 0;

    % ---- RBMO group-size parameters ----
    p_min = 2;   p_max = 5;      % Small-group size range
    q_min = 10;  q_max = nPop;   % Cluster size range
    CF_init = 0.8;               % Initial cooperation factor

    % ---- ERBMO adaptation weights (Table 4 of the paper) ----
    alpha_cf = 0.95;
    beta_cf  = 0.05;
    gamma_cf = 0.01;

    % ---- Periodic pattern search parameters (Eq. 12) ----
    delta_step_init = 0.5;   % Initial step size
    b_acc           = 1.2;   % Acceleration factor for pattern move
    a_step          = 0.5;   % Step-size reduction rate
    alpha_min       = 1e-3;  % Minimum step size

    % ---- Evolutionary combinatorial mutation (Eq. 13-15) ----
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
    [BestScore, idx] = min(Fit);
    BestPos = X(idx, :);

    CF              = CF_init;
    PrevBest        = BestScore;
    delta_local     = delta_step_init;
    ConvergenceCurve = zeros(1, MaxIter);

    % ===== Main loop =====
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % ===== Diversity-adaptive CF update (Eq. 7) =====
        improvement = abs(PrevBest - BestScore);
        diversity   = std(Fit);
        CF = alpha_cf * CF ...
           + beta_cf  * (1 / (1 + improvement^(1 + rand()))) ...
           + gamma_cf * diversity;
        CF = min(max(CF, 0.01), 1.5);

        do_pattern = (mod(t, 100) == 0);
        p_current  = min(p_mut_init + dp_mut * (t - 1), p_mut_max);

        % ===== RBMO core: food searching + prey attacking =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            Xnew = X(i, :);
            if rand() < 0.5
                % ----- Food searching phase -----
                if rand() < 0.5
                    % Small-group searching (Eq. 2)
                    p_size = randi([p_min, min(p_max, nPop)]);
                    idxs   = randperm(nPop, p_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xrs    = X(randi(nPop), :);
                    Xnew   = X(i, :) + (Xm - Xrs) * rand();
                else
                    % Cluster searching (Eq. 3)
                    q_size = randi([q_min, min(q_max, nPop)]);
                    idxs   = randperm(nPop, q_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xrs    = X(randi(nPop), :);
                    Xnew   = X(i, :) + (Xm - Xrs) * rand();
                end
            else
                % ----- Prey attacking phase -----
                if rand() < 0.5
                    % Small-group attacking (Eq. 4)
                    p_size = randi([p_min, min(p_max, nPop)]);
                    idxs   = randperm(nPop, p_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xnew   = BestPos + CF * (Xm - X(i, :)) * randn();
                else
                    % Cluster attacking (Eq. 5)
                    q_size = randi([q_min, min(q_max, nPop)]);
                    idxs   = randperm(nPop, q_size);
                    Xm     = mean(X(idxs, :), 1);
                    Xnew   = BestPos + CF * (Xm - X(i, :)) * randn();
                end
            end

            Xnew   = min(max(Xnew, lb), ub);
            NewFit = fobj(Xnew);  FEs = FEs + 1;

            % ----- Food caching (Eq. 6) -----
            if NewFit < Fit(i)
                X(i, :) = Xnew;  Fit(i) = NewFit;
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = Xnew;
                end
            end
        end

        % ===== Evolutionary probabilistic combinatorial mutation (Eq. 13-15) =====
        for i = 1:nPop
            if FEs >= MaxFEs, break; end
            if rand() < p_current
                if rand() < p_gauss
                    X_mut = X(i, :) + sigma_gauss * randn(1, dim);   % Eq. 14
                else
                    X_mut = lb + (ub - lb) .* rand(1, dim);           % Eq. 15
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

        % ===== Periodic pattern search (Eq. 8-12) =====
        if do_pattern && FEs < MaxFEs
            x1     = BestPos;   % Eq. 9: store pre-exploration position
            x_curr = BestPos;
            improved_explore = false;

            % ----- Step 1: Exploratory moves along each axis (Eq. 8) -----
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

            % ----- Step 2: Pattern move (Eq. 10-11) -----
            if improved_explore && FEs < MaxFEs
                x = x_curr;                                    % Eq. 10
                y_pattern = x + b_acc * (x - x1);              % Eq. 11
                y_pattern = min(max(y_pattern, lb), ub);
                f_pattern = fobj(y_pattern);  FEs = FEs + 1;
                if f_pattern < BestScore
                    BestPos = y_pattern;  BestScore = f_pattern;
                end
            end

            % ----- Step 3: Step-size adjustment (Eq. 12) -----
            if delta_local < alpha_min
                delta_local = delta_local * a_step;
            end
        end

        PrevBest = BestScore;
        ConvergenceCurve(t) = BestScore;
    end
end