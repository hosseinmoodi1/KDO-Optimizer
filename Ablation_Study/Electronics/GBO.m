function [BestScore, BestPos, ConvergenceCurve] = GBO( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% GBO: Gradient-Based Optimizer
% Reference: Ahmadianfar, I., Bozorg-Haddad, O., & Chu, X. (2020).
%           Gradient-based optimizer: A new metaheuristic optimization
%           algorithm. Information Sciences.
%
% Input:
%   nPop      : Population size
%   MaxIter   : Maximum number of iterations
%   lb        : Lower bounds (scalar or 1 x dim vector)
%   ub        : Upper bounds (scalar or 1 x dim vector)
%   dim       : Dimension of the problem
%   fobj      : Objective function handle (minimization)
%   MaxFEs    : Maximum number of function evaluations
%
% Output:
%   BestScore         : Best fitness value found
%   BestPos           : Best solution found
%   ConvergenceCurve  : Best fitness value at each iteration

    % Convert bounds to row vectors
    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';
        ub = ub(:)';
    end

    % GBO parameters
    pr = 0.5;               % Probability for LEO
    beta_min = 0.2;         % Minimum beta (Eq. 6)
    beta_max = 1.2;         % Maximum beta (Eq. 6)

    % Function evaluation counter
    FEs = 0;

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;
    Fit = zeros(nPop, 1);

    for i = 1:nPop
        Fit(i) = fobj(X(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    % Find initial best
    [BestScore, best_idx] = min(Fit);
    BestPos = X(best_idx, :);

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Adaptive beta parameter (Eq. 6)
        beta = beta_min + (beta_max - beta_min) * (1 - (iter/MaxIter)^3)^2;

        % Sort for best and worst identification
        [~, idx_sort] = sort(Fit);
        X_best = X(idx_sort(1), :);
        X_worst = X(idx_sort(end), :);
        X_mean = mean(X, 1);

        for i = 1:nPop

            if FEs >= MaxFEs, break; end

            % Select four random indices different from i
            idxs = randperm(nPop, 4);
            while any(idxs == i)
                idxs = randperm(nPop, 4);
            end
            r1 = idxs(1); r2 = idxs(2); r3 = idxs(3); r4 = idxs(4);

            % Gradient Search Rule (GSR) - Eq. 4 & 5
            delta = 2 * rand() * ((X(r1,:) + X(r2,:) + X(r3,:) + X(r4,:)) / 4 - X(i,:));
            step = (X_best - X(r1,:) + delta) / 2;
            DeltaX = beta * rand(1, dim) .* step;

            alpha = beta * sin(3*pi/2 + sin(3*pi/2*beta));
            rho1 = 2 * rand() * alpha - alpha;
            rho2 = 2 * rand() * alpha - alpha;

            GSR = randn() * rho1 * (2 * DeltaX .* X(i,:)) ./ (X_worst - X_best + eps);

            X1 = X(i,:) - GSR + rand() * rho2 * (X_best - X(i,:));

            % Random operator (Eq. 8)
            if rand() < 0.5
                rA = randi(nPop);
                rB = randi(nPop);
                while rA == rB
                    rB = randi(nPop);
                end
                X1 = X1 + 1e-3 * (X(rA,:) - X(rB,:));
            end

            % Calculate X2 and X3 (Eq. 10-11)
            X2 = X_best - rho1 * (2 * DeltaX .* X(i,:)) ./ (X_worst - X_best + eps) + ...
                 rand() * rho2 * (X(r1,:) - X(r2,:));
            X3 = X(i,:) - rho1 * (X2 - X1);

            % Combine to get new solution (Eq. 12)
            if rand() < 0.5
                X_new = X1 + rand() * (X2 - X1);
            else
                X_new = X1 + rand() * (X3 - X1);
            end

            % Local Escaping Operator (LEO) - Eq. 13-17
            if rand() < pr
                u1 = rand(1, dim);
                u2 = rand(1, dim);

                if rand() < 0.5
                    us = u1 .* abs(tan(rand(1, dim) * pi/2));
                else
                    us = u1 .* (2 * rand(1, dim) - 1);
                end

                k = randi(nPop);
                r = randi(nPop);

                L1 = exp(-u2 .* abs(tan(rand(1, dim) * pi/2)));
                L2 = 1 - L1;

                X_leo = X(k,:) + us .* (X(r,:) - L1 .* X_best - L2 .* X_mean) + ...
                        us .* (L1 .* X_worst - L2 .* X(r,:));

                if rand() < 0.5
                    X_new = X_new + rand() * (X_leo - X_new);
                else
                    X_new = X_leo;
                end
            end

            % Boundary control
            X_new = max(min(X_new, ub), lb);

            % Evaluate new solution
            NewFit = fobj(X_new);
            FEs = FEs + 1;

            % Greedy selection
            if NewFit < Fit(i)
                X(i,:) = X_new;
                Fit(i) = NewFit;

                if NewFit < BestScore
                    BestScore = NewFit;
                    BestPos = X_new;
                end
            end
        end

        ConvergenceCurve(iter) = BestScore;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, BestScore, FEs);
        end
    end
end