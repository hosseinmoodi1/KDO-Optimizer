function [BestScore, BestPos, ConvergenceCurve] = WOA( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% WOA: Whale Optimization Algorithm
% Reference: Mirjalili, S., & Lewis, A. (2016). The Whale Optimization
%           Algorithm. Advances in Engineering Software.
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

    % Function evaluation counter
    FEs = 0;

    % Initialize population
    X = rand(nPop, dim) .* (ub - lb) + lb;
    fitness = zeros(nPop, 1);

    for i = 1:nPop
        fitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    % Find initial best (prey position)
    [BestScore, bestIdx] = min(fitness);
    BestPos = X(bestIdx, :);

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Linearly decreasing a from 2 to 0 (Eq. 2 in WOA paper)
        a = 2 - iter * (2 / MaxIter);

        for i = 1:nPop

            % Update parameters
            r1 = rand(); r2 = rand();
            A = 2 * a * r1 - a;     % Eq. 2.3
            C = 2 * r2;             % Eq. 2.4
            b = 1;                  % Spiral constant
            l = -1 + 2 * rand();    % Random in [-1, 1]
            p = rand();             % Probability

            if p < 0.5
                % Shrinking encircling mechanism
                if abs(A) < 1
                    % Exploitation: update towards prey (Eq. 2.1)
                    D = abs(C * BestPos - X(i, :));
                    X(i, :) = BestPos - A * D;
                else
                    % Exploration: search randomly (Eq. 2.8)
                    rand_idx = randi(nPop);
                    X_rand = X(rand_idx, :);
                    D = abs(C * X_rand - X(i, :));
                    X(i, :) = X_rand - A * D;
                end
            else
                % Spiral updating position (Eq. 2.7)
                D_prime = abs(BestPos - X(i, :));
                X(i, :) = D_prime * exp(b * l) * cos(2 * pi * l) + BestPos;
            end

            % Boundary control
            X(i, :) = max(min(X(i, :), ub), lb);

            % Evaluate new position
            new_fitness = fobj(X(i, :));
            FEs = FEs + 1;

            % Greedy selection
            if new_fitness < fitness(i)
                X(i, :) = X(i, :);
                fitness(i) = new_fitness;

                if new_fitness < BestScore
                    BestScore = new_fitness;
                    BestPos = X(i, :);
                end
            end

            if FEs >= MaxFEs, break; end
        end

        ConvergenceCurve(iter) = BestScore;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, BestScore, FEs);
        end
    end
end