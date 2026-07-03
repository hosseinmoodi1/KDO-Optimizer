function [BestScore, BestPos, ConvergenceCurve] = GWO( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% GWO: Grey Wolf Optimizer
% Reference: Mirjalili, S., Mirjalili, S. M., & Lewis, A. (2014).
%           Grey Wolf Optimizer. Advances in Engineering Software.
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

    % Initialize alpha, beta, delta wolves
    [~, idx] = sort(fitness);
    alpha_pos = X(idx(1), :);
    alpha_score = fitness(idx(1));
    beta_pos = X(idx(2), :);
    beta_score = fitness(idx(2));
    delta_pos = X(idx(3), :);
    delta_score = fitness(idx(3));

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Linearly decreasing a from 2 to 0 (Eq. 3)
        a = 2 - iter * (2 / MaxIter);

        for i = 1:nPop

            for j = 1:dim

                % Update around alpha (Eq. 5-7)
                r1 = rand(); r2 = rand();
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = abs(C1 * alpha_pos(j) - X(i, j));
                X1 = alpha_pos(j) - A1 * D_alpha;

                % Update around beta
                r1 = rand(); r2 = rand();
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * beta_pos(j) - X(i, j));
                X2 = beta_pos(j) - A2 * D_beta;

                % Update around delta
                r1 = rand(); r2 = rand();
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_delta = abs(C3 * delta_pos(j) - X(i, j));
                X3 = delta_pos(j) - A3 * D_delta;

                % Position update (Eq. 11)
                X(i, j) = (X1 + X2 + X3) / 3;
            end

            % Boundary control
            X(i, :) = max(min(X(i, :), ub), lb);

            % Evaluate new position
            new_fitness = fobj(X(i, :));
            FEs = FEs + 1;

            % Update alpha, beta, delta
            if new_fitness < alpha_score
                delta_score = beta_score;
                delta_pos = beta_pos;
                beta_score = alpha_score;
                beta_pos = alpha_pos;
                alpha_score = new_fitness;
                alpha_pos = X(i, :);
            elseif new_fitness < beta_score
                delta_score = beta_score;
                delta_pos = beta_pos;
                beta_score = new_fitness;
                beta_pos = X(i, :);
            elseif new_fitness < delta_score
                delta_score = new_fitness;
                delta_pos = X(i, :);
            end

            if FEs >= MaxFEs, break; end
        end

        ConvergenceCurve(iter) = alpha_score;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, alpha_score, FEs);
        end
    end

    BestScore = alpha_score;
    BestPos = alpha_pos;
end