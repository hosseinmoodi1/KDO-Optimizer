function [BestScore, BestPos, ConvergenceCurve] = DE( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% DE: Differential Evolution
% Reference: Storn, R., & Price, K. (1997). Differential evolution -
%           a simple and efficient heuristic for global optimization
%           over continuous spaces. Journal of Global Optimization.
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

    % DE parameters (standard: DE/rand/1/bin)
    F = 0.8;        % Mutation factor (range [0,2])
    CR = 0.9;       % Crossover probability (range [0,1])

    % Function evaluation counter
    FEs = 0;

    % Initialize population
    pop = rand(nPop, dim) .* (ub - lb) + lb;
    fit = zeros(nPop, 1);

    for i = 1:nPop
        fit(i) = fobj(pop(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    % Find initial best
    [BestScore, bestIdx] = min(fit);
    BestPos = pop(bestIdx, :);

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        for i = 1:nPop

            if FEs >= MaxFEs
                break;
            end

            % Mutation: DE/rand/1
            % Select three distinct random indices different from i
            idx = randperm(nPop, 3);
            while any(idx == i)
                idx = randperm(nPop, 3);
            end
            r1 = idx(1); r2 = idx(2); r3 = idx(3);

            mutant = pop(r1, :) + F * (pop(r2, :) - pop(r3, :));

            % Boundary control after mutation
            mutant = max(min(mutant, ub), lb);

            % Crossover: binomial
            trial = pop(i, :);
            j_rand = randi(dim);
            for j = 1:dim
                if rand() < CR || j == j_rand
                    trial(j) = mutant(j);
                end
            end

            % Evaluate trial vector
            trialFit = fobj(trial);
            FEs = FEs + 1;

            % Selection (greedy)
            if trialFit < fit(i)
                pop(i, :) = trial;
                fit(i) = trialFit;

                if trialFit < BestScore
                    BestScore = trialFit;
                    BestPos = trial;
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