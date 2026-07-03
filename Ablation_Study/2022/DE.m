function [Best_score, Best_pos, Convergence_curve] = DE( ...
    N, Max_iter, lb, ub, dim, fobj, MaxFEs)
% DE: Differential Evolution Algorithm
% Based on: Storn, R., & Price, K. (1997). Differential evolution – a simple 
%           and efficient heuristic for global optimization over continuous spaces.
%
% Input:
%   N         : Population size
%   Max_iter  : Maximum number of iterations
%   lb        : Lower bounds (1 x dim vector or scalar)
%   ub        : Upper bounds (1 x dim vector or scalar)
%   dim       : Dimension of the problem
%   fobj      : Objective function handle (minimization)
%   MaxFEs    : Maximum number of function evaluations
%
% Output:
%   Best_score     : Best fitness value found
%   Best_pos       : Best solution found
%   Convergence_curve : Best fitness value at each iteration

    % DE Parameters
    F = 0.8;   % Mutation factor
    CR = 0.9;  % Crossover probability

    % Initialize positions of agents
    Position = rand(N, dim) .* (ub - lb) + lb;
    Fitness = inf(N, 1);

    % Function evaluations counter
    FEs = 0;

    % Evaluate initial population
    for i = 1:N
        Fitness(i) = fobj(Position(i, :));
        FEs = FEs + 1;
    end

    % Initialize best solution
    [Best_score, best_idx] = min(Fitness);
    Best_pos = Position(best_idx, :);
    Convergence_curve = zeros(1, Max_iter);

    % Main optimization loop
    for iter = 1:Max_iter

        % Stop if MaxFEs reached
        if FEs >= MaxFEs
            Convergence_curve = Convergence_curve(1:iter-1);
            return;
        end

        for i = 1:N

            if FEs >= MaxFEs
                break;
            end

            % Mutation: Generate donor vector (DE/rand/1)
            idxs = randperm(N, 3);
            while any(idxs == i)
                idxs = randperm(N, 3);
            end
            donor = Position(idxs(1), :) + ...
                     F * (Position(idxs(2), :) - Position(idxs(3), :));

            % Crossover: Generate trial vector (binomial crossover)
            trial = Position(i, :);
            rand_idx = randi(dim);
            for j = 1:dim
                if rand < CR || j == rand_idx
                    trial(j) = donor(j);
                end
            end

            % Boundary control
            trial = max(min(trial, ub), lb);

            % Selection
            trial_fitness = fobj(trial);
            FEs = FEs + 1;

            if trial_fitness < Fitness(i)
                Position(i, :) = trial;
                Fitness(i) = trial_fitness;
            end

            % Update global best
            if trial_fitness < Best_score
                Best_score = trial_fitness;
                Best_pos = trial;
            end
        end

        % Store convergence information
        Convergence_curve(iter) = Best_score;

        % Display progress
        if mod(iter, ceil(Max_iter/10)) == 0
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                    iter, Best_score, FEs);
        end
    end
end
