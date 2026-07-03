function [BestScore, BestPos, ConvergenceCurve] = L_SHADE( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% L-SHADE: Linear Population Size Reduction - Success-History Adaptive DE
% Reference: Tanabe, R., & Fukunaga, A. S. (2014). Improving the search
%           performance of SHADE using linear population size reduction.
%           IEEE Congress on Evolutionary Computation (CEC).
%
% Input:
%   nPop      : Initial population size
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

    % L-SHADE parameters
    H = 6;                  % Historical memory size
    minPop = 4;             % Minimum population size
    p_best_rate = 0.11;     % p-best selection rate
    archive_rate = 2.6;     % Archive size factor

    % Function evaluation counter
    FEs = 0;
    currentPop = nPop;

    % Initialize population
    pop = lb + rand(currentPop, dim) .* (ub - lb);
    fitness = zeros(currentPop, 1);

    for i = 1:currentPop
        fitness(i) = fobj(pop(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    % Find initial best
    [BestScore, best_idx] = min(fitness);
    BestPos = pop(best_idx, :);

    % Initialize archive (empty at start)
    archive = [];

    % Initialize historical memories for CR and F
    M_CR = 0.5 * ones(H, 1);
    M_F = 0.5 * ones(H, 1);
    k = 1;  % Memory index

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Linear Population Size Reduction (Eq. 4)
        newPop = round(minPop + (nPop - minPop) * (MaxFEs - FEs) / MaxFEs);
        newPop = max(newPop, minPop);

        % Sort population for p-best selection
        [~, sorted_idx] = sort(fitness);
        p_best_size = max(1, round(p_best_rate * currentPop));

        % Storage for successful parameters
        S_CR = [];
        S_F = [];
        S_df = [];

        for i = 1:currentPop

            if FEs >= MaxFEs, break; end

            % Select memory index randomly
            r = randi(H);

            % Generate CR_i (normal distribution)
            CR_i = M_CR(r) + 0.1 * randn();
            CR_i = max(0, min(1, CR_i));

            % Generate F_i (Cauchy distribution)
            F_i = M_F(r) + 0.1 * tan(pi * (rand() - 0.5));
            while F_i <= 0
                F_i = M_F(r) + 0.1 * tan(pi * (rand() - 0.5));
            end
            F_i = min(F_i, 1);

            % Select p-best individual
            p_best_idx = sorted_idx(randi(p_best_size));

            % Select random index r1 (different from i)
            r1 = randi(currentPop);
            while r1 == i
                r1 = randi(currentPop);
            end

            % Select random index r2 from population U archive
            P_A = [pop; archive];
            r2 = randi(size(P_A, 1));
            while r2 == i || r2 == r1
                r2 = randi(size(P_A, 1));
            end

            % Mutation: current-to-pbest/1 (Eq. 1)
            mutant = pop(i, :) + ...
                     F_i * (pop(p_best_idx, :) - pop(i, :)) + ...
                     F_i * (pop(r1, :) - P_A(r2, :));

            % Boundary control
            mutant = max(min(mutant, ub), lb);

            % Binomial crossover
            trial = pop(i, :);
            j_rand = randi(dim);
            for j = 1:dim
                if j == j_rand || rand() <= CR_i
                    trial(j) = mutant(j);
                end
            end

            % Evaluate trial vector
            trial_fitness = fobj(trial);
            FEs = FEs + 1;

            % Selection
            if trial_fitness < fitness(i)
                % Store successful parameters
                S_CR = [S_CR; CR_i];
                S_F = [S_F; F_i];
                S_df = [S_df; (fitness(i) - trial_fitness)];

                % Update archive
                if size(archive, 1) < round(archive_rate * nPop)
                    archive = [archive; pop(i, :)];
                else
                    idx = randi(size(archive, 1));
                    archive(idx, :) = pop(i, :);
                end

                % Update population
                pop(i, :) = trial;
                fitness(i) = trial_fitness;

                % Update global best
                if trial_fitness < BestScore
                    BestScore = trial_fitness;
                    BestPos = trial;
                end
            end
        end

        % Update historical memory (Eq. 6-8)
        if ~isempty(S_CR)
            w = S_df / sum(S_df);
            M_CR(k) = sum(w .* S_CR);
            M_F(k) = sum(w .* (S_F.^2)) / sum(w .* S_F);
            k = mod(k, H) + 1;
        end

        % Linear Population Size Reduction
        if newPop < currentPop
            [~, idx] = sort(fitness);
            pop = pop(idx(1:newPop), :);
            fitness = fitness(idx(1:newPop));
            currentPop = newPop;
        end

        ConvergenceCurve(iter) = BestScore;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, BestScore, FEs);
        end
    end
end