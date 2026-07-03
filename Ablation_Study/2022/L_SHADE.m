function [Best_score, Best_pos, Convergence_curve] = L_SHADE( ...
    N, Max_iter, lb, ub, dim, fobj, MaxFEs)
% L_SHADE: Linear Population Size Reduction - Success-History Adaptive DE
% Reference: Tanabe, R., & Fukunaga, A. S. (2014).
%           Improving the search performance of SHADE using linear population
%           size reduction. IEEE CEC 2014.
%
% Input:
%   N         : Initial population size
%   Max_iter  : Maximum number of iterations
%   lb        : Lower bounds (scalar or 1 x dim vector)
%   ub        : Upper bounds (scalar or 1 x dim vector)
%   dim       : Dimension of the problem
%   fobj      : Objective function handle (minimization)
%   MaxFEs    : Maximum number of function evaluations
%
% Output:
%   Best_score         : Best fitness value found
%   Best_pos           : Best solution found
%   Convergence_curve  : Best fitness value at each iteration

    % L-SHADE parameters
    H = 6;                  % Size of historical memory
    minN = 4;               % Minimum population size
    p = 0.11;               % p-best selection rate
    archive_rate = 2.6;     % Archive size factor (2.6 * N)
    
    % Convert bounds to row vectors
    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';
        ub = ub(:)';
    end
    
    % Initialize population
    currentN = N;
    pop = lb + rand(currentN, dim) .* (ub - lb);
    fitness = zeros(currentN, 1);
    
    % Function evaluation counter
    FEs = 0;
    
    % Evaluate initial population
    for i = 1:currentN
        fitness(i) = fobj(pop(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs
            break;
        end
    end
    
    % Find initial global best
    [Best_score, best_idx] = min(fitness);
    Best_pos = pop(best_idx, :);
    
    % Initialize archive (empty at start)
    archive = [];
    
    % Initialize historical memory for CR and F
    memory_CR = 0.5 * ones(H, 1);
    memory_F = 0.5 * ones(H, 1);
    k = 1;  % Memory index
    
    % Convergence curve
    Convergence_curve = zeros(1, Max_iter);
    
    % Main loop
    for iter = 1:Max_iter
        
        if FEs >= MaxFEs
            Convergence_curve = Convergence_curve(1:iter-1);
            break;
        end
        
        % Linear Population Size Reduction (Eq. 4 in L-SHADE paper)
        % Based on remaining function evaluations, not iterations
        newN = round(minN + (N - minN) * (MaxFEs - FEs) / MaxFEs);
        newN = max(newN, minN);
        
        % Sort population by fitness (for p-best selection)
        [~, sorted_idx] = sort(fitness);
        pbest_size = max(1, round(p * currentN));
        
        % Store successful parameters and their fitness improvements
        S_CR = [];
        S_F = [];
        S_df = [];
        
        for i = 1:currentN
            
            if FEs >= MaxFEs
                break;
            end
            
            % Select memory index randomly
            r = randi(H);
            
            % Generate CR_i (normal distribution with mean memory_CR(r), std 0.1)
            CR_i = memory_CR(r) + 0.1 * randn();
            CR_i = max(0, min(1, CR_i));
            
            % Generate F_i (Cauchy distribution with location memory_F(r), scale 0.1)
            F_i = memory_F(r) + 0.1 * tan(pi * (rand() - 0.5));
            while F_i <= 0
                F_i = memory_F(r) + 0.1 * tan(pi * (rand() - 0.5));
            end
            F_i = min(F_i, 1);
            
            % Select p-best individual
            pbest_idx = sorted_idx(randi(pbest_size));
            
            % Select random index r1 (different from i)
            r1 = randi(currentN);
            while r1 == i
                r1 = randi(currentN);
            end
            
            % Select random index r2 from population U archive
            P_A = [pop; archive];
            r2 = randi(size(P_A, 1));
            while r2 == i || r2 == r1
                r2 = randi(size(P_A, 1));
            end
            
            % Mutation: current-to-pbest/1 (Eq. 1 in L-SHADE paper)
            mutant = pop(i, :) ...
                + F_i * (pop(pbest_idx, :) - pop(i, :)) ...
                + F_i * (pop(r1, :) - P_A(r2, :));
            
            % Boundary control after mutation
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
                
                % Archive management
                if size(archive, 1) < round(archive_rate * N)
                    archive = [archive; pop(i, :)];
                else
                    idx = randi(size(archive, 1));
                    archive(idx, :) = pop(i, :);
                end
                
                % Update population
                pop(i, :) = trial;
                fitness(i) = trial_fitness;
                
                % Update global best
                if trial_fitness < Best_score
                    Best_score = trial_fitness;
                    Best_pos = trial;
                end
            end
        end
        
        % Update historical memory (Eq. 6, 7, 8 in L-SHADE paper)
        if ~isempty(S_CR)
            % Calculate weights based on fitness improvement
            w = S_df / sum(S_df);
            
            % Update memory for CR (weighted arithmetic mean)
            memory_CR(k) = sum(w .* S_CR);
            
            % Update memory for F (weighted Lehmer mean)
            memory_F(k) = sum(w .* (S_F.^2)) / sum(w .* S_F);
            
            % Move to next memory cell
            k = mod(k, H) + 1;
        end
        
        % Linear Population Size Reduction (execute after parameter update)
        if newN < currentN
            [~, idx] = sort(fitness);
            pop = pop(idx(1:newN), :);
            fitness = fitness(idx(1:newN));
            currentN = newN;
        end
        
        % Store convergence
        Convergence_curve(iter) = Best_score;
        
        % Display progress
        if mod(iter, ceil(Max_iter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', iter, Best_score, FEs);
        end
    end
end