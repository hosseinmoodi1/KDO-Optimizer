function [Best_score, Best_pos, Convergence_curve] = GWO( ...
    N, Max_iter, lb, ub, dim, fobj, MaxFEs)
% GWO: Grey Wolf Optimizer
% Based on: Mirjalili et al., "Grey Wolf Optimizer," 
%           Advances in Engineering Software, vol. 69, pp. 46-61, 2014.
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

    % Initialize alpha, beta, delta wolves
    Alpha_pos = zeros(1, dim);
    Alpha_score = inf;

    Beta_pos = zeros(1, dim);
    Beta_score = inf;

    Delta_pos = zeros(1, dim);
    Delta_score = inf;

    % Initialize population
    Positions = lb + rand(N, dim) .* (ub - lb);

    % Initialize convergence curve
    Convergence_curve = zeros(1, Max_iter);

    % Function evaluations counter
    FEs = 0;

    % Main optimization loop
    for t = 1:Max_iter

        % Stop if MaxFEs reached
        if FEs >= MaxFEs
            Convergence_curve = Convergence_curve(1:t-1);
            break;
        end
        
        % Evaluate all wolves
        for i = 1:N

            if FEs >= MaxFEs
                break;
            end

            % Boundary check
            Positions(i, :) = min(max(Positions(i, :), lb), ub);
            
            % Calculate fitness
            fitness = fobj(Positions(i, :));
            FEs = FEs + 1;
            
            % Update Alpha, Beta, Delta wolves
            if fitness < Alpha_score
                Delta_score = Beta_score;
                Delta_pos = Beta_pos;
                
                Beta_score = Alpha_score;
                Beta_pos = Alpha_pos;
                
                Alpha_score = fitness;
                Alpha_pos = Positions(i, :);
                
            elseif fitness < Beta_score
                Delta_score = Beta_score;
                Delta_pos = Beta_pos;
                
                Beta_score = fitness;
                Beta_pos = Positions(i, :);
                
            elseif fitness < Delta_score
                Delta_score = fitness;
                Delta_pos = Positions(i, :);
            end
        end
        
        % a decreases linearly from 2 to 0
        a = 2 - t * (2 / Max_iter);
        
        % Update positions of all wolves
        for i = 1:N
            for j = 1:dim
                % Update positions based on alpha wolf
                r1 = rand(); r2 = rand();
                A1 = 2 * a * r1 - a;
                C1 = 2 * r2;
                D_alpha = abs(C1 * Alpha_pos(j) - Positions(i, j));
                X1 = Alpha_pos(j) - A1 * D_alpha;
                
                % Update positions based on beta wolf
                r1 = rand(); r2 = rand();
                A2 = 2 * a * r1 - a;
                C2 = 2 * r2;
                D_beta = abs(C2 * Beta_pos(j) - Positions(i, j));
                X2 = Beta_pos(j) - A2 * D_beta;
                
                % Update positions based on delta wolf
                r1 = rand(); r2 = rand();
                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;
                D_delta = abs(C3 * Delta_pos(j) - Positions(i, j));
                X3 = Delta_pos(j) - A3 * D_delta;
                
                % Update wolf position
                Positions(i, j) = (X1 + X2 + X3) / 3;
            end
        end
        
        % Store convergence
        Convergence_curve(t) = Alpha_score;
        
        % Display progress
        if mod(t, ceil(Max_iter/10)) == 0 || t == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                    t, Alpha_score, FEs);
        end
    end
    
    % Final results
    Best_score = Alpha_score;
    Best_pos = Alpha_pos;
end
