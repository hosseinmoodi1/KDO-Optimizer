function [Best_score, Best_pos, Convergence_curve] = HHO( ...
    N, Max_iter, lb, ub, dim, fobj, MaxFEs)
% HHO: Harris Hawks Optimization
% Reference: Heidari, A. A., Mirjalili, S., Faris, H., et al. (2019).
%           Harris hawks optimization: Algorithm and applications.
%           Future Generation Computer Systems, 97, 849-872.
%
% Input:
%   N         : Population size
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

    % Convert bounds to row vectors
    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';
        ub = ub(:)';
    end

    % Initialize population
    X = lb + rand(N, dim) .* (ub - lb);
    fitness = zeros(N, 1);

    % Function evaluation counter
    FEs = 0;

    % Evaluate initial population
    for i = 1:N
        fitness(i) = fobj(X(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs
            break;
        end
    end

    % Initialize rabbit (best solution)
    [rabbit_fitness, rabbit_idx] = min(fitness);
    rabbit_location = X(rabbit_idx, :);

    % Convergence curve
    Convergence_curve = zeros(1, Max_iter);

    % Main loop
    for t = 1:Max_iter

        if FEs >= MaxFEs
            Convergence_curve = Convergence_curve(1:t-1);
            break;
        end

        % Escaping energy (Eq. 3)
        E0 = 2 * rand() - 1;                    % Initial energy [-1, 1]
        E = 2 * E0 * (1 - t / Max_iter);        % Decreasing energy
        
        % Jump strength (Eq. 5)
        J = 2 * (1 - rand());

        for i = 1:N

            if FEs >= MaxFEs
                break;
            end

            X_old = X(i, :);
            fit_old = fitness(i);

            % ---------- Exploration Phase (|E| >= 1) ----------
            if abs(E) >= 1
                % Random hawk index
                idx_rand = randi(N);
                while idx_rand == i
                    idx_rand = randi(N);
                end

                if rand() >= 0.5
                    % Eq. 1: Perching based on random location
                    X(i, :) = X(idx_rand, :) - rand() * ...
                              abs(X(idx_rand, :) - 2 * rand() * X(i, :));
                else
                    % Eq. 2: Perching based on average position
                    X_mean = mean(X, 1);
                    X(i, :) = (rabbit_location - X_mean) - ...
                              rand() * (lb + rand() .* (ub - lb));
                end

            % ---------- Exploitation Phase (|E| < 1) ----------
            else
                r = rand();

                % Case 1: Soft besiege (|E| >= 0.5, r >= 0.5)
                if abs(E) >= 0.5 && r >= 0.5
                    X(i, :) = rabbit_location - E * ...
                              abs(J * rabbit_location - X(i, :));

                % Case 2: Hard besiege (|E| < 0.5, r >= 0.5)
                elseif abs(E) < 0.5 && r >= 0.5
                    delta_X = rabbit_location - X(i, :);
                    X(i, :) = rabbit_location - E * abs(delta_X);

                % Case 3: Soft besiege with progressive rapid dives (|E| >= 0.5, r < 0.5)
                elseif abs(E) >= 0.5 && r < 0.5
                    % First dive (Eq. 11)
                    Y = rabbit_location - E * abs(J * rabbit_location - X(i, :));
                    Y = max(min(Y, ub), lb);
                    fY = fobj(Y);
                    FEs = FEs + 1;
                    
                    if fY < fit_old
                        X(i, :) = Y;
                        fitness(i) = fY;
                    else
                        % Second dive (Eq. 12)
                        S = levy_flight(dim);
                        Z = Y + S .* rand(1, dim);
                        Z = max(min(Z, ub), lb);
                        fZ = fobj(Z);
                        FEs = FEs + 1;
                        if fZ < fit_old
                            X(i, :) = Z;
                            fitness(i) = fZ;
                        else
                            X(i, :) = X_old;
                            fitness(i) = fit_old;
                        end
                    end

                % Case 4: Hard besiege with progressive rapid dives (|E| < 0.5, r < 0.5)
                else
                    X_mean = mean(X, 1);
                    
                    % First dive (Eq. 15)
                    Y = rabbit_location - E * abs(J * rabbit_location - X_mean);
                    Y = max(min(Y, ub), lb);
                    fY = fobj(Y);
                    FEs = FEs + 1;
                    
                    if fY < fit_old
                        X(i, :) = Y;
                        fitness(i) = fY;
                    else
                        % Second dive (Eq. 16)
                        S = levy_flight(dim);
                        Z = Y + S .* rand(1, dim);
                        Z = max(min(Z, ub), lb);
                        fZ = fobj(Z);
                        FEs = FEs + 1;
                        if fZ < fit_old
                            X(i, :) = Z;
                            fitness(i) = fZ;
                        else
                            X(i, :) = X_old;
                            fitness(i) = fit_old;
                        end
                    end
                end
            end

            % Final boundary control
            X(i, :) = max(min(X(i, :), ub), lb);
            
            % Evaluate new position (only if not already evaluated in rapid dives)
            % Check if position changed and not already evaluated
            if ~isequal(X(i, :), X_old) && (abs(E) >= 1 || (abs(E) < 0.5 && r >= 0.5) || (abs(E) >= 0.5 && r >= 0.5))
                f_new = fobj(X(i, :));
                FEs = FEs + 1;
                
                if f_new < fitness(i)
                    fitness(i) = f_new;
                    if f_new < rabbit_fitness
                        rabbit_fitness = f_new;
                        rabbit_location = X(i, :);
                    end
                else
                    X(i, :) = X_old;
                    fitness(i) = fit_old;
                end
            else
                % Update rabbit if improved
                if fitness(i) < rabbit_fitness
                    rabbit_fitness = fitness(i);
                    rabbit_location = X(i, :);
                end
            end
        end

        % Store convergence
        Convergence_curve(t) = rabbit_fitness;

        % Display progress
        if mod(t, ceil(Max_iter/10)) == 0 || t == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', t, rabbit_fitness, FEs);
        end
    end

    Best_score = rabbit_fitness;
    Best_pos = rabbit_location;
end

function L = levy_flight(d)
% Levy flight function for HHO
% Generates a random step vector with Levy distribution
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
             (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    L = 0.01 * u ./ (abs(v).^(1/beta));
end