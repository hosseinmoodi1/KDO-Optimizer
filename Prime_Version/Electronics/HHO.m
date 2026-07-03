function [BestScore, BestPos, ConvergenceCurve] = HHO( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% HHO: Harris Hawks Optimization
% Reference: Heidari, A. A., Mirjalili, S., Faris, H., et al. (2019).
%           Harris hawks optimization: Algorithm and applications.
%           Future Generation Computer Systems.
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

    % Find initial rabbit (best solution)
    [BestScore, rabbit_idx] = min(fitness);
    rabbit_pos = X(rabbit_idx, :);

    % Convergence curve
    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for iter = 1:MaxIter

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Escaping energy (Eq. 3)
        E0 = 2 * rand() - 1;                    % Initial energy [-1, 1]
        E = 2 * E0 * (1 - iter / MaxIter);      % Decreasing energy

        % Jump strength (Eq. 5)
        J = 2 * (1 - rand());

        for i = 1:nPop

            % Store old position for possible rollback
            X_old = X(i, :);
            fit_old = fitness(i);

            % ---------- Exploration Phase (|E| >= 1) ----------
            if abs(E) >= 1
                % Random hawk index
                rand_idx = randi(nPop);
                while rand_idx == i
                    rand_idx = randi(nPop);
                end

                if rand() >= 0.5
                    % Eq. 1: Perching based on random location
                    X(i, :) = X(rand_idx, :) - rand() * ...
                              abs(X(rand_idx, :) - 2 * rand() * X(i, :));
                else
                    % Eq. 2: Perching based on average position
                    X_mean = mean(X, 1);
                    X(i, :) = (rabbit_pos - X_mean) - ...
                              rand() * (lb + rand() .* (ub - lb));
                end

            % ---------- Exploitation Phase (|E| < 1) ----------
            else
                r = rand();

                % Soft besiege (|E| >= 0.5, r >= 0.5) - Eq. 6
                if abs(E) >= 0.5 && r >= 0.5
                    X(i, :) = rabbit_pos - E * ...
                              abs(J * rabbit_pos - X(i, :));

                % Hard besiege (|E| < 0.5, r >= 0.5) - Eq. 7
                elseif abs(E) < 0.5 && r >= 0.5
                    delta_X = rabbit_pos - X(i, :);
                    X(i, :) = rabbit_pos - E * abs(delta_X);

                % Soft besiege with rapid dives (|E| >= 0.5, r < 0.5)
                elseif abs(E) >= 0.5 && r < 0.5
                    % First dive (Eq. 11)
                    Y = rabbit_pos - E * abs(J * rabbit_pos - X(i, :));
                    Y = max(min(Y, ub), lb);
                    fY = fobj(Y);
                    FEs = FEs + 1;

                    if fY < fit_old
                        X(i, :) = Y;
                        fitness(i) = fY;
                    else
                        % Second dive with Levy flight (Eq. 12)
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

                % Hard besiege with rapid dives (|E| < 0.5, r < 0.5)
                else
                    X_mean = mean(X, 1);
                    % First dive (Eq. 15)
                    Y = rabbit_pos - E * abs(J * rabbit_pos - X_mean);
                    Y = max(min(Y, ub), lb);
                    fY = fobj(Y);
                    FEs = FEs + 1;

                    if fY < fit_old
                        X(i, :) = Y;
                        fitness(i) = fY;
                    else
                        % Second dive with Levy flight (Eq. 16)
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

            % Evaluate if not already evaluated in rapid dives
            if ~isequal(X(i, :), X_old) && ...
               (abs(E) >= 1 || (abs(E) < 0.5 && r >= 0.5) || ...
                (abs(E) >= 0.5 && r >= 0.5))
                new_fit = fobj(X(i, :));
                FEs = FEs + 1;
                if new_fit < fitness(i)
                    fitness(i) = new_fit;
                else
                    X(i, :) = X_old;
                    fitness(i) = fit_old;
                end
            end

            % Update rabbit (global best)
            if fitness(i) < BestScore
                BestScore = fitness(i);
                rabbit_pos = X(i, :);
            end

            if FEs >= MaxFEs, break; end
        end

        ConvergenceCurve(iter) = BestScore;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, BestScore, FEs);
        end
    end

    BestPos = rabbit_pos;
end

% Levy flight function for HHO
function L = levy_flight(d)
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
             (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    L = 0.01 * u ./ (abs(v).^(1/beta));
end