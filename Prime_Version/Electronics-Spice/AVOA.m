function [BestScore, BestPos, ConvergenceCurve] = AVOA( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% AVOA: African Vultures Optimization Algorithm (2021)
% Reference: Abdollahzadeh, B., Gharehchopogh, F. S., & Mirjalili, S. (2021).
% African vultures optimization algorithm: A new nature-inspired
% metaheuristic algorithm for global optimization problems.
% Computers & Industrial Engineering, 158, 107408.

    % Convert bounds to row vectors
    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';
        ub = ub(:)';
    end

    FEs = 0;

    % Initialization (جایگزین تابع initialization)
    if numel(lb) == 1
        lb_vec = lb * ones(1, dim);
        ub_vec = ub * ones(1, dim);
    else
        lb_vec = lb(:)';
        ub_vec = ub(:)';
    end
    Positions = rand(nPop, dim) .* (ub_vec - lb_vec) + lb_vec;
    
    Fitness = zeros(1, nPop);
    for i = 1:nPop
        Fitness(i) = fobj(Positions(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    % Sort and identify best two vultures
    [~, sorted_idx] = sort(Fitness);
    BestPos = Positions(sorted_idx(1), :);
    SecondBest = Positions(sorted_idx(2), :);
    BestScore = Fitness(sorted_idx(1));

    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % Parameters (Eq. 2-4)
        R = -2 * rand * (1 - t/MaxIter) + 2 * rand + 1;  % → [-2, 2]
        F = (2 * rand + 1) * (1 - t/MaxIter)^(2*rand) + R;  % Eq. (3)

        % Pre-calculate probabilities (Eq. 5-6)
        p1 = 0.6 * (1 - t/MaxIter)^0.5;
        p2 = 0.4 * (1 - t/MaxIter)^0.5;
        p3 = 0.1 * (1 - t/MaxIter)^0.5;

        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % Choose target between Best and SecondBest (Eq. 7)
            L1 = 0.8; L2 = 0.2;
            if Fitness(sorted_idx(1)) / (Fitness(sorted_idx(2)) + eps) < 0.9
                R_target = BestPos;
            else
                if rand < L1
                    R_target = BestPos;
                else
                    R_target = SecondBest;
                end
            end

            if abs(F) >= 1  % Exploration phase (Eq. 8-9)
                if rand < p1
                    % Eq. (8)
                    random_idx = randi(nPop);
                    NewPos = R_target - abs(R_target - Positions(i, :)) .* F .* Levy(dim);
                else
                    % Eq. (9)
                    NewPos = R_target - F + rand(1, dim) .* ...
                        ((ub - lb) .* rand(1, dim) + lb);
                end

            else  % Exploitation phase (|F| < 1)
                if abs(F) >= 0.5  % Phase 1: Competition for food (Eq. 10-11)
                    if rand < p2
                        % Eq. (10): Rotating flight
                        s1 = R_target .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* cos(Positions(i, :));
                        s2 = R_target .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* sin(Positions(i, :));
                        NewPos = R_target - (s1 + s2);
                    else
                        % Eq. (11): Aggressive competition
                        NewPos = (R_target - abs(R_target - Positions(i, :))) .* F .* Levy(dim);
                    end
                else  % Phase 2: Aggressive competition (Eq. 12-13)
                    if rand < p3
                        % Eq. (12): Accumulation
                        s1 = BestPos .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* cos(Positions(i, :));
                        s2 = BestPos .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* sin(Positions(i, :));
                        NewPos = (R_target + s1 + s2) / 3;
                    else
                        % Eq. (13): Siege fight
                        rp = randperm(nPop, 1);
                        NewPos = R_target - abs(R_target - Positions(i, :)) .* F .* Levy(dim) + ...
                            (Positions(rp, :) - rand(1, dim) .* Positions(i, :));
                    end
                end
            end

            % Boundary control
            NewPos = min(max(NewPos, lb), ub);
            NewFit = fobj(NewPos);
            FEs = FEs + 1;

            % Greedy selection
            if NewFit < Fitness(i)
                Positions(i, :) = NewPos;
                Fitness(i) = NewFit;
                if NewFit < BestScore
                    BestScore = NewFit;
                    BestPos = NewPos;
                end
            end
        end

        % Re-sort and update second best
        [~, sorted_idx] = sort(Fitness);
        BestPos = Positions(sorted_idx(1), :);
        SecondBest = Positions(sorted_idx(2), :);
        BestScore = Fitness(sorted_idx(1));

        ConvergenceCurve(t) = BestScore;

        if mod(t, 50) == 0
            fprintf('AVOA Iter %4d | Best = %.4e | FEs = %d\n', t, BestScore, FEs);
        end
    end
end

% Levy flight function (Eq. 14)
function L = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / ...
        (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v).^(1 / beta);
    L = 0.01 * step;
end