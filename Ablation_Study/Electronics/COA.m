function [BestScore, BestPos, ConvergenceCurve] = COA( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% COA: Crayfish Optimization Algorithm (2023)
% Reference: Jia, H., Rao, H., Wen, C., & Mirjalili, S. (2023).
% Crayfish optimization algorithm. Artificial Intelligence Review.

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

    [BestScore, best_idx] = min(Fitness);
    BestPos = Positions(best_idx, :);

    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % Temperature: Eq. (3) in paper → Temp = rand × 15 + 20
        Temp = rand * 15 + 20;

        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % Adaptive step k (Eq. 4)
            k = 2 * rand - 1;

            if Temp > 30  % Summer resort → Exploration
                % Random cave location (Eq. 5-6)
                X_cave = (BestPos + Positions(i, :)) / 2;

                if rand < 0.5
                    % Enter cave → Eq. (7)
                    NewPos = Positions(i, :) + ...
                        (2 - t/MaxIter) * rand(1, dim) .* (X_cave - Positions(i, :));
                else
                    % Competition for cave → Eq. (8)
                    Idx = randperm(nPop, 2);
                    NewPos = Positions(i, :) - Positions(Idx(1), :) + X_cave;
                end
            else  % Temp ≤ 30 → Exploitation
                % Foraging amount (Eq. 9-10)
                p = 0.2 + 0.3 * (1 - t/MaxIter);
                Q = 2 * rand(1, dim) - 1;

                if rand < p
                    % Foraging behavior → Eq. (11)
                    NewPos = Positions(i, :) + k * (BestPos - Positions(i, :)) .* Q;
                else
                    % Fighting behavior → Eq. (12)
                    z = round(rand(1, dim));
                    NewPos = Positions(i, :) + k * (BestPos - z .* Positions(i, :));
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

        ConvergenceCurve(t) = BestScore;

        if mod(t, 50) == 0
            fprintf('COA Iter %4d | Best = %.4e | FEs = %d\n', t, BestScore, FEs);
        end
    end
end