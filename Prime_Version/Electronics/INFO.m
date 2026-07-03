function [BestScore, BestPos, ConvergenceCurve] = INFO( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% INFO: Information-based Optimization
% Source: Ahmadianfar, I., Bozorgi-Amiri, A., & Zare, M. (2022).
% Information-based optimization: A novel metaheuristic algorithm.
%
% Modified for Electronics Benchmark (Hard Problems)
% - Added MaxFEs support
% - Adaptive alpha and delta parameters
% - Improved local search for electronics problems

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
        Fitness(i) = fobj(Positions(i,:));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end

    [BestScore, best_idx] = min(Fitness);
    BestPos = Positions(best_idx,:);

    ConvergenceCurve = zeros(1, MaxIter);

    % Main loop
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end

        % Adaptive parameters for electronics problems
        alpha = 0.5 + 0.5 * (1 - t/MaxIter);  % Gradually decreases
        delta = 2 * (1 - t/MaxIter);          % Decreases over time

        % Sort solutions
        [~, sorted_indices] = sort(Fitness);
        BestSol = Positions(sorted_indices(1),:);
        BetterSol = Positions(sorted_indices(2),:);
        WorstSol = Positions(sorted_indices(end),:);

        % Updating Rule Stage
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % Select three random agents
            Idx = randperm(nPop);
            Idx(Idx==i) = [];
            a = Idx(1); b = Idx(2); c = Idx(3);

            % Calculate weighted mean
            if Fitness(a) < Fitness(b)
                u1 = Positions(a,:);
                u2 = Positions(b,:);
            else
                u1 = Positions(b,:);
                u2 = Positions(a,:);
            end

            Mean_rule = (Positions(a,:) + Positions(b,:) + Positions(c,:)) / 3;
            z = rand(1, dim);

            if rand > 0.5
                NewPos1 = Positions(i,:) + delta * (u1 - u2) + alpha * (Mean_rule - Positions(i,:)) .* z;
            else
                NewPos1 = BestSol + delta * (u1 - u2) + alpha * (Mean_rule - Positions(i,:)) .* z;
            end

            % Boundary checking
            NewPos1 = min(max(NewPos1, lb), ub);
            NewFit1 = fobj(NewPos1);
            FEs = FEs + 1;

            if NewFit1 < Fitness(i)
                Positions(i,:) = NewPos1;
                Fitness(i) = NewFit1;
                if NewFit1 < BestScore
                    BestScore = NewFit1;
                    BestPos = NewPos1;
                end
            end
        end

        % Vector Combining Stage
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            if rand > 0.5
                Idx = randperm(nPop);
                Idx(Idx==i) = [];
                a = Idx(1); b = Idx(2);

                if Fitness(a) < Fitness(b)
                    u1 = Positions(a,:); u2 = Positions(b,:);
                else
                    u1 = Positions(b,:); u2 = Positions(a,:);
                end

                if rand > 0.5
                    NewPos2 = Positions(i,:) + rand * (BestSol - u1) + rand * (BetterSol - u2);
                else
                    NewPos2 = BestSol + rand * (BestSol - u1) + rand * (BetterSol - u2);
                end

                NewPos2 = min(max(NewPos2, lb), ub);
                NewFit2 = fobj(NewPos2);
                FEs = FEs + 1;

                if NewFit2 < Fitness(i)
                    Positions(i,:) = NewPos2;
                    Fitness(i) = NewFit2;
                    if NewFit2 < BestScore
                        BestScore = NewFit2;
                        BestPos = NewPos2;
                    end
                end
            end
        end

        % Local Search Stage (Enhanced for Electronics)
        pr = 0.1 + (0.3 * (1 - (t/MaxIter)^2));
        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            if rand < pr
                L = rand < 0.5;
                I1 = round(1 + rand);
                I2 = round(1 + rand);

                Idx = randperm(nPop);
                x_r1 = Positions(Idx(1),:);
                x_r2 = Positions(Idx(2),:);

                % Adaptive local search
                if L == 1
                    NewPos3 = BestSol + randn(1, dim) .* (Mean_rule - rand * BestSol);
                else
                    NewPos3 = Positions(i,:) + randn(1, dim) .* ((I1 * BestSol - I2 * Mean_rule) + rand * (x_r1 - x_r2));
                end

                NewPos3 = min(max(NewPos3, lb), ub);
                NewFit3 = fobj(NewPos3);
                FEs = FEs + 1;

                if NewFit3 < Fitness(i)
                    Positions(i,:) = NewPos3;
                    Fitness(i) = NewFit3;
                    if NewFit3 < BestScore
                        BestScore = NewFit3;
                        BestPos = NewPos3;
                    end
                end
            end
        end

        ConvergenceCurve(t) = BestScore;

        if mod(t, ceil(MaxIter/10)) == 0 || t == 1
            fprintf('INFO Iter %4d | Best = %.4e | FEs = %d\n', t, BestScore, FEs);
        end
    end
end