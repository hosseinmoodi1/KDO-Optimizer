function [Best_score, Best_pos, cg_curve] = COA(PopSize, MaxIter, lb, ub, dim, fobj, MaxFEs)
    % COA: Crayfish Optimization Algorithm (2023) - Pure Version
    % Reference: Jia, H., Rao, H., Wen, C., & Mirjalili, S. (2023).
    % Crayfish optimization algorithm. Artificial Intelligence Review.
    
    % Initialization
    Positions = initialization(PopSize, dim, ub, lb);
    Fitness = zeros(1, PopSize);
    FEs = 0;
    
    for i = 1:PopSize
        Fitness(i) = fobj(Positions(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end
    
    [Best_score, best_idx] = min(Fitness);
    Best_pos = Positions(best_idx, :);
    
    cg_curve = zeros(1, MaxIter);
    
    % Main loop
    for t = 1:MaxIter
        if FEs >= MaxFEs
            cg_curve = cg_curve(1:t-1);
            break;
        end
        
        % Eq. (3): Temperature calculation
        Temp = rand * 15 + 20;
        
        for i = 1:PopSize
            if FEs >= MaxFEs, break; end
            
            % Eq. (4): Adaptive step parameter k
            k = 2 * rand - 1;
            
            if Temp > 30
                % ===== Summer Resort (Exploration) =====
                % Eq. (6): Cave location
                X_cave = (Best_pos + Positions(i, :)) / 2;
                
                if rand < 0.5
                    % Eq. (7): Craw into cave
                    C = 2 - t / MaxIter;
                    NewPos = Positions(i, :) + C * rand(1, dim) .* (X_cave - Positions(i, :));
                else
                    % Eq. (8): Competition for cave
                    Idx = randperm(PopSize, 2);
                    NewPos = Positions(i, :) - Positions(Idx(1), :) + X_cave;
                end
                
            else
                % ===== Competition (Exploitation) =====
                % Eq. (9): Foraging amount p
                p = 0.2 + 0.3 * (1 - t / MaxIter);
                
                % Eq. (10): Food attraction Q
                Q = 2 * rand(1, dim) - 1;
                
                if rand < p
                    % Eq. (11): Foraging behavior
                    NewPos = Positions(i, :) + k * (Best_pos - Positions(i, :)) .* Q;
                else
                    % Eq. (12): Fighting behavior
                    z = round(rand(1, dim));
                    NewPos = Positions(i, :) + k * (Best_pos - z .* Positions(i, :));
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
                
                if NewFit < Best_score
                    Best_score = NewFit;
                    Best_pos = NewPos;
                end
            end
        end
        
        cg_curve(t) = Best_score;
    end
end

function Positions = initialization(SearchAgents_no, dim, ub, lb)
    Boundary_no = size(ub, 2);
    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
    else
        Positions = zeros(SearchAgents_no, dim);
        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            Positions(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
        end
    end
end