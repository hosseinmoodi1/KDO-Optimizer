function [Best_score, Best_pos, cg_curve] = AVOA(PopSize, MaxIter, lb, ub, dim, fobj, MaxFEs)
    % AVOA: African Vultures Optimization Algorithm (2021) - Pure Version
    % Reference: Abdollahzadeh, B., Gharehchopogh, F. S., & Mirjalili, S. (2021).
    % African vultures optimization algorithm: A new nature-inspired
    % metaheuristic algorithm for global optimization problems.
    % Computers & Industrial Engineering, 158, 107408.
    
    % Initialization
    Positions = initialization(PopSize, dim, ub, lb);
    Fitness = zeros(1, PopSize);
    FEs = 0;
    
    for i = 1:PopSize
        Fitness(i) = fobj(Positions(i, :));
        FEs = FEs + 1;
        if FEs >= MaxFEs, break; end
    end
    
    % Sort and identify best two vultures
    [~, SortOrder] = sort(Fitness);
    Best_pos = Positions(SortOrder(1), :);
    Best_pos2 = Positions(SortOrder(2), :);
    Best_score = Fitness(SortOrder(1));
    
    cg_curve = zeros(1, MaxIter);
    
    % Main loop
    for t = 1:MaxIter
        if FEs >= MaxFEs
            cg_curve = cg_curve(1:t-1);
            break;
        end
        
        % Eq. (2): Parameter a (decreases from 2 to 0)
        a = 2 * (1 - t / MaxIter);
        
        % Pre-calculate probabilities (Eq. 5-6)
        p1 = 0.6 * (1 - t / MaxIter);
        p2 = 0.4 * (1 - t / MaxIter);
        p3 = 0.1 * (1 - t / MaxIter);
        
        for i = 1:PopSize
            if FEs >= MaxFEs, break; end
            
            % Eq. (3): Flight speed F
            F = (2 * rand + 1) * (1 - t / MaxIter)^(2*rand) + a * (rand - 1);
            
            % Eq. (7): Select target (R) between Best and Second Best
            L1 = 0.8; L2 = 0.2;
            if rand <= L1
                R = Best_pos;
            else
                R = Best_pos2;
            end
            
            if abs(F) >= 1
                % ===== EXPLORATION PHASE =====
                if rand < p1
                    % Eq. (8): Random selection
                    X = 2 * rand;
                    Di = abs(X * R - Positions(i, :));
                    NewPos = R - Di .* F;
                else
                    % Eq. (9): Random long-distance movement
                    NewPos = R - F + rand(1, dim) .* ((ub - lb) .* rand(1, dim) + lb);
                end
                
            elseif abs(F) >= 0.5
                % ===== EXPLOITATION PHASE 1 =====
                if rand < p2
                    % Eq. (10): Rotating flight (sine-cosine)
                    s1 = R .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* cos(Positions(i, :));
                    s2 = R .* (rand(1, dim) .* Positions(i, :) / (2*pi)) .* sin(Positions(i, :));
                    NewPos = R - (s1 + s2);
                else
                    % Eq. (11): Aggressive competition with Levy flight
                    Di = abs(R - Positions(i, :));
                    L = Levy(dim);
                    NewPos = Di - F .* L;
                end
                
            else
                % ===== EXPLOITATION PHASE 2 =====
                if rand < p3
                    % Eq. (12): Accumulation around food (mean of 3 best)
                    NewPos = (Best_pos + Best_pos2 + Positions(i, :)) / 3;
                else
                    % Eq. (13): Siege fight with Levy and displacement
                    d_t = Positions(i, :) - Best_pos;
                    Di = abs(R - Positions(i, :));
                    L = Levy(dim);
                    NewPos = Di .* F .* L + d_t;
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
        
        % Re-sort population and update best two vultures
        [~, SortOrder] = sort(Fitness);
        Best_pos = Positions(SortOrder(1), :);
        Best_pos2 = Positions(SortOrder(2), :);
        Best_score = Fitness(SortOrder(1));
        
        cg_curve(t) = Best_score;
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