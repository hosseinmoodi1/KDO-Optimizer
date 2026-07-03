function [Best_score, Best_pos, cg_curve] = RUN(PopSize, MaxIter, lb, ub, dim, fobj, MaxFEs)
    % RUN: Runge Kutta Optimizer (Corrected Version)
    % Reference: Ahmadianfar, I., Heidari, A. A., Gandomi, A. H., Chu, X., & Chen, H. (2021).
    % RUN beyond the metaphor: An efficient optimization algorithm 
    % based on Runge Kutta method. Expert Systems with Applications, 181, 115079.
    %
   
    
    % Initialize boundaries
    if length(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    end

    % Initialization
    X = zeros(PopSize, dim);
    for i = 1:dim
        X(:, i) = rand(PopSize, 1) .* (ub(i) - lb(i)) + lb(i);
    end

    Fitness = zeros(PopSize, 1);
    for i = 1:PopSize
        Fitness(i) = fobj(X(i, :));
    end
    FEs = PopSize;

    % Find initial best
    [Best_score, ind] = min(Fitness);
    Best_pos = X(ind, :);
    
    cg_curve = zeros(1, MaxIter);
    cg_curve(1) = Best_score;
    it = 1;

    % Algorithmic parameters (as per original paper)
    a = 20;  
    b = 12;  
    c = 2;

    % Main optimization loop
    while FEs < MaxFEs
        % Update adaptive parameters
        f = 0.2 + 0.3 * rand(); % Adaptive scaling factor (Eq. 14)
        
        % Sort population to find best and worst
        [~, SortIdx] = sort(Fitness);
        Xb = X(SortIdx(1), :);   % Best solution
        
        % Calculate population average
        Xavg = mean(X, 1);
        
        for i = 1:PopSize
            if FEs >= MaxFEs
                break;
            end
            
            % Select three random distinct individuals
            Id = randperm(PopSize);
            Id(Id == i) = [];
            r1 = Id(1);
            r2 = Id(2);
            r3 = Id(3);
            
            % Calculate XC (center of three random points) - Eq. (15)
            XC = (X(i, :) + X(r1, :) + X(r2, :)) / 3;
            
            % Determine Delta X - Eq. (16)
            r = rand();
            if r < 0.5
                DX = XC - Xavg;       % Exploration: move toward average
            else
                DX = XC - Xb;         % Exploitation: move toward best
            end
            
            % Generate random coefficients for RK mechanism
            A = randn(1, dim);           % Normal distribution
            B = exp(-rand(1, dim));      % Exponential decay
            
            % Search mechanism (SM) based on difference of random solutions
            SM = f * (X(r1, :) - X(r2, :));  % Eq. (17)
            
            % Runge-Kutta coefficients (Eq. 18-21)
            k1 = (1 / (2 * a)) * (A .* DX + B .* SM);
            k2 = (1 / (2 * a)) * (A .* (DX + c * k1) + B .* SM);
            k3 = (1 / (2 * a)) * (A .* (DX + c * k2) + B .* SM);
            k4 = (1 / (2 * a)) * (A .* (DX + c * k3) + B .* SM);
            
            % Calculate search factor SF (Eq. 22)
            SF = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            
            % Update position - Eq. (23)
            X_new = X(i, :) + SF .* X(i, :);
            
            % Boundary control
            X_new = min(max(X_new, lb), ub);
            
            % Evaluate new solution
            Fit_new = fobj(X_new);
            FEs = FEs + 1;
            
            % Greedy selection
            if Fit_new < Fitness(i)
                X(i, :) = X_new;
                Fitness(i) = Fit_new;
            end
            
            % Enhanced Solution Quality (ESQ) Phase - Eq. (24-26)
            if rand() < 0.5 && FEs < MaxFEs
                W = rand();
                if W < 0.5
                    % Local search around best solution
                    X_new2 = Xb + randn(1, dim) .* (X(r1, :) - X(r2, :));
                else
                    % Search around current solution
                    X_new2 = X(i, :) + randn(1, dim) .* (Xb - X(i, :));
                end
                
                % Boundary control for ESQ
                X_new2 = min(max(X_new2, lb), ub);
                
                % Evaluate ESQ solution
                Fit_new2 = fobj(X_new2);
                FEs = FEs + 1;
                
                % Greedy selection for ESQ
                if Fit_new2 < Fitness(i)
                    X(i, :) = X_new2;
                    Fitness(i) = Fit_new2;
                end
            end
        end
        
        % Update global best
        [current_best, current_ind] = min(Fitness);
        if current_best < Best_score
            Best_score = current_best;
            Best_pos = X(current_ind, :);
        end
        
        % Store convergence
        it = it + 1;
        if it <= MaxIter
            cg_curve(it) = Best_score;
        else
            break;
        end
    end
    
    % Trim convergence curve if necessary
    cg_curve = cg_curve(1:min(it, MaxIter));
end