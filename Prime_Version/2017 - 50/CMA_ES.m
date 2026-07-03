function [Best_score, Best_pos, cg_curve] = CMA_ES(PopSize, MaxIter, lb, ub, dim, fobj, MaxFEs)
    % CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
    % Adapted for standard metaheuristic function signatures
    
    % --- Parameters ---
    lambda = PopSize; % Population size
    mu = round(lambda / 2); % Number of parents/points for recombination
    weights = log(mu + 1/2) - log(1:mu)'; % Mu recombination weights
    weights = weights / sum(weights); % Normalize recombination weights
    mueff = sum(weights)^2 / sum(weights.^2); % Variance-effectiveness
    
    % Strategy parameter setting
    cc = (4 + mueff/dim) / (dim + 4 + 2*mueff/dim);
    cs = (mueff + 2) / (dim + mueff + 5);
    c1 = 2 / ((dim + 1.3)^2 + mueff);
    cmu = min(1 - c1, 2 * (mueff - 2 + 1/mueff) / ((dim + 2)^2 + mueff));
    damps = 1 + 2 * max(0, sqrt((mueff - 1)/(dim + 1)) - 1) + cs;
    
    % Initialization
    if length(lb) == 1
        lb = repmat(lb, 1, dim);
        ub = repmat(ub, 1, dim);
    end
    
    pc = zeros(dim, 1);
    ps = zeros(dim, 1);
    B = eye(dim, dim);
    D = ones(dim, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^-1) * B';
    eigeneval = 0; % Track update of B and D
    chiN = dim^0.5 * (1 - 1/(4*dim) + 1/(21*dim^2));
    
    % Initial solution
    xmean = (lb + ub)' / 2; % Start in the middle of bounds
    sigma = 0.3 * (ub(1) - lb(1)); % Initial step size (approx 1/3 of range)
    
    Best_score = inf;
    Best_pos = zeros(1, dim);
    cg_curve = zeros(1, MaxIter);
    
    % --- Main Loop ---
    for iter = 1:MaxIter
        % Generate and evaluate lambda offspring
        arx = zeros(dim, lambda);
        fitness = zeros(1, lambda);
        
        for k = 1:lambda
            arx(:, k) = xmean + sigma * B * (D .* randn(dim, 1));
            
            % Boundary checking
            arx(:, k) = min(max(arx(:, k), lb'), ub');
            
            % Evaluation
            fitness(k) = fobj(arx(:, k)');
            
            % Update Global Best
            if fitness(k) < Best_score
                Best_score = fitness(k);
                Best_pos = arx(:, k)';
            end
        end
        
        % Sort by fitness and compute weighted mean into xmean
        [fitness, arindex] = sort(fitness);
        xold = xmean;
        arx = arx(:, arindex(1:mu));
        xmean = arx * weights;
        
        % Cumulation: Update evolution paths
        ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean - xold) / sigma;
        hsig = sum(ps.^2) / (1-(1-cs)^(2*iter/lambda)) / dim < 2 + 4/(dim+1);
        pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean - xold) / sigma;
        
        % Adapt covariance matrix C
        artmp = (1/sigma) * (arx - repmat(xold, 1, mu));
        C = (1-c1-cmu) * C ...
            + c1 * (pc * pc' ...
            + (1-hsig) * cc*(2-cc) * C) ...
            + cmu * artmp * diag(weights) * artmp';
        
        % Adapt step size sigma
        sigma = sigma * exp((cs/damps) * (norm(ps)/chiN - 1));
        
        % Update B and D from C
        if iter - eigeneval > lambda/(c1+cmu)/dim/10
            eigeneval = iter;
            C = triu(C) + triu(C,1)'; % Enforce symmetry
            [B, D] = eig(C);
            D = sqrt(diag(D));
            invsqrtC = B * diag(D.^-1) * B';
        end
        
        cg_curve(iter) = Best_score;
    end
end
