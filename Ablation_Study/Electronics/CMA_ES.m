function [BestScore, BestPos, ConvergenceCurve] = CMA_ES( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% CMA-ES: Covariance Matrix Adaptation Evolution Strategy
% Adapted for Electronics Benchmark (Hard Problems)
%
% Modified for MATLAB R2017b compatibility
% - Added MaxFEs support
% - Improved boundary handling
% - Adaptive sigma for electronics problems
% FIX: Corrected dimension handling for R2017b

    % Convert bounds to column vectors (FIX for R2017b)
    if numel(lb) == 1
        lb = lb * ones(dim, 1);
        ub = ub * ones(dim, 1);
    else
        lb = lb(:);
        ub = ub(:);
    end

    % Function evaluation counter
    FEs = 0;

    % --- CMA-ES Parameters ---
    lambda = nPop;
    mu = round(lambda / 2);
    weights = log(mu + 1/2) - log(1:mu)';
    weights = weights / sum(weights);
    mueff = sum(weights)^2 / sum(weights.^2);

    % Strategy parameters (Hansen's defaults)
    cc = (4 + mueff/dim) / (dim + 4 + 2*mueff/dim);
    cs = (mueff + 2) / (dim + mueff + 5);
    c1 = 2 / ((dim + 1.3)^2 + mueff);
    cmu = min(1 - c1, 2 * (mueff - 2 + 1/mueff) / ((dim + 2)^2 + mueff));
    damps = 1 + 2 * max(0, sqrt((mueff - 1)/(dim + 1)) - 1) + cs;

    % Initialization
    pc = zeros(dim, 1);
    ps = zeros(dim, 1);
    B = eye(dim);
    D = ones(dim, 1);
    C = B * diag(D.^2) * B';
    invsqrtC = B * diag(D.^-1) * B';
    eigeneval = 0;
    chiN = dim^0.5 * (1 - 1/(4*dim) + 1/(21*dim^2));

    % Initial solution (midpoint of bounds) - FIX: Use column vectors
    xmean = (lb + ub) / 2;
    
    % Adaptive initial sigma for electronics problems
    sigma = 0.3 * min(ub - lb);
    sigma = max(sigma, 0.001);

    BestScore = inf;
    BestPos = zeros(1, dim);
    ConvergenceCurve = zeros(1, MaxIter);

    % --- Main Loop ---
    for iter = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Generate and evaluate lambda offspring
        arx = zeros(dim, lambda);
        fitness = zeros(1, lambda);

        for k = 1:lambda
            arx(:, k) = xmean + sigma * B * (D .* randn(dim, 1));
            
            % Boundary checking
            arx(:, k) = min(max(arx(:, k), lb), ub);
            
            % Evaluation
            fitness(k) = fobj(arx(:, k)');
            FEs = FEs + 1;

            % Update Global Best
            if fitness(k) < BestScore
                BestScore = fitness(k);
                BestPos = arx(:, k)';
            end

            if FEs >= MaxFEs, break; end
        end

        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:iter-1);
            break;
        end

        % Sort by fitness and compute weighted mean
        [fitness, arindex] = sort(fitness);
        xold = xmean;
        arx_selected = arx(:, arindex(1:mu));
        xmean = arx_selected * weights;

        % Cumulation: Update evolution paths
        ps = (1-cs) * ps + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean - xold) / sigma;
        hsig = sum(ps.^2) / (1-(1-cs)^(2*iter/lambda)) / dim < 2 + 4/(dim+1);
        pc = (1-cc) * pc + hsig * sqrt(cc*(2-cc)*mueff) * (xmean - xold) / sigma;

        % Adapt covariance matrix C
        artmp = (1/sigma) * (arx_selected - repmat(xold, 1, mu));
        C = (1-c1-cmu) * C ...
            + c1 * (pc * pc' + (1-hsig) * cc*(2-cc) * C) ...
            + cmu * artmp * diag(weights) * artmp';

        % Adapt step size sigma
        sigma = sigma * exp((cs/damps) * (norm(ps)/chiN - 1));
        
        % Sigma bounds for stability
        sigma = max(sigma, 1e-6 * min(ub - lb));
        sigma = min(sigma, 0.5 * min(ub - lb));

        % Update B and D from C (every few iterations)
        if iter - eigeneval > lambda/(c1+cmu)/dim/10
            eigeneval = iter;
            C = triu(C) + triu(C,1)';
            [B, D] = eig(C);
            D = sqrt(max(diag(D), 1e-20));
            invsqrtC = B * diag(D.^-1) * B';
        end

        ConvergenceCurve(iter) = BestScore;

        if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
            fprintf('CMA_ES Iter %4d | Best = %.4e | FEs = %d\n', iter, BestScore, FEs);
        end
    end
end