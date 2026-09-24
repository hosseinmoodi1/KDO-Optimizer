function [BestScore, BestPos, ConvergenceCurve] = SFO(nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)
% SFO - Swift Flight Optimizer
% -------------------------------------------------------------------------
% Reference:
%   Kareem, A. A., Abid, A. J., Hammood, D. A., Yaqoob, S. J., Husein, A.,
%   & Bereznychenko, V. (2025). Swift Flight Optimizer: a novel bio-inspired
%   optimization algorithm based on swift bird behavior. Scientific Reports,
%   15, 41597. https://doi.org/10.1038/s41598-025-25490-7
%
% Inputs / Outputs: standard SFO.m / INFO.m interface convention.
% -------------------------------------------------------------------------

    if numel(lb) == 1
        lb = lb * ones(1, dim);
        ub = ub * ones(1, dim);
    else
        lb = lb(:)';  ub = ub(:)';
    end

    FEs = 0;

    % ---- Algorithm parameters ----
    w         = 0.5;
    alpha_g   = 0.5;
    beta_l    = 0.57;
    psi1      = 1.45;
    psi2      = 1.45;
    StagLimit = 20;

    % ---- Initialization ----
    X   = rand(nPop, dim) .* (ub - lb) + lb;
    V   = zeros(nPop, dim);
    Fit = inf(nPop, 1);
    for i = 1:nPop
        if FEs >= MaxFEs, break; end
        Fit(i) = fobj(X(i,:));  FEs = FEs + 1;
    end

    if all(~isfinite(Fit))
        BestScore = inf;  BestPos = X(1, :);
        ConvergenceCurve = [];
        return;
    end

    Pbest    = X;
    PbestVal = Fit;
    [BestScore, idx] = min(Fit);
    BestPos  = X(idx, :);

    StallCounter     = zeros(nPop, 1);
    ConvergenceCurve = zeros(1, MaxIter);

    % ---- Main loop ----
    for t = 1:MaxIter
        if FEs >= MaxFEs
            ConvergenceCurve = ConvergenceCurve(1:t-1);
            break;
        end
        progress = t / MaxIter;

        for i = 1:nPop
            if FEs >= MaxFEs, break; end

            % ----- Stagnation-aware reinitialization -----
            if StallCounter(i) >= StagLimit
                X(i,:) = rand(1, dim) .* (ub - lb) + lb;
                V(i,:) = zeros(1, dim);
                StallCounter(i) = 0;
                if FEs < MaxFEs
                    Fit(i) = fobj(X(i,:));  FEs = FEs + 1;
                    if Fit(i) < PbestVal(i)
                        Pbest(i,:) = X(i,:);  PbestVal(i) = Fit(i);
                    end
                    if Fit(i) < BestScore
                        BestScore = Fit(i);  BestPos = X(i,:);
                    end
                end
                continue;
            end

            % ----- Mode selection -----
            r_mode = rand();
            if r_mode < (1 - progress)
                mode = 1;   % Glide
            elseif r_mode < (1 - progress) + 0.6 * progress
                mode = 2;   % Target
            else
                mode = 3;   % Micro
            end

            switch mode
                case 1  % Glide Mode
                    L = levy_step(dim, 1.5);
                    V(i,:) = w * V(i,:) + alpha_g * randn(1, dim) + beta_l * L;
                    Xnew   = X(i,:) + V(i,:);

                case 2  % Target Mode
                    r1 = rand(1, dim);  r2 = rand(1, dim);
                    V(i,:) = w * V(i,:) ...
                           + psi1 * r1 .* (Pbest(i,:) - X(i,:)) ...
                           + psi2 * r2 .* (BestPos   - X(i,:));
                    Xnew   = X(i,:) + V(i,:);

                case 3  % Micro Mode
                    r2 = rand(1, dim);
                    V(i,:) = w * V(i,:) ...
                           + psi2 * r2 .* (BestPos - X(i,:)) ...
                           + 0.01 * randn(1, dim);
                    Xnew   = X(i,:) + V(i,:);
            end

            Xnew   = min(max(Xnew, lb), ub);
            NewFit = fobj(Xnew);  FEs = FEs + 1;

            % ----- Greedy selection -----
            if NewFit < Fit(i)
                X(i,:) = Xnew;  Fit(i) = NewFit;
                StallCounter(i) = 0;
                if NewFit < PbestVal(i)
                    Pbest(i,:) = Xnew;  PbestVal(i) = NewFit;
                end
                if NewFit < BestScore
                    BestScore = NewFit;  BestPos = Xnew;
                end
            else
                StallCounter(i) = StallCounter(i) + 1;
            end
        end
        ConvergenceCurve(t) = BestScore;
    end
end

% -------------------------------------------------------------------------
function L = levy_step(d, beta)
    sigma = (gamma(1+beta) * sin(pi*beta/2) / ...
             (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    L = u ./ (abs(v).^(1/beta));
end