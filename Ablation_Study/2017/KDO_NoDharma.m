function [BestScore, BestPos, ConvergenceCurve] = KDO_NoDharma( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)

%% Parameters
MemorySize = 25;      % Size of the Cosmic Memory (Updates but never used for Dharma)

% --- Adaptive Parameters Setup ---
H = 5;                % Historical Memory Size
M_KR = 0.7 * ones(H, 1); % Memory of Karma Rate means
M_TP = 0.5 * ones(H, 1); % Memory of Transformation Probability means
k_mem = 1;            % Counter for memory update

NL_max = 100;         % Maximum Nirvana Limit
NL_min = 10;          % Minimum Nirvana Limit

% Ensure boundaries are row vectors
lb = lb(:)'; 
ub = ub(:)';

FEs = 0; % Function Evaluations Counter

%% Initialization
% Initialize population positions randomly within boundaries
X = rand(nPop, dim) .* (ub - lb) + lb;
Fit = zeros(nPop, 1);
StallCounter = zeros(nPop, 1); % Tracks consecutive non-improvements for each agent

% Evaluate initial population
for i = 1:nPop
    Fit(i) = fobj(X(i,:)');
    FEs = FEs + 1;
    if FEs >= MaxFEs
        break;
    end
end

% Find initial global best
[BestScore, idx] = min(Fit);
BestPos = X(idx, :);

% Initialize personal bests
Pbest    = X;
PbestVal = Fit;

%% Cosmic Memory Structure
Memory.Delta = zeros(MemorySize, dim);
Memory.Gain  = zeros(MemorySize, 1);
Memory.Count = 0;

% Array to store the convergence curve
ConvergenceCurve = zeros(1, MaxIter);

% Pre-allocated buffer for successful parameters (avoid dynamic growth)
S_KR_buffer = zeros(nPop, 1);
S_TP_buffer = zeros(nPop, 1);

%% Main Loop
for iter = 1:MaxIter

    % Stop if maximum function evaluations are reached
    if FEs >= MaxFEs
        ConvergenceCurve = ConvergenceCurve(1:iter-1);
        break;
    end
    
    % --- Dynamic Nirvana Limit (decreasing for exploitation focus) ---
    NL = round(NL_max - (NL_max - NL_min) * (iter / MaxIter));
    
    % --- Adaptive Parameters for all agents ---
    mem_rand_idx = randi(H, nPop, 1);
    
    % Generate Karma Rate using randn (faster than normrnd)
    KR_pop = M_KR(mem_rand_idx) + 0.1 * randn(nPop, 1);
    KR_pop = min(max(KR_pop, 0), 1);
    
    % Generate Transformation Probability using randn
    TP_pop = M_TP(mem_rand_idx) + 0.1 * randn(nPop, 1);
    TP_pop = min(max(TP_pop, 0), 1);
    
    % Reset success counter for this iteration
    success_count = 0;
    
    % Find best and worst fitness of the current population (once per iteration)
    f_best  = min(Fit);
    f_worst = max(Fit);
    
    for i = 1:nPop
        
        % Stop evaluating if FEs limit is reached
        if FEs >= MaxFEs
            break;
        end
        
        % ================= ONLY KARMA PHASE (NO DHARMA) =================
        % Adaptive coefficients based on current fitness ranking
        if f_worst ~= f_best
            alpha = (Fit(i) - f_best) / (f_worst - f_best);
        else
            alpha = 0.5;
        end
        beta = 1 - alpha;
        
        % Generate Karma Step using personal and global bests
        Step = alpha .* (Pbest(i,:) - X(i,:)) + ...
               beta  .* (BestPos - X(i,:));
           
        % Apply random perturbation
        Step = Step .* randn(1, dim);
        proposed_pos = X(i,:) + Step;
        
        % --- VECTORIZED KARMIC INHERITANCE (Crossover) ---
        j_rand = randi(dim);
        crossover_mask = (rand(1, dim) < TP_pop(i));
        crossover_mask(j_rand) = 1;
        
        NewX = X(i,:);
        NewX(crossover_mask) = proposed_pos(crossover_mask);
        
        % ================= SINGLE-LINE BOUNDARY CONTROL =================
        NewX = min(max(NewX, lb), ub);
        
        % ================= FITNESS EVALUATION =================
        NewFit = fobj(NewX');
        FEs = FEs + 1;
        
        % ================= SELECTION & UPDATE =================
        if NewFit < Fit(i)
            
            improvement = abs(Fit(i) - NewFit);
            Delta = NewX - X(i,:);
            
            X(i,:) = NewX;
            Fit(i) = NewFit;
            
            StallCounter(i) = 0;
            
            success_count = success_count + 1;
            S_KR_buffer(success_count) = KR_pop(i);
            S_TP_buffer(success_count) = TP_pop(i);
            
            if Memory.Count < MemorySize
                Memory.Count = Memory.Count + 1;
                idx_mem = Memory.Count;
            else
                idx_mem = randi(MemorySize);
            end
            Memory.Delta(idx_mem, :) = Delta;
            Memory.Gain(idx_mem) = improvement;
            
            if NewFit < PbestVal(i)
                Pbest(i,:)  = NewX;
                PbestVal(i) = NewFit;
            end
            
            if NewFit < BestScore
                BestScore = NewFit;
                BestPos   = NewX;
            end
            
        else
            StallCounter(i) = StallCounter(i) + 1;
        end
        
        % ================= OPTIMIZED NIRVANA RESET MECHANISM =================
        if StallCounter(i) >= NL
            X(i,:) = rand(1, dim) .* (ub - lb) + lb;
            if FEs < MaxFEs
                Fit(i) = fobj(X(i,:)');
                FEs = FEs + 1;
                
                Pbest(i,:)  = X(i,:);
                PbestVal(i) = Fit(i);
                StallCounter(i) = 0;
                
                if Fit(i) < BestScore
                    BestScore = Fit(i);
                    BestPos   = X(i,:);
                end
            end
        end
        % ===========================================================
    end
    
    % --- Update historical memory ---
    if success_count > 0
        M_KR(k_mem) = mean(S_KR_buffer(1:success_count));
        M_TP(k_mem) = mean(S_TP_buffer(1:success_count));
        
        k_mem = k_mem + 1;
        if k_mem > H
            k_mem = 1;
        end
    end
    
    ConvergenceCurve(iter) = BestScore;
    
    if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
        fprintf('Iter %4d | Best = %.4e | FEs = %d\n', iter, BestScore, FEs);
    end
end
end
