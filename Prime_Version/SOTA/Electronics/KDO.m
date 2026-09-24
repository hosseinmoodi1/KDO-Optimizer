function [BestScore, BestPos, ConvergenceCurve] = KDO_Optimized( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)

% ==========================================================
% Karma-Dharma Optimization (KDO) - FULLY OPTIMIZED VERSION
% 
% Optimizations applied:
%   1. Vectorized Crossover (removed inner loop)
%   2. Pre-allocation for S_KR and S_TP arrays
%   3. Replaced normrnd with faster randn
%   4. Replaced min search with random replacement (O(1))
%   5. Single-line Boundary Control
%   6. Reduced redundant f_best/f_worst calculations
%   7. Optimized Nirvana Reset mechanism
%   8. Logical masking instead of array multiplication
%
% Input:
%   nPop      : Population size
%   MaxIter   : Maximum number of iterations
%   lb        : Lower bounds (scalar or 1 x dim vector)
%   ub        : Upper bounds (scalar or 1 x dim vector)
%   dim       : Dimension of the problem
%   fobj      : Objective function handle (minimization)
%   MaxFEs    : Maximum number of function evaluations
%
% Output:
%   BestScore         : Best fitness value found
%   BestPos           : Best solution found (1 x dim vector)
%   ConvergenceCurve  : Best fitness value at each iteration
% ==========================================================

%% Parameters
MemorySize = 25;      % Size of the Cosmic Memory

% --- Adaptive Parameters Setup ---
H = 5;                % Historical Memory Size
M_KR = 0.5 * ones(H, 1); % Memory of Karma Rate means
M_TP = 0.5 * ones(H, 1); % Memory of Transformation Probability means
k_mem = 1;            % Counter for memory update

NL_max = 100;         % Maximum Nirvana Limit
NL_min = 10;          % Minimum Nirvana Limit

%% Initialization
if numel(lb) == 1
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
else
    lb = lb(:)';
    ub = ub(:)';
end

FEs = 0;

% Vectorized population initialization
X = rand(nPop, dim) .* (ub - lb) + lb;
Fit = zeros(nPop, 1);
StallCounter = zeros(nPop, 1);

% Initial fitness evaluation
for i = 1:nPop
    Fit(i) = fobj(X(i,:)');
    FEs = FEs + 1;
    if FEs >= MaxFEs
        break;
    end
end

[BestScore, idx] = min(Fit);
BestPos = X(idx, :);

Pbest = X;
PbestVal = Fit;

% Global best tracking
GlobalBest = BestScore;
GlobalBestPos = BestPos;

%% Cosmic Memory Structure (Optimized with random replacement)
Memory.Delta = zeros(MemorySize, dim);
Memory.Gain = zeros(MemorySize, 1);
Memory.Count = 0;

ConvergenceCurve = zeros(1, MaxIter);

% Pre-allocated arrays for successful parameters (avoid dynamic growth)
S_KR_buffer = zeros(nPop, 1);
S_TP_buffer = zeros(nPop, 1);

%% Main Loop
for iter = 1:MaxIter

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
    
    % Get global best and worst (once per iteration)
    f_best = min(Fit);
    f_worst = max(Fit);
    
    for i = 1:nPop
        
        if FEs >= MaxFEs
            break;
        end
        
        r = rand;
        
        % ================= KARMA PHASE =================
        if r < KR_pop(i)
            
            if f_worst ~= f_best
                alpha = (Fit(i) - f_best) / (f_worst - f_best);
            else
                alpha = 0.5;
            end
            beta = 1 - alpha;
            
            Step = alpha .* (Pbest(i,:) - X(i,:)) + ...
                   beta  .* (GlobalBestPos - X(i,:));
               
            Step = Step .* randn(1, dim);
            proposed_pos = X(i,:) + Step;
            
        % ================= DHARMA PHASE =================
        else
            
            if Memory.Count >= 2
                % Random selection from memory (O(1) - no min search)
                idx1 = randi(Memory.Count);
                idx2 = randi(Memory.Count);
                
                w1 = Memory.Gain(idx1);
                w2 = Memory.Gain(idx2);
                
                Transform = (w1 * Memory.Delta(idx1,:) + ...
                             w2 * Memory.Delta(idx2,:)) ...
                             / (w1 + w2 + eps);
            else
                Transform = randn(1, dim);
            end
            
            gamma = 0.5 + rand/2;
            proposed_pos = GlobalBestPos + gamma * Transform;
        end
        
        % --- VECTORIZED KARMIC INHERITANCE (Crossover) ---
        % No for-loop! Uses logical masking for O(1) complexity
        j_rand = randi(dim);
        crossover_mask = (rand(1, dim) < TP_pop(i));
        crossover_mask(j_rand) = 1;  % Ensure at least one dimension changes
        
        % Apply mask efficiently
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
            
            % Store successful parameters (using pre-allocated buffer)
            success_count = success_count + 1;
            S_KR_buffer(success_count) = KR_pop(i);
            S_TP_buffer(success_count) = TP_pop(i);
            
            % Update Cosmic Memory (O(1) random replacement)
            if Memory.Count < MemorySize
                Memory.Count = Memory.Count + 1;
                idx_mem = Memory.Count;
            else
                % Random replacement instead of min search (much faster)
                idx_mem = randi(MemorySize);
            end
            Memory.Delta(idx_mem, :) = Delta;
            Memory.Gain(idx_mem) = improvement;
            
            % Update personal best
            if NewFit < PbestVal(i)
                Pbest(i,:) = NewX;
                PbestVal(i) = NewFit;
            end
            
            % Update global best
            if NewFit < GlobalBest
                GlobalBest = NewFit;
                GlobalBestPos = NewX;
            end
            
        else
            StallCounter(i) = StallCounter(i) + 1;
        end
        
        % ================= OPTIMIZED NIRVANA RESET =================
        if StallCounter(i) >= NL
            
            % Random reset
            X(i,:) = rand(1, dim) .* (ub - lb) + lb;
            
            if FEs < MaxFEs
                Fit(i) = fobj(X(i,:)');
                FEs = FEs + 1;
                
                Pbest(i,:) = X(i,:);
                PbestVal(i) = Fit(i);
                StallCounter(i) = 0;
                
                if Fit(i) < GlobalBest
                    GlobalBest = Fit(i);
                    GlobalBestPos = X(i,:);
                end
            end
        end
        
        % Update current global best for next iteration
        BestScore = GlobalBest;
        BestPos = GlobalBestPos;
        
    end % End of population loop
    
    % --- Update historical memory using pre-allocated buffers ---
    if success_count > 0
        M_KR(k_mem) = mean(S_KR_buffer(1:success_count));
        M_TP(k_mem) = mean(S_TP_buffer(1:success_count));
        
        k_mem = k_mem + 1;
        if k_mem > H
            k_mem = 1;
        end
    end
    
    ConvergenceCurve(iter) = BestScore;
    
    % Progress reporting (only at 10% intervals)
    if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
        fprintf('Iter %4d | Best = %.4e | FEs = %d\n', iter, BestScore, FEs);
    end
end

end