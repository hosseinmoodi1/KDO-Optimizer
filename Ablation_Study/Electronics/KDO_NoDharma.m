function [BestScore, BestPos, ConvergenceCurve] = KDO_NoDharma( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)

%% Parameters
H = 5;                
M_KR = 0.5 * ones(H, 1); 
M_TP = 0.5 * ones(H, 1); 
k_mem = 1;            

NL_max = 100;         
NL_min = 10;          

%% Initialization
if numel(lb) == 1
    lb = lb * ones(1, dim);
    ub = ub * ones(1, dim);
else
    lb = lb(:)';
    ub = ub(:)';
end

FEs = 0;
X = rand(nPop, dim) .* (ub - lb) + lb;
Fit = zeros(nPop, 1);
StallCounter = zeros(nPop, 1);

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

GlobalBest = BestScore;
GlobalBestPos = BestPos;

ConvergenceCurve = zeros(1, MaxIter);

S_KR_buffer = zeros(nPop, 1);
S_TP_buffer = zeros(nPop, 1);

%% Main Loop
for iter = 1:MaxIter

    if FEs >= MaxFEs
        ConvergenceCurve = ConvergenceCurve(1:iter-1);
        break;
    end
    
    NL = round(NL_max - (NL_max - NL_min) * (iter / MaxIter));
    
    mem_rand_idx = randi(H, nPop, 1);
    
    KR_pop = M_KR(mem_rand_idx) + 0.1 * randn(nPop, 1);
    KR_pop = min(max(KR_pop, 0), 1);
    
    TP_pop = M_TP(mem_rand_idx) + 0.1 * randn(nPop, 1);
    TP_pop = min(max(TP_pop, 0), 1);
    
    success_count = 0;
    
    f_best = min(Fit);
    f_worst = max(Fit);
    
    for i = 1:nPop
        
        if FEs >= MaxFEs
            break;
        end
        
        % ================= KARMA PHASE ONLY =================
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
        
        % --- Crossover ---
        j_rand = randi(dim);
        crossover_mask = (rand(1, dim) < TP_pop(i));
        crossover_mask(j_rand) = 1;  
        
        NewX = X(i,:);
        NewX(crossover_mask) = proposed_pos(crossover_mask);
        
        NewX = min(max(NewX, lb), ub);
        
        NewFit = fobj(NewX');
        FEs = FEs + 1;
        
        % ================= SELECTION & UPDATE =================
        if NewFit < Fit(i)
            X(i,:) = NewX;
            Fit(i) = NewFit;
            StallCounter(i) = 0;
            
            success_count = success_count + 1;
            S_KR_buffer(success_count) = KR_pop(i);
            S_TP_buffer(success_count) = TP_pop(i);
            
            if NewFit < PbestVal(i)
                Pbest(i,:) = NewX;
                PbestVal(i) = NewFit;
            end
            
            if NewFit < GlobalBest
                GlobalBest = NewFit;
                GlobalBestPos = NewX;
            end
            
        else
            StallCounter(i) = StallCounter(i) + 1;
        end
        
        % ================= NIRVANA RESET =================
        if StallCounter(i) >= NL
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
        
        BestScore = GlobalBest;
        BestPos = GlobalBestPos;
        
    end 
    
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
