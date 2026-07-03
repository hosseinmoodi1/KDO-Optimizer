function [BestScore, BestPos, ConvergenceCurve] = WOA( ...
    nPop, MaxIter, lb, ub, dim, fobj, MaxFEs)

% ==========================================================
% Whale Optimization Algorithm (WOA)
% Standardized Format (Compatible with KDO Framework)
% Based on: Mirjalili & Lewis (2016)
% ==========================================================

%% Bounds
lb = lb(:)';
ub = ub(:)';

%% Function evaluations
FEs = 0;

%% Initialization
X = rand(nPop,dim).*(ub-lb) + lb;
Fit = zeros(nPop,1);

for i = 1:nPop
    Fit(i) = fobj(X(i,:)');
    FEs = FEs + 1;
    
    if FEs >= MaxFEs
        break;
    end
end

[BestScore, idx] = min(Fit);
BestPos = X(idx,:);

ConvergenceCurve = zeros(1,MaxIter);

%% Main Loop
for iter = 1:MaxIter
    
    if FEs >= MaxFEs
        ConvergenceCurve = ConvergenceCurve(1:iter-1);
        break;
    end
    
    % a decreases linearly from 2 to 0
    a = 2 - iter*(2/MaxIter);
    
    for i = 1:nPop
        
        if FEs >= MaxFEs
            break;
        end
        
        r1 = rand;
        r2 = rand;
        
        A = 2*a*r1 - a;
        C = 2*r2;
        
        b = 1;
        l = -1 + 2*rand;
        p = rand;
        
        if p < 0.5
            
            if abs(A) >= 1
                
                % Exploration
                rand_idx = randi(nPop);
                X_rand = X(rand_idx,:);
                
                D = abs(C*X_rand - X(i,:));
                NewX = X_rand - A*D;
                
            else
                
                % Exploitation
                D = abs(C*BestPos - X(i,:));
                NewX = BestPos - A*D;
                
            end
            
        else
            
            % Spiral update
            D = abs(BestPos - X(i,:));
            NewX = D .* exp(b*l) .* cos(2*pi*l) + BestPos;
            
        end
        
        % Boundary control
        NewX = max(NewX, lb);
        NewX = min(NewX, ub);
        
        % Fitness evaluation
        NewFit = fobj(NewX');
        FEs = FEs + 1;
        
        % Selection
        if NewFit < Fit(i)
            X(i,:) = NewX;
            Fit(i) = NewFit;
        end
        
        % Update global best
        if NewFit < BestScore
            BestScore = NewFit;
            BestPos   = NewX;
        end
    end
    
    ConvergenceCurve(iter) = BestScore;
    
    if mod(iter, ceil(MaxIter/10)) == 0 || iter == 1
        fprintf('Iter %4d | Best = %.4e | FEs = %d\n', ...
                iter, BestScore, FEs);
    end
end

end
