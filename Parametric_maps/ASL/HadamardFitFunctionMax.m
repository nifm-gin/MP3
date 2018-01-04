function maxX=HadamardFitFunctionMax(xfix,xvar)

f=xvar(2); % mL/g/s
deltaA=xvar(1); %s

if numel(xvar)>2    
    %T1t = xvar(3); %s
    T1app = xvar(3); %s
else
    T1app = xfix.T1app; %s
end
% if strcmp(xfix.model,'Buxton') == 1
%     deltaFromAtoTissue = xvar(3);
% end

X = [];
delta = deltaA; %s
M0 = xfix.M0;
SliceDecM0 = xfix.SliceDecM0;
Effective_PLD = xfix.step/10:xfix.step/10:max(xfix.Effective_PLD); %s Effective_PLD corrected for the interslice time;
PostLabelTime = xfix.PostLabelTime; %s
alpha=xfix.alpha; % Inversion efficiency
Labeling_Duration1SB = xfix.SubboliDuration'; %s Labeling duration of one subbolus
tL = Labeling_Duration1SB(1);  % (ms) labeling time of the subbolus (each SB has the same duration)
T1blood = xfix.T1sang; %s
lambda = 0.9; % blood:brain partition coefficient for water


% if numel(xvar)>2    
%     T1app = lambda * T1t/(lambda + T1t * f);
% end

for num = 1:numel(Effective_PLD)
    current_pld = Effective_PLD(num);
    
    
    %% Thomas 2006
    
    if strcmp(xfix.model,'Thomas') == 1
        if (current_pld < deltaA)
            X = [X 0];
        elseif current_pld > deltaA     
           X = [X   2 * alpha * f * M0 * SliceDecM0 /lambda *  T1app *exp(-delta/T1blood) * (exp(min(delta-current_pld,0)/T1app)-...
               exp(-current_pld/T1app)) + T1blood * (exp((min(deltaA-current_pld,0)-deltaA)/T1blood) - exp((min(delta-current_pld,0)-delta)/T1blood))];
        end
    elseif strcmp(xfix.model,'Wells') == 1
        %% Wells 2010
            if tL + current_pld > delta
                X = [X    2 / lambda * f * M0 * SliceDecM0 * alpha * (exp(-delta/T1blood)* T1app * ...
                    (exp(min(delta-current_pld,0)/T1app) - exp((delta - tL - current_pld)/T1app)) ...
                    + T1blood * (exp((min(deltaA-current_pld,0)-deltaA)/T1blood) - exp((min(delta-current_pld,0)-delta)/T1blood)))];
            elseif tL + current_pld <= delta
                X  = [X    2 / lambda * f* M0* SliceDecM0 * alpha *T1blood* (exp((min(deltaA-current_pld,0)-deltaA)/T1blood) ...
                    - exp(-(tL+current_pld)/T1blood)) * min(deltaA-current_pld-tL,0)/(deltaA-current_pld-tL)  ];
            end
        
    elseif strcmp(xfix.model,'Buxton') == 1
        %% Buxton
        if (current_pld < deltaA)
            X = [X 0];
        elseif current_pld >= deltaA && current_pld < deltaA + tL
            X = [X   2 * M0 * SliceDecM0 /lambda * f * T1app * alpha * exp(-deltaA/T1blood) * (1- exp(-(current_pld-deltaA)/T1app))];
        elseif tL + deltaA  <= current_pld %&& current_pld <= deltaA + deltaFromAtoTissue
            X = [X   2 * M0 * SliceDecM0 /lambda * f * T1app * alpha * exp(-deltaA/T1blood) * exp(-(current_pld-tL-deltaA)/T1app) *...
                (1- exp(-tL/T1app))];
%         elseif current_pld> deltaA+deltaFromAtoTissue
%             X = [X   2 * M0 * SliceDecM0 /lambda * f * T1app * alpha * exp(-deltaA/T1blood) * exp(-(current_pld-tL-deltaA)/T1app) *...
%                 (1- exp(-tL/T1app))];
        end
        
    end
end

maxX = max(X);
%figure(100), hold on , plot((Effective_PLD)*1000,X,'g');
