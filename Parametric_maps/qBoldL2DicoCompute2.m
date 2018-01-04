function DataOut = qBoldL2DicoCompute2(S,ADC,Dico)

% Process the signal by comparing the input (S,ADC) to the dictionary of
% curves (Dico). Raw signal S is compared in term of normalized root mean
% squared deviation (nRMSD) to the dictionary.
% Input: - S: Array of signals. N x T where length of T must the same than
% size(Dico.dat,2)
%        - ADC: Array of corresponding ADC values. N x 1
%        - Dico: Structure containing the information of the dictionary Dico.dat and Dico.par
%
% Output: - DataOut: N x 4 containing the estimates in this order: R, Vf, Dk, r2
%
%  Nicolas Pannetier, Fev 28th 2013, UCSF


if numel(ADC) ~= size(S,1),
    fprintf('qBoldL2DicoCompute Erro: Data and ADC size mismatch');
end


% Match physio label
% DH2Oix = find(strcmp(Dico.Label,'Model.phy.DH2O'));
% Rix = find(strcmp(Dico.Label,'Model.geo.vasc.R'));
% Vfix = find(strcmp(Dico.Label,'Model.geo.vasc.Vf'));
% Dkix = find(strcmp(Dico.Label,'Model.phy.vasc.khi'));
% B0thetaix = find(strcmp(Dico.Label,'Model.vox.B0theta'));
DH2Oix = strcmp(Dico.Label,'Model.phy.DH2O');
Rix = strcmp(Dico.Label,'Model.geo.vasc.R');
Vfix = strcmp(Dico.Label,'Model.geo.vasc.Vf');
Dkix = strcmp(Dico.Label,'Model.phy.vasc.khi');
B0thetaix = strcmp(Dico.Label,'Model.vox.B0theta');

% Initialization
ADCval  = unique(Dico.par(:,DH2Oix))*1e12;
ADCvalD = diff(ADCval)/2;
ADCvalD1= [Inf; ADCvalD];
ADCvalD2= [ADCvalD; Inf];
ticI    = tic;tic_dt=tic;str = '';
ADCix   = zeros([1 numel(ADC(:))]);
DataOut = zeros([numel(ADC) 5]);


% Parse ADC value of the table to speed up estimation latter

NbLn = sum(Dico.par(:,DH2Oix)*1e12 == ADCval(1));
Idx  = zeros([NbLn numel(ADCval)]);
for a=1:numel(ADCval)
    MapDH2O(a).Idx = find(Dico.par(:,DH2Oix)*1e12 == ADCval(a));
    MapDH2O(a).Val = ADCval(a);
end



% Round ADC map to the closest the ADC value of the table
parfor a=1:numel(ADC(:))
    ADCix(a) = find(ADC(a) > ADCval-ADCvalD1 & ADC(a) <= ADCval+ADCvalD2);
end

% Compute chi2 and estimate parameters from thedictionary
ticI = tic;
tic_dt=tic;
SInd = find(~isnan(S(:,1)) & S(:,1)>0 & ADCix(:) ~= 0); % Out of the mask (may be adapted)
DataOut(:,5) = 2;

nbr=0;
for a=1:numel(SInd) 
    if ~mod(a,20),[str, tic_dt] = DisplayProgress(a,size(SInd,1),str,tic_dt,ticI);end    

    R2 = 1 - sum(bsxfun(@minus,Dico.dat(MapDH2O(ADCix(SInd(a))).Idx,:),S(SInd(a),:)).^2,2)/sum((S(SInd(a),:) - mean(S(SInd(a),:))).^2);
    [R2m,ind] = max(R2(:));
    %     [R2sorted, indices] = sort(R2, 'descend');
        Jdx             = MapDH2O(ADCix(SInd(a))).Idx(ind);

     DataOut(SInd(a),1)    = Dico.par(Jdx,Vfix); % Vf
    DataOut(SInd(a),2)    = Dico.par(Jdx,Rix); % Radius
   

    DataOut(SInd(a),3)    = Dico.par(Jdx,Dkix); %Dk
    DataOut(SInd(a),4)    = Dico.par(Jdx,B0thetaix);  % B0theta
    DataOut(SInd(a),5)    = R2m;
%     if R2m > 0.8
%         
%         figure;plot(S(SInd(a),:), 'Color', 'r'); hold on
%         plot(smooth(S(SInd(a),:)), 'Color', 'g');
%          plot(smooth(S(SInd(a),:), 'lowess'), 'Color', 'b');
%         
%     end
    
    
    loop = 0;
    if R2m>0.8
        nbr = nbr+1;
        if loop == 1
             
            clear adc r bv so2 R20
            [tt, indices] = sort(R2, 'descend');
            indices             = MapDH2O(ADCix(SInd(a))).Idx(indices);
            
                    if abs(Dico.par(indices(1),1)*1e6 - Dico.par(indices(2),1)*1e6) >10 || ...
                             abs(Dico.par(indices(1),2)*100 - Dico.par(indices(2),2)*100) > 10000
                         if loop == 1
            figure;plot(aaa', 'Color', 'r'); hold on
            Dk0 = 0.264 * 4 * pi * 1e-6; Hct=0.375;
            color_curves = [0, 1, 0; 0,0,1];
            leg = {'SI'};
            for i = 1:1
                plot(Dico.dat(indices(i),:)); hold on;%,'Color',color_curves(i,:)); hold on
                adc(i) = Dico.par(indices(i),3)*1e12;
                r(i) = Dico.par(indices(i),1)*1e6;
                bv(i) = Dico.par(indices(i),2)*100;
                so2(i) = (1- Dico.par(indices(i),4)/(Dk0 * Hct))*100;
                R20(i) = tt(i);
%                 leg(size(leg,1)+1,1) ={sprintf('R²=%f\nR=%f\nVf=%f\nStO2=%f\n',R20(i), r(i),bv(i),so2(i) )};
            end
            figure('Name', 'adc');plot(adc)
            figure('Name', 'r');plot(r);
            figure('Name', 'bv');plot(bv);
            figure('Name', 'so2');plot(so2)
            figure('Name', 'R2');plot(R20);
%             legend(leg)
            drawnow
                         end
                    end
        end
    end
    
    
    
    %
    % %         if R2m > 0.8 && (1-  Dico.par(Jdx,Dkix)/(Dk0 * Hct))*100 < 0
    %              close all
    %            plot(S(SInd(a),:),'.'),hold on,%plot(Dico.dat(ind,:),'r-')
    %            plot(S(SInd(a),:));hold on;
    
    %         Dk0 = 0.27 * 4 * pi * 1e-6;
    %         Hct = 0.4;
    %             plot(Dico.dat(Jdx,:));hold on;
    %             leg = sprintf('R²=%f\nR=%f\nVf=%f\nStO2=%f\n',R2m, Dico.par(Jdx,Rix)*1e6,Dico.par(Jdx,Vfix)*100,(1-  Dico.par(Jdx,Dkix)/(Dk0 * Hct))*100);
    
    %              legend(leg)
    
    %             drawnow
    
    
    %         end
    %
end
%
%         tt = DataOut(:,1)*100 ;
%         tt(tt ==0) = [];
%          figure;plot(tt)
[~, ~] = DisplayProgress(a,size(SInd,1),str,tic_dt,ticI);

