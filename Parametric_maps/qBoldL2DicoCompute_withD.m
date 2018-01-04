function DataOut = qBoldL2DicoCompute_withD(S,Dico)

% Process the signal by comparing the input S to the dictionary of
% curves (Dico)
% Input: - S: Array of signals. N x T where length of T must the same than
% size(Dico.dat,2)
%        - Dico: Structure containing the information of the dictionary Dico.dat and Dico.par
%
% Output: - DataOut: N x 6 containing the estimates in this order: R, Vf, Dk, B0theta, DH20, r2
%
%  Nicolas Pannetier, Fev 28th 2013, UCSF


% if numel(ADC) ~= size(S,1),
%     fprintf('qBoldL2DicoCompute Erro: Data and ADC size mismatch');
% end


% Match physio label
DH2Oix = strcmp(Dico.Label,'Model.phy.DH2O');
Rix = strcmp(Dico.Label,'Model.geo.vasc.R');
Vfix = strcmp(Dico.Label,'Model.geo.vasc.Vf');
Dkix = strcmp(Dico.Label,'Model.phy.vasc.khi');
B0thetaix = strcmp(Dico.Label,'Model.vox.B0theta');


% Initialization
ticI    = tic;tic_dt=tic;str = '';



% Compute chi2 and estimate parameters from thedictionary
ticI = tic;
tic_dt=tic;
SInd = find(~isnan(S(:,1)) & S(:,1)>0); % Out of the mask (may be adapted)

R2m = zeros([size(SInd,1) 1]);
ind =zeros([size(SInd,1) 1]);

dico_data = Dico.dat;
sub_data = S(SInd,:);


DataOut = zeros([size(S,1) 6]);

tic
for a=1:numel(SInd)
    
    % loop for
                if ~mod(a,100),[str, tic_dt] = DisplayProgress(a,size(SInd,1),str,tic_dt,ticI);end
    
                 R2 = 1 - sum(bsxfun(@minus,Dico.dat,S(SInd(a),:)).^2,2)/sum((S(SInd(a),:) - mean(S(SInd(a),:))).^2);
                [R2m,ind] = max(R2(:));
    
                Jdx             = ind;
    
                DataOut(SInd(a),1)    =  Dico.par(Jdx,Vfix); % Vf
                DataOut(SInd(a),2)    = Dico.par(Jdx,Rix); % Radius
                DataOut(SInd(a),3)    = Dico.par(Jdx,Dkix); % Dk
                DataOut(SInd(a),4)    = Dico.par(Jdx,B0thetaix); % B0theta
                DataOut(SInd(a),5)    = Dico.par(Jdx,DH2Oix); % DH2O
                DataOut(SInd(a),6)    = R2m;
    
    
%     %% loop parfor
%     R2 = 1 - sum(bsxfun(@minus,dico_data,sub_data(a,:)).^2,2)/sum((sub_data(a,:) - mean(sub_data(a,:))).^2);
%     [R2m(a),ind(a)] = max(R2(:));
    
    %     R2 = 1 - sum(bsxfun(@minus,Dico.dat(MapDH2O(ADCix(SInd(a))).Idx,:),S(SInd(a),:)).^2,2)/sum((S(SInd(a),:) - mean(S(SInd(a),:))).^2);
    %         [R2m,ind] = max(R2(:));
    %
    %         Jdx             = MapDH2O(ADCix(SInd(a))).Idx(ind);
    
%     
%     tt = R2;
%     tt(:,2) = 1:size(R2,1)';
%      % sort the data manualy in matlab
%      figure;plot(sub_data(a,:)); hold on
%       Dk0 = 0.264 * 4 * pi * 1e-6; Hct=0.375;
%      % bv*100
%     for i = 1:100
%         plot(dico_data(tt(i,2),:)); hold on
%         adc(i) = Dico.par(tt(i,2),3)*1e12;
%         r(i) = Dico.par(tt(i,2),1)*1e6;
%         bv(i) = Dico.par(tt(i,2),2)*100;
%         so2(i) = (1- Dico.par(tt(i,2),4)/(Dk0 * Hct))*100;
%     end
%     figure('Name','adc');plot(adc)
%     figure('Name','r');plot(r)
%     figure('Name','bv');plot(bv)
%     figure('Name','so2');plot(so2)
end
toc

% DataOut = zeros([size(S,1) 6]);
% DataOut(:,6) = 2;
% for a=1:numel(SInd)
%     DataOut(SInd(a),1)  =  Dico.par(ind(a),Vfix); % Vf
%     DataOut(SInd(a),2)  = Dico.par(ind(a),Rix); % Radius
%     DataOut(SInd(a),3)   = Dico.par(ind(a),Dkix); % Dk
%     DataOut(SInd(a),4)   = Dico.par(ind(a),B0thetaix); % B0theta
%     DataOut(SInd(a),5)    = Dico.par(ind(a),DH2Oix)*1e12; % DH2O
%     DataOut(SInd(a),6)    =  R2m(a);
% end

% [~, ~] = DisplayProgress(a,size(S,1),str,tic_dt,ticI);

