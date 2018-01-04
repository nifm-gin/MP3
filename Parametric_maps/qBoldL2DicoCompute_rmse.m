function DataOut = qBoldL2DicoCompute_rmse(S,ADC,Dico)

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


% Initialization
ADCval  = unique(Dico.par(:,4))*1e12;
ADCvalD = diff(ADCval)/2;
ADCvalD1= [Inf; ADCvalD];
ADCvalD2= [ADCvalD; Inf];
ticI    = tic;tic_dt=tic;str = '';
ADCix   = zeros([1 numel(ADC(:))]);

DataOut = NaN([numel(ADC) 4]);
DataOut1 = NaN([numel(ADC) 1]);
DataOut2 = NaN([numel(ADC) 1]);
DataOut3 = NaN([numel(ADC) 1]);
DataOut4 = NaN([numel(ADC) 1]);

% Parse ADC value of the table to speed up estimation latter
NbLn = sum(Dico.par(:,4)*1e12 == ADCval(1));%NbLn = sum(Dico.par(:,5)*1e12 == ADCval(1))
Idx  = zeros([NbLn numel(ADCval)]);
for a=1:numel(ADCval)
    Idx(:,a) = find(Dico.par(:,4)*1e12 == ADCval(a));%Idx(:,a) = find(Dico.par(:,5)*1e12 == ADCval(a));
end

% Round ADC map to the closest the ADC value of the table
ADC(isnan(ADC)) = 0;
parfor a=1:numel(ADC(:))
    ADCix(a) = find(ADC(a) > ADCval-ADCvalD1 & ADC(a) <= ADCval+ADCvalD2);
end

% Compute chi2 and estimate parameters from thedictionary
% ticI = tic;
% tic_dt=tic;
SInd = find(~isnan(S(:,1)) & S(:,1)>0 & ADCix(:) ~= 0); % Out of the mask (may be adapted)
DataOut(:,4) = 2;

% parfor a=1:numel(SInd)
for a=1:size(SInd,1)

    if ~mod(a,20),[str, tic_dt] = DisplayProgress(a,size(SInd,1),str,tic_dt,ticI);end
%    R2 = 1 - sum(bsxfun(@minus,Dico.dat(Idx(:,ADCix(SInd(a))),:),S(SInd(a),:)).^2,2)/sum((S(SInd(a),:) - mean(S(SInd(a),:))).^2);

	rsme = sum(bsxfun(@minus,Dico.dat(Idx(:,ADCix(SInd(a))),:),S(SInd(a),:)).^2,2);

	[R2m,ind] = min(rsme(:));
    
    Jdx             = Idx(ind,ADCix(SInd(a)));
    DataOut1(SInd(a)) =Dico.par(Jdx,1); % Vf
    DataOut2(SInd(a)) =Dico.par(Jdx,2); % Radius
    DataOut3(SInd(a)) = Dico.par(Jdx,3);% Dk
    DataOut4(SInd(a)) = R2m;
   
end
DataOut(:,1) = DataOut1;
DataOut(:,2) = DataOut2;
DataOut(:,3) = DataOut3;
DataOut(:,4) = DataOut4;
% [~, ~] = DisplayProgress(a,size(S,1),str,tic_dt,ticI);

