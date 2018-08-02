function interSliceTime = InterSliceTimeCaslEpi_Adapted(Map, Json)
%% Function for calculate interslice delay for CASL or pCASL labeling module with EPI acquisition module
% written by C.Debacker & Lydiane Hirschler for PV6

freq = Json.ImagingFrequency.value*1e6; % Mz
Rap_Gyr = 2.68e8;
B0 = 2*pi*freq/Rap_Gyr;
B0 = round(B0*10)/10; %Arrondi au 10eme

%B0 = Map.acq.champ_magnetique;
%PVM_MagTransOnOff=scan_acqp('##$PVM_MagTransOnOff=',Map.texte,0);
PVM_MagTransOnOff = Json.MagTransOnOff.value{1};
% PVM_FovSatOnOff=scan_acqp('##$PVM_FovSatOnOff=',Map.texte,0);        % a remettre si saturation remis a sa place par defaut (apres labeling)
PVM_FovSatOnOff = Json.FovSatOnOff.value{1};
% PVM_FatSupOnOff=scan_acqp('##$PVM_FatSupOnOff=',Map.texte,0);
PVM_FatSupOnOff = Json.FatSupOnOff.value{1};
% PVM_TriggerModule=scan_acqp('##$PVM_TriggerModule=',Map.texte,0);
PVM_TriggerModule = Json.TriggerModule.value{1};
% PVM_TriggerOutOnOff=scan_acqp('##$PVM_TriggerOutOnOff=',Map.texte,0);
PVM_TriggerOutOnOff = Json.TriggerOutOnOff.value{1};

% PVM_MagTransModuleTime=scan_acqp('##$PVM_MagTransModuleTime=',Map.texte,1);
PVM_MagTransModuleTime = Json.MagTransModuleTime.value;
% PVM_FovSatModuleTime=scan_acqp('##$PVM_FovSatModuleTime=',Map.texte,1);   % a remettre si saturation remis a sa place par defaut (apres labeling)
PVM_FovSatModuleTime = Json.FovSatModuleTime.value;
% PVM_FatSupModuleTime=scan_acqp('##$PVM_FatSupModuleTime=',Map.texte,1);
PVM_FatSupModuleTime = Json.FatSupModuleTime.value;
% SliceSpoilerDuration=scan_acqp('##$SliceSpoilerDuration=',Map.texte,1); % PAS TROUVE %AJOUTE QUAND MEME
if isfield(Json, 'SliceSpoilerDuration') && ~isempty(Json.SliceSpoilerDuration)
    SliceSpoilerDuration = Json.SliceSpoilerDuration.value;
end
PVM_RiseTime = 0.114;
if B0 >= 9.4 && B0<11.7  % correspond au rise time du 9.4T à Grenoble
    PVM_RiseTime = 0.14656;
end
CFG_AmplifierEnable = 0.05; % 50us
% ExcPulse_Length= scan_acqp('##$ExcPulse=(',Map.texte,1); % PAS TROUVE %AJOUTE QUAND MEME
if isfield(Json, 'ExcPulse') && ~isempty(Json.ExcPulse.value{1})
    ExcPulse_Length = Json.ExcPulse.value{1};
end

PVM_InterGradientWaitTime=0.01;
%RephaseTime=scan_acqp('##$RephaseTime=',Map.texte,1);
RephaseTime = Json.RephaseTime.value;
% SignalType=scan_acqp('##$SignalType=',Map.texte,0);
SignalType = Json.SignalType.value{1};
% EchoTime=scan_acqp('##$EchoTime=',Map.texte,1); % DEJA RECUPERE
EchoTime = Json.EchoTime.value;
% RefPulse_Length= scan_acqp('##$RefPulse=(',Map.texte,1); % PAS TROUVE %AJOUTE QUAND MEME
if isfield(Json, 'RefPulse') && ~isempty(Json.RefPulse.value{1})
    RefPulse_Length = Json.RefPulse.value;
    RfcSpoilerDuration = max(RefPulse_Length, 2*PVM_RiseTime);
else
    RfcSpoilerDuration = 2*PVM_RiseTime;
end

% PVM_RepetitionTime=scan_acqp('##$PVM_RepetitionTime=',Map.texte,1); % DEJA RECUPERE
PVM_RepetitionTime = Json.RepetitionTime.value;
% nSlices=scan_acqp('##$PVM_SPackArrNSlices=( 1 )',Map.texte,1);
nSlices = Json.SPackArrNSlices.value;
% PVM_EpiModuleTime=scan_acqp('##$PVM_EpiModuleTime=',Map.texte,1);
PVM_EpiModuleTime = Json.EpiModuleTime.value;
% PVM_EpiEchoDelay=scan_acqp('##$PVM_EpiEchoDelay=',Map.texte,1);
PVM_EpiEchoDelay = Json.EpiEchoDelay.value;
% PVM_TaggingModuleTime=scan_acqp('##$PVM_TaggingModuleTime=',Map.texte,1);
PVM_TaggingModuleTime = Json.TaggingModuleTime.value;
% ASL_ModuleTime=scan_acqp('ASL_ModuleTime=',Map.texte,1);
ASL_ModuleTime = Json.ASL_ModuleTime.value;
% PackDel=scan_acqp('##$PackDel=',Map.texte,1);
PackDel = Json.PackDel.value;
% ExcPulseTemp = scan_acqp('##$ExcPulse=',Map.texte,2); % PAS TROUVE
if isfield(Json, 'ExcPulse') && ~isempty(Json.ExcPulse.value{1})
    ExcPulseTemp = Json.ExcPulse.value;
end

PvVersion = Json.SoftwareVersions.value{1};
%PvVersion = scan_acqp('##TITLE=Parameter List, ParaVision ',Map.texte,0);
if  ~sum(isnan(PvVersion)) && ~isempty(regexpi(PvVersion,'\w*6.\w*'))
    %SliceSpoiler=strsplit(scan_acqp('##$SliceSpoiler=(',Map.texte,0),', ');
    SliceSpoiler = strsplit(Json.SliceSpoiler.value{1},', ');
    SliceSpoilerDuration = str2double(SliceSpoiler(2));
%     ExcPulse_Length= scan_acqp('##$ExcPul=(',Map.texte,1);
    ExcPulseTemp = strsplit(Json.ExcPul.value{1}, ', ');
    ExcPulse_Length = str2double(ExcPulseTemp{1});
    ExcPulseTemp = strjoin(ExcPulseTemp(1:end-1), ', ');
    %ExcPulse_Length = Json.ExcPul.value{1}; Clement
%     RefPulse_Length= scan_acqp('##$RefPul=(',Map.texte,1);
    RefPulse_Length = strsplit(Json.RefPul.value{1}, ', ');
    RefPulse_Length = str2double(RefPulse_Length{1});
    %RefPulse_Length = Json.RefPul.value{1}; Clement
%     SignalType=scan_acqp('##$PVM_SignalType=',Map.texte,0); % A peu près le meme que plus haut
    SignalType = Json.SignalType.value{1};
%     ASL_ModuleTime=scan_acqp('ASL_ModuleTime=',Map.texte,1); % Le meme que plus haut ...
    ASL_ModuleTime = Json.ASL_ModuleTime.value;
%     ExcPulseTemp = scan_acqp('##$ExcPul=(',Map.texte,0); % Putain c'est quoi cette merde
    %ExcPulseTemp = Json.ExcPul.value{1}; Clement
end                

ExcPulse_RephaseFactorNum = strfind(ExcPulseTemp,',');
ExcPulse_RephaseFactor = str2double(ExcPulseTemp(ExcPulse_RephaseFactorNum(end-1)-3:ExcPulse_RephaseFactorNum(end-1)-1));

if ~strcmp(PVM_MagTransOnOff,'On')
    PVM_MagTransModuleTime = 0;
end
if(~strcmp(PVM_FovSatOnOff,'On'))                 % a remettre si saturation remis a sa place par defaut (apres labeling)
    PVM_FovSatModuleTime = 0;
end
if ~strcmp(PVM_FatSupOnOff,'On')
    PVM_FatSupModuleTime = 0;
end
if strcmp(PVM_TriggerModule,'On')
    trigger = 1;
%     PVM_TriggerMode=scan_acqp('##$PVM_TriggerMode=',Map.texte,0);
    PVM_TriggerMode = Json.TriggerMode.value;
%     PVM_TriggerDelay=scan_acqp('##$PVM_TriggerDelay=',Map.texte,1);
    PVM_TriggerDelay = Json.TriggerDelay.value;
else
    trigger = 0;
    PVM_TriggerMode = 'empty';
end
if strcmp(PVM_TriggerMode,'per_PhaseStep') % /* per volume */
    trigger_v=trigger;
    trigger=0.0;
    PVM_TriggerModuleTime = 0;
elseif  strcmp(PVM_TriggerMode,'per_Slice') 
    trigger_v=0.0;
    %TrigD0 = scan_acqp('##$PVM_TrigD0=',Map.texte,1);
    TrigD0 = Json.TrigD0.value;
    PVM_TriggerModuleTime = PVM_TriggerDelay + 0.04 + TrigD0;
else
    trigger_v=0.0;
    PVM_TriggerModuleTime = 0;
end

if strcmp(PVM_TriggerOutOnOff,'On')
    %TrigOutD0 = scan_acqp('##$PVM_TrigOutD0=',Map.texte,1);
    TrigOutD0 = Json.TrigOutD0.value;
    PVM_TriggerOutModuleTime = 0.01 + TrigOutD0;
    if isnan(PVM_TriggerOutModuleTime)
        %PVM_TriggerOutModuleTime = scan_acqp('##$PVM_TriggerOutModuleTime=',Map.texte,1);
        PVM_TriggerOutModuleTime = Json.TriggerOutModuleTime.value;
    end
    %PVM_TriggerOutMode=scan_acqp('##$PVM_TriggerOutMode=',Map.texte,0);
    PVM_TriggerOutMode = Json.TriggerOutMode.value{1};
    if strcmp(PVM_TriggerOutMode,'PER_SLICE')
        trigOutSlice = PVM_TriggerOutModuleTime;
        trigOutVol = 0.0;
    elseif strcmp(PVM_TriggerOutMode,'PER_VOLUME')
        trigOutSlice = 0.0;
        trigOutVol = PVM_TriggerOutModuleTime;
    elseif strcmp(PVM_TriggerOutMode,'AT_START') 
        % modif Lydiane - pour le cas PVM_TriggerOutMode=AT_START
        % (le trigger est fait une seule fois au début donc ne doit pas intervenir dans le interslice time)
        trigOutSlice = 0.0;
        trigOutVol = 0;
        PVM_TriggerOutModuleTime = 0;
    end
else
    trigOutVol = 0.0;
    trigOutSlice = trigOutVol;
    PVM_TriggerOutModuleTime = trigOutVol;
end

PVM_MinRepetitionTime = ...
    nSlices * ( 0.01 +  ...
    PVM_FatSupModuleTime + ...
    PVM_MagTransModuleTime + ...
    trigger + ...
    trigOutSlice + ...
    SliceSpoilerDuration + ...
    CFG_AmplifierEnable + ...
    ExcPulse_Length/2 + ...
    EchoTime + ...
    PVM_EpiModuleTime - PVM_EpiEchoDelay + ...
    0.01 +  ... % /* 10u */
    0.1  ...  %/* 100u: min d0, includes islice */
    ) +  ...
    PVM_FovSatModuleTime + ...
    PVM_TaggingModuleTime + ...
    ASL_ModuleTime + ...
    trigOutVol + ...
    trigger_v + ...
    PackDel + ...
    0.010;

if( strcmp(SignalType,'FID_signal') )
    PVM_MinEchoTime =  ...    % /* min gradient echo time */
        ExcPulse_Length * ExcPulse_RephaseFactor/100 + ...
        PVM_RiseTime               + ...
        2 * PVM_InterGradientWaitTime  + ...
        RephaseTime       + ...
        PVM_EpiEchoDelay;
    MinTE_left = 0.0; % /* not used */
    MinTE_right  = MinTE_left; % /* not used */
else
    MinTE_left =  ... % /* min half spinecho-time given by left hand side of pi */
        ExcPulse_Length * ExcPulse_RephaseFactor/100  + ...
        PVM_RiseTime               + ...
        2 * PVM_InterGradientWaitTime  +  ...
        RephaseTime       + ...
        RfcSpoilerDuration +  ...
        RefPulse_Length/2.0;
    
    MinTE_right =  ... %/* min half spinecho-time given by right hand side of pi */
        RefPulse_Length/2.0 + ...
        RfcSpoilerDuration + ...
        PVM_InterGradientWaitTime  + ...
        PVM_EpiEchoDelay;
    
    PVM_MinEchoTime = 2 * max(MinTE_left, MinTE_right);
end


if strcmp(SignalType,'SPINECHO_signal')
    Term1 = EchoTime/2 - MinTE_right + PVM_InterGradientWaitTime;
else
    Term1 = EchoTime - PVM_MinEchoTime + PVM_InterGradientWaitTime;
end

if strcmp(SignalType,'SPINECHO_signal') %spin echo
    Term2 = EchoTime/2 - MinTE_left + PVM_InterGradientWaitTime + ...
        RfcSpoilerDuration+ ...
        RefPulse_Length+ ...
        RfcSpoilerDuration;
else
    Term2 = 0;
end

interSliceTime = 0.02 + PVM_TriggerModuleTime + PVM_TriggerOutModuleTime + PVM_MagTransModuleTime + PVM_FatSupModuleTime + SliceSpoilerDuration + ...
    PVM_RiseTime + CFG_AmplifierEnable + ExcPulse_Length + PVM_InterGradientWaitTime + ...
    RephaseTime + Term2 + Term1 + PVM_EpiModuleTime + (PVM_RepetitionTime - PVM_MinRepetitionTime)/nSlices + 0.1;

if isnan(interSliceTime)
    interSliceTime = 50;
    disp('Warning error on interslice time of %s',Map.acq.ppl_name)
end