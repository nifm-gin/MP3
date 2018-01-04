function interSliceTime = InterSliceTimeFairEpi(Map)
%% Function for calculate interslice delay for FAIR labeling module with EPI acquisition module
% wrtitten by C.Debacker & Lydiane Hirschler for PV6

%         PVM_FovSatOnOff=scan_acqp('##$PVM_FovSatOnOff=',Map.texte,0);
PVM_FatSupOnOff=scan_acqp('##$PVM_FatSupOnOff=',Map.texte,0);

PVM_FovSatModuleTime=scan_acqp('##$PVM_FovSatModuleTime=',Map.texte,1);
PVM_FatSupModuleTime=scan_acqp('##$PVM_FatSupModuleTime=',Map.texte,1);
SliceSpoilerDuration=scan_acqp('##$SliceSpoilerDuration=',Map.texte,1);
PVM_RampTime = 0.113;
PVM_RiseTime = 0.114;
PVM_GradDelayTime=0.001; % 1us
CFG_AmplifierEnable = 0.05; % 50us
ExcPulse_Length= scan_acqp('##$ExcPulse=(',Map.texte,1);
PVM_ExSliceRephaseTime=scan_acqp('##$PVM_ExSliceRephaseTime=',Map.texte,1);
PVM_ppgFlag1=scan_acqp('##$PVM_ppgFlag1=',Map.texte,0);
PVM_ppgFlag1=char2logical(PVM_ppgFlag1);
EchoTime=scan_acqp('##$EchoTime=',Map.texte,1);
PVM_InterGradientWaitTime=0.01;
RefPulse_Length= scan_acqp('##$RefPulse=(',Map.texte,1);
SignalType=scan_acqp('##$SignalType=',Map.texte,0);
PVM_MinEchoTime=scan_acqp('##$PVM_MinEchoTime=',Map.texte,1);
PVM_EpiEchoDelay=scan_acqp('##$PVM_EpiEchoDelay=',Map.texte,1);
PVM_EpiModuleTime=scan_acqp('##$PVM_EpiModuleTime=',Map.texte,1);
PVM_TriggerModule=scan_acqp('##$PVM_TriggerModule=',Map.texte,0);

PvVersion = scan_acqp('##TITLE=Parameter List, ParaVision ',Map.texte,0);
if ~isnan(PvVersion(1)) && ~isempty(regexpi(PvVersion,'\w*6.\w*'))
    SliceSpoiler=strsplit(scan_acqp('##$SliceSpoiler=(',Map.texte,0),', ');
    SliceSpoilerDuration = str2double(SliceSpoiler(2));
    ExcPulse_Length= scan_acqp('##$ExcPulse1=(',Map.texte,1);
    PVM_ExSliceRephaseTime=scan_acqp('##$ExSliceRephaseTime=',Map.texte,1);
    RefPulse_Length= scan_acqp('##$RefPulse1=(',Map.texte,1);
    SignalType=scan_acqp('##$PVM_SignalType=',Map.texte,0);
end
if strcmp(PVM_TriggerModule,'On')
    PVM_TriggerMode=scan_acqp('##$PVM_TriggerMode=',Map.texte,0);
    PVM_TriggerDelay=scan_acqp('##$PVM_TriggerDelay=',Map.texte,1);
else
    PVM_TriggerMode = 'empty';
end
if strcmp(PVM_TriggerMode,'empty')
    PVM_TriggerModuleTime = 0;
else
    TrigD0 = scan_acqp('##$PVM_TrigD0=',Map.texte,1);
    PVM_TriggerModuleTime = PVM_TriggerDelay + 0.04 + TrigD0;
    if isnan(PVM_TriggerModuleTime)
       PVM_TriggerModuleTime =  scan_acqp('##$PVM_TriggerModuleTime=',Map.texte,1);
    end
    
end

if( strcmp(SignalType,'FID_signal') )
    MinTE_left = 0.0;
    MinTE_right  = MinTE_left;
else
    MinTE_left = ...                                % /* min half spinecho-time given by left hand side of pi */
        ExcPulse_Length / 2.0 + ...
        PVM_RampTime               + ...
        2 * PVM_InterGradientWaitTime  + ...
        PVM_ExSliceRephaseTime       + ...
        SliceSpoilerDuration + ...
        PVM_RiseTime + ...
        RefPulse_Length/2.0;
    
    MinTE_right = ...                              % /* min half spinecho-time given by right hand side of pi */
        RefPulse_Length/2.0 + ...
        PVM_RampTime + ...
        SliceSpoilerDuration + ...
        PVM_InterGradientWaitTime  + ...
        PVM_EpiEchoDelay;
end

if(~strcmp(PVM_FatSupOnOff,'On'))
    PVM_FatSupModuleTime = 0;
end

if( strcmp(SignalType,'FID_signal') )
    Term1 = EchoTime - PVM_MinEchoTime + PVM_InterGradientWaitTime;
else
    Term1 = EchoTime/2 - MinTE_right + PVM_InterGradientWaitTime;
end

if PVM_ppgFlag1 %spin echo
    Term2 = EchoTime/2 - MinTE_left + PVM_InterGradientWaitTime + ...
        2 * SliceSpoilerDuration + ...
        2 * PVM_RampTime + 2*PVM_GradDelayTime + ...
        RefPulse_Length;
else
    Term2 = 0;
end

interSliceTime = 0.04 + PVM_TriggerModuleTime + PVM_FovSatModuleTime + PVM_FatSupModuleTime + SliceSpoilerDuration + ...
    2 * PVM_RampTime + PVM_GradDelayTime + CFG_AmplifierEnable + ExcPulse_Length + ...
    PVM_ExSliceRephaseTime + Term1 + Term2 + PVM_EpiModuleTime;

if isnan(interSliceTime)
    interSliceTime = 50;
    disp('Warning error on interslice time of %s',Map.acq.ppl_name)
end
