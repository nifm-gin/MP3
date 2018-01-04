function [output1, output2,output3]= ASL_HadamardFit(it_sl,ASL,VoiOutlier_filename,M0,SliceDecM0,T1app,PVM_ObjOrderList,interSliceTime,alpha,T1sang,lambda,progression_Hadamard,fit_CBF_result,fit_transit_result,fit_T1app_result,ASL_filename)
fitT1 = 1;
%if progression_Hadamard == 0
%% Init variables
[Voi, VoiYn] = extract_map(VoiOutlier_filename,'ASL-dyn','VOI');
if isempty(VoiOutlier_filename)
    VoiYn = 0;
end

if VoiYn
    factorResize = size(ASL.reco.data,1) / size(Voi(it_sl).value,1);
    VoiTemp = imresize(Voi(it_sl).value,factorResize,'nearest');
    ASL.reco.data(:,:,:,it_sl,:) = ASL.reco.data(:,:,:,it_sl,:).*repmat(VoiTemp,[1 1 size(ASL.reco.data,3) 1 size(ASL.reco.data,5)]);
    if numel(M0) ~= 1
        M0= M0.*VoiTemp;
    end
    if numel(T1app) ~= 1
        T1app= T1app.*VoiTemp;
    end
end

data_in_vector = reshape(ASL.reco.data, [size(ASL.reco.data,1)*size(ASL.reco.data,2),...
    size(ASL.reco.data,3), size(ASL.reco.data,4),size(ASL.reco.data,5)]);
m0_in_vector = reshape(M0, [size(M0,1)*size(M0,2),...
    size(M0,3), size(M0,4),size(M0,5)]);
if numel(T1app) ~= 1
    t1_in_vector = reshape(T1app, [size(T1app,1)*size(T1app,2),...
        size(T1app,3), size(T1app,4),size(T1app,5)]);
end
sliceDecM0_in_vector = reshape(SliceDecM0, [size(SliceDecM0,1)*size(SliceDecM0,2),...
    size(SliceDecM0,3), size(SliceDecM0,4),size(SliceDecM0,5)]);


%% Retrieve parameter values
PostLabelTime=scan_acqp('ASL_PostLabelTime=',ASL.texte,1)/1000;
Subboli_Duration = scan_acqp('ASL_SubboliTime=',ASL.texte,1)/1000;
TotalLabeling_Duration = scan_acqp('ASL_LabelTime=',ASL.texte,1)/1000;
Effective_PLD = [];
for subbolus = 1:size(ASL.reco.data,3)
    Effective_PLD(subbolus) = TotalLabeling_Duration - sum(Subboli_Duration(1:subbolus)) + PostLabelTime ;
end

NumSlice = find(PVM_ObjOrderList==it_sl-1)-1;
if size(PVM_ObjOrderList) == 1 % cas 1 slice ou 3D
    NumSlice=0;
end
xfix.Effective_PLD=Effective_PLD(end:-1:1)+ interSliceTime/1000 * NumSlice;
step =  Subboli_Duration(1);
xfix.alpha = alpha;
xfix.TotalLabeling_Duration =TotalLabeling_Duration;
xfix.step = step;
xfix.T1sang  = T1sang/1000;
xfix.PostLabelTime  = PostLabelTime ;
xfix.SubboliDuration = Subboli_Duration;
xfix.model = 'Buxton';
CBFi = 0.05;
N = size(data_in_vector,1);
parfor_progress(N); % Initialize
for voxel_nbr=1:N
    %         parfor_progress; % Count
    
    if sum(isnan(data_in_vector(voxel_nbr,:,it_sl))) ~= length(data_in_vector(voxel_nbr,:,it_sl)) && sum(data_in_vector(voxel_nbr,:,it_sl) == 0) ~= length(data_in_vector(voxel_nbr,:,it_sl)) & ...
            m0_in_vector(voxel_nbr,:) ~= 0 && t1_in_vector(voxel_nbr,:) ~= 0
        
        
        % Conditions initiales fit
        %             if numel(T1app) ~= 1
        xfix.T1app = t1_in_vector(voxel_nbr)/1000;
        xfix.SliceDecM0 = sliceDecM0_in_vector(voxel_nbr);
        %xfix.SliceDecM0 = 1;
        T1init = t1_in_vector(voxel_nbr)/1000;
        %             else
        %                 xfix.T1app = T1app/1000;
        %                 T1init =T1app/1000;
        %                 xfix.SliceDecM0 = sliceDecM0_in_vector;
        %             end
        
        %             if numel(M0) ~= 1
        xfix.M0 =  m0_in_vector(voxel_nbr,:);
        %             else
        %                 xfix.M0 = M0 ;
        %             end
        
        
        [maxData, maxTime] = max(data_in_vector(voxel_nbr,end:-1:3,it_sl)); %temps de transit
        Ttransiti = max((maxTime-2)*step+xfix.Effective_PLD(1),0.1); %maxy-2 car le max est a Ttransit+LabelingTime et l'indice de maxy commence a 1 et pas 0
        %CBFi = maxData*lambda / (2*  xfix.M0 * xfix.SliceDecM0 /lambda * T1init * alpha * exp(-Ttransiti/(T1sang/1000)));
        
        % FIT
        if fitT1 == 1
            %Ttransiti=Ttransitinit(it_sl);
            [aaa, bbb,  ccc]=levenbergmarquardt('HadamardFitFunction',xfix,data_in_vector(voxel_nbr,end:-1:1,it_sl),[Ttransiti CBFi 0.2]);
            maxCBF = HadamardFitFunctionMax(xfix,aaa);
        else
            [aaa, bbb,  ccc]=levenbergmarquardt('HadamardFitFunction',xfix,data_in_vector(voxel_nbr,end:-1:1,it_sl),[Ttransiti CBFi]);
        end
        %         if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 & ccc==-1 %#ok<AND2>
        %fit_CBF_result(voxel_nbr,it_sl)=aaa(2);
        %fit_CBF_result(voxel_nbr,it_sl)=maxCBF/(2.*xfix.M0.*xfix.SliceDecM0/lambda*alpha.* xfix.T1app *exp(-aaa(1)./xfix.T1sang).*(1-exp(-Subboli_Duration(1)./ xfix.T1app)));
        fit_CBF_result(voxel_nbr,it_sl)=maxCBF/(2.*xfix.M0.*xfix.SliceDecM0/lambda*alpha.* xfix.T1app *exp(-aaa(1)./xfix.T1sang).*(1-exp(-Subboli_Duration(1)./ xfix.T1app)));
        %fit_CBF_result(voxel_nbr,it_sl)=maxCBF;
        %fit_CBF_err(voxel_nbr,it_sl)=bbb(2);
        fit_transit_result(voxel_nbr,it_sl)=aaa(1); %+ interSliceTime/1000 * NumSlice;
        %fit_transit_err(voxel_nbr,it_sl)=bbb(1);
        if fitT1 == 1
            fit_T1app_result(voxel_nbr,it_sl)=aaa(3);
            % fit_T1app_err(voxel_nbr,it_sl)=bbb(3);
        end
        %        end
    end
end
parfor_progress(0); % Clean up
output1 = fit_CBF_result(:,it_sl);
output2 = fit_transit_result(:,it_sl);
if fitT1 == 1
    output3 = fit_T1app_result(:,it_sl);
end

