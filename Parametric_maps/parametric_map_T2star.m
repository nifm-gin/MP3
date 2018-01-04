function T2star_map = parametric_map_T2star(MGE_map_filename, Mask_filename, add_parameters)
% generate a T2star_map.reco from a RAW MultiGradientEchoes scan
% this code come from the VSI_fit_t2s function
% additional_parameters correspond to the method of fit used and the
% threshold used

mask_to_use='';
seuil_du_fit=str2double(add_parameters{:}(1));
fitmethod=add_parameters{:}(2);
sequence_used = add_parameters{:}(3);
trash_below=str2double(add_parameters{:}(4));
trash_after = str2double(add_parameters{:}(5));

% load the MGE file
fid=fopen(MGE_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(MGE_map_filename);
    reco = data.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the T2star map because there is\n##$ Somthing wrong with the data\n##$T2*=%s\n##$',...
        MGE_map_filename);
    msgbox(warning_text, 'T2star map warning') ;
    T2star_map = [];
    return
end

% Load mask (if selected and/or exist)
if ~isempty(Mask_filename)
    fid=fopen(Mask_filename ,'r');
    if fid>0
        fclose(fid);
        Mask = load(Mask_filename);
        Mask = Mask.uvascroi;
    else
        Mask = []; 
    end
else
    Mask = [];
end


% save imformation
T2star_map=data.uvascim.image;
T2star_map.reco_number = 1;
T2star_map.filename=data.uvascim.image.filename;


% Empty memory
clear data

% for later
mask_to_use(1:reco.no_samples,1:reco.no_views)=1;

% initializing variables
% T2star_map.reco='';
T2star_map.reco.data=NaN(reco.no_samples,reco.no_views,2,reco.no_slices);     % parameters(T2, M0)
T2star_map.reco.err=zeros(reco.no_samples,reco.no_views,2,reco.no_slices);     % uncertainties on the parameters
T2star_map.reco.wrongpix=zeros(reco.no_samples,reco.no_views,reco.no_slices);  % map of excluded pixels
T2star_map.reco.ki2=zeros(reco.no_samples,reco.no_views,reco.no_slices);       % Adjustment quality
echotime=reco.echotime(:);                                                % echo time in ms

% if MGESE or MGEFIDSE
switch sequence_used{:}
    case {'MGEFIDSE', 'MGESE'}
        SpinEcho = scan_acqp('##$SpinEchoTime=',T2star_map.texte,1);
        if isnan(SpinEcho)
             SpinEcho = scan_acqp('##$PVM_EchoTime=',T2star_map.texte,1);
        end
        lastecho = sum(echotime<SpinEcho/2);
        if trash_after ~= Inf
            lastecho = sum(echotime < trash_after);
        end
    otherwise
        if trash_after == Inf
            lastecho = length(echotime);
        else
            lastecho = sum(echotime < trash_after);
        end
end
% Human data
% lastecho = 9;

firstecho = sum(echotime < trash_below) + 1;


tmp_data = reco.data;
if ~isempty(Mask)
    [~, ~, E, Z] = size(tmp_data);
    for i = 1:Z
        for j = 1:numel(Mask)
            if abs(reco.fov_offsets(3,1,i) - Mask(j).fov_offsets(3)) < 1e-5
                 tmp_data(:,:,:,i) = tmp_data(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);          
            end
        end
    end
end
tmp_data = permute(tmp_data, [1 2 4 3]);
data_in_vector = reshape(tmp_data, [size(reco.data,1)*size(reco.data,2)*size(reco.data,4), size(reco.data,3)]);
maxim=max(data_in_vector(:,1)) * seuil_du_fit / 100;
first_echotime=echotime(firstecho);
last_echotime= echotime(lastecho);
echotime_used = echotime(firstecho:lastecho);
fit_result = NaN(size(data_in_vector,1),1);
fit_err= NaN(size(data_in_vector,1),1);
wrongpix= NaN(size(data_in_vector,1),1);
tic
% %% smooth data
% for i=1:size(data_in_vector,1)
%     if sum(data_in_vector(i,:)) ~= 0
%         data_in_vector(i,:) = smooth(data_in_vector(i,:), 'lowess')';
%     end
% end

% y=AB_t2s(1,[1 1]); % used in order to do that AB_t2s exists when compiling matlab
parfor voxel_nbr=1:size(data_in_vector,1)
    tempydata=data_in_vector(voxel_nbr,:);
    if tempydata(1)>maxim
        % fit initializing
        t2s=(first_echotime-last_echotime)/log(tempydata(lastecho)/tempydata(firstecho));
        if t2s<=0 || isnan(t2s),
            t2s=30;
        end
        % apply the fit
        [aaa, bbb,  ~]=levenbergmarquardt('AB_t2s',echotime_used, tempydata(firstecho:lastecho)',[t2s max(tempydata(firstecho:lastecho))*1.5]);
        %%%%%%%%%%%%%%%%
        %         out= MERA_1D( tempydata(firstecho:lastecho)',echotime_used,0,'none',-1,1,100);
        %         T2MeraS(voxel_nbr,:) = out.S;
        %         T2MeraT(voxel_nbr,:) = out.T;
        %         aaa(1) = sum(out.T.*out.S)./sum(out.S);
        %         aaa(2)= sum(out.T.*out.S)./sum(out.S);
        %         bbb =0;
        %         voxel_nbr
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 %#ok<AND2>
            fit_result(voxel_nbr)=aaa(1);
            fit_err(voxel_nbr)=bbb(1);
        else
            wrongpix(voxel_nbr)=2; % the fit does not work
        end
    else % else below the fit threshold
        wrongpix(voxel_nbr)=1;
    end
end
tmp=reshape(fit_result,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
T2star_map.reco.data=permute(tmp, [1 2 4 3]);
tmp=reshape(fit_err,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
T2star_map.reco.err=permute(tmp, [1 2 4 3]);
tmp=reshape(wrongpix,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
T2star_map.reco.wrongpix=permute(tmp, [1 2 4 3]); % the fit does not work
toc

% echoindex=1:reco.no_echoes;
%
% %Adjustment loop
% h = waitbar(0, 'T2* map in progress...');
% steps_tot= reco.no_slices*reco.no_samples*reco.no_views;
% current_step = 1;
% tic
% for m_slice=1:reco.no_slices,
%     tempim=squeeze(reco.data(:,:,1,m_slice));
%     maxim=max(tempim(:)) * seuil_du_fit / 100;
%     index=find(tempim > maxim);
%     nbp=length(index);
%     curnbp=0;
%
%     tempydata=squeeze(reco.data(:,:,echoindex,m_slice));
%     y=AB_t2s(1,[1 1]); % used in order to do that AB_t2s exists when compiling matlab
%     for m=1:reco.no_samples,
%         for n=1:reco.no_views,
%             waitbar(current_step/steps_tot, h)
%             if (tempim(m,n)*double(mask_to_use(m,n)))>maxim,
%                 % fit initializing
%                 ydata=squeeze(tempydata(m,n,:));
%                 t2s=(echotime(firstecho)-echotime(lastecho))/log(ydata(lastecho)/ydata(firstecho));
%                 if t2s<=0 || isnan(t2s),
%                     ts2=30;
%                 end
%                 % apply the fit
%                 [aaa, bbb,  convergence]=levenbergmarquardt('AB_t2s',echotime(firstecho:lastecho), ydata(firstecho:lastecho),[t2s max(ydata(firstecho:lastecho))*1.5]);
%                 if aaa(1)>0 & aaa(2)>0 & imag(aaa)==0 %#ok<AND2>
%                     T2star_map.reco.data(m,n,1:2,m_slice)=aaa;
%                     T2star_map.reco.err(m,n,1:2,m_slice)=bbb;
%                 else
%                     T2star_map.reco.wrongpix(m,n,m_slice)=2; % the fit does not work
%                 end
%                 curnbp=curnbp+1;
%             else % else below the fit threshold
%                 T2star_map.reco.wrongpix(m,n,m_slice)=1;
%             end
%             current_step = current_step+1;
%         end
%     end
% end
% close(h);
% % Complete the T2star_map.reco structure
% toc
T2star_map.reco.no_echoes=1; %Number of parameters stored here
T2star_map.reco.no_expts=1;

%Adapt the fov offsets and orientation infos
tmpfov=reco.fov_offsets(:,1,:,1);
tmpori=reco.fov_orientation(:,1,:,1);
% tmpphaselabel=reco.phaselabel(1,:,1);
T2star_map.reco.fov_offsets = [];
T2star_map.reco.fov_orientation = [];
T2star_map.reco.scaling_factor = [];
T2star_map.reco.scaling_offset = [];
T2star_map.reco.label = {};
T2star_map.reco.phaselabel = {};

for m_slice=1:reco.no_slices,
    T2star_map.reco.fov_offsets(:,:,m_slice)=squeeze(tmpfov(:,m_slice));
    T2star_map.reco.fov_orientation(:,:,m_slice)=squeeze(tmpori(:,m_slice));
    T2star_map.reco.label(1,m_slice)=reco.label(1,m_slice,1);
    T2star_map.reco.phaselabel(1,m_slice)=reco.phaselabel(1,m_slice,1);
    T2star_map.reco.scaling_factor(1,m_slice) = reco.scaling_factor(1,m_slice,1);
    T2star_map.reco.scaling_offset(1,m_slice) = reco.scaling_offset(1,m_slice,1);
end
T2star_map.reco.fov_phase_orientation=reco.fov_phase_orientation(1,:);
T2star_map.reco.mask_to_use=mask_to_use;
T2star_map.reco.texte='T2*map';
T2star_map.reco.unit='ms';
T2star_map.reco.date=date;

if isfield(reco, 'paramQuantif')
      ParamConfig=sprintf('##$QuantifMethod=%s\n##$Fit threshold=%s\n##$Sequence used=%s\n##$Echo trash below (ms)=%s\n##$Echo trash after (ms)=%s\n##$Raw scan used=%s\n\n##$MGE scan info\n%s',...
        fitmethod{:},...
        add_parameters{:}{1},...
        sequence_used{:},...
        add_parameters{:}{4},...
        add_parameters{:}{5},...
        MGE_map_filename,...
        reco.paramQuantif);
else
    ParamConfig=sprintf('##$QuantifMethod=%s\n##$Fit threshold=%s\n##$Sequence used=%s\n##$Echo trash below (ms)=%s\n##$Echo trash after (ms)=%s\n##$Raw scan used=%s',...
        fitmethod{:},...
        add_parameters{:}{1},...
        sequence_used{:},...
        add_parameters{:}{4},...
        add_parameters{:}{5},...
        MGE_map_filename);
end
T2star_map.reco.paramQuantif = ParamConfig;

T2star_map.reco.data=T2star_map.reco.data(:,:,1,:);
T2star_map.reco.globalmax=max(T2star_map.reco.data(:));
T2star_map.reco.globalmin=min(T2star_map.reco.data(:));
T2star_map.reco.echo_label ='';
T2star_map.reco_number = 1;
T2star_map.reco=orderfields(T2star_map.reco);

%delete useless field
if isfield(T2star_map.reco, 'angAP')
   T2star_map.reco = rmfield(T2star_map.reco, 'angAP'); 
end
if isfield(T2star_map.reco, 'angFH')
   T2star_map.reco = rmfield(T2star_map.reco, 'angFH'); 
end
if isfield(T2star_map.reco, 'angRL')
   T2star_map.reco = rmfield(T2star_map.reco, 'angRL'); 
end
if isfield(T2star_map.reco, 'angulation')
   T2star_map.reco = rmfield(T2star_map.reco, 'angulation'); 
end
if isfield(T2star_map.reco, 'echotime')
   T2star_map.reco = rmfield(T2star_map.reco, 'echotime'); 
end
if isfield(T2star_map.reco, 'iminfos')
   T2star_map.reco = rmfield(T2star_map.reco, 'iminfos'); 
end
if isfield(T2star_map.reco, 'reco_meth')
   T2star_map.reco = rmfield(T2star_map.reco, 'reco_meth'); 
end

% T2star_map.reco = rmfield(T2star_map.reco, 'angFH');
% T2star_map.reco = rmfield(T2star_map.reco, 'angRL');
% T2star_map.reco = rmfield(T2star_map.reco, 'angulation');
% T2star_map.reco = rmfield(T2star_map.reco, 'echotime');
% T2star_map.reco = rmfield(T2star_map.reco, 'iminfos');
% T2star_map.reco = rmfield(T2star_map.reco, 'reco_meth');


T2star_map.scan_number = 104;
T2star_map.clip = [0 100 1];



