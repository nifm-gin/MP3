function [T2map_a, T2map_b] = parametric_map_T2(MSE_map_filename, json_map_filename, Mask_filename, add_parameters)
%% (xdata,ydata,n,seuil_du_fit)
% generate a T2_mapfrom a RAW MSE scan
% this code come from the SO2_fig function
% additional_parameters correspond to the threshold used and the echo used


seuil_du_fit=str2double(add_parameters{:}(1));
trash_below=str2double(add_parameters{:}(2));
trash_after = str2double(add_parameters{:}(3));
method = add_parameters{:}(4);
remove_last_echo = add_parameters{:}(5) ;
% load the MSME file
fid=fopen(MSE_map_filename ,'r');
if fid>0
    fclose(fid);
    data = load(MSE_map_filename);
    reco = data.uvascim.image.reco;
else
    warning_text = sprintf('##$ Can not calculate the T2 map because there is\n##$ Somthing wrong with the data\n##$MSE=%s\n##$',...
        MSE_map_filename);
    msgbox(warning_text, 'T2 map warning') ;
    T2map_a = [];
    T2map_b = [];
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
T2map_a.acq=data.uvascim.image.acq;
T2map_a.filename=data.uvascim.image.filename;

% Empty memory
clear data

%  human data
if strcmp(remove_last_echo, 'Yes')
    reco.data= reco.data(:,:,1:end-1,:);
    reco.echo_label = reco.echo_label(1:end-1);
    reco.echotime = reco.echotime(1:end-1);
    reco.fov_offsets = reco.fov_offsets(:,1:end-1,:);
    reco.fov_orientation = reco.fov_orientation(:,1:end-1,:);
    reco.fov_phase_orientation = reco.fov_phase_orientation(1:end-1,:);
    reco.label = reco.label(1:end-1,:);
    reco.no_echoes = reco.no_echoes-1;
    reco.phaselabel = reco.phaselabel(1:end-1,:);
    reco.scaling_factor = reco.scaling_factor(1:end-1,:);
    reco.scaling_offset = reco.scaling_offset(1:end-1,:);
    reco.unit = reco.unit(1:end-1,1:end-1);
end

% initializing variables

echotime=reco.echotime(:);
% Selectet echoes used
if trash_after == Inf
    lastecho = length(echotime);
else
    lastecho = sum(echotime < trash_after);
end

firstecho = sum(echotime < trash_below) + 1;
T2map_a.reco='';
T2map_a.reco.data=NaN(reco.no_samples,reco.no_views,1,reco.no_slices);     %
data_MSME=squeeze(reco.data(:,:,:,:));
maxim=max(data_MSME(:)) * seuil_du_fit / 100;
echotimes_used = echotime(firstecho:lastecho)';

tmp_data = reco.data(:,:,firstecho:lastecho,:);

if ~isempty(Mask)
    tmp2_data = zeros(size(tmp_data));
    [~, ~, E, Z] = size(tmp_data);
    for i = 1:Z
        for j = 1:numel(Mask)
            if abs(reco.fov_offsets(3,1,i) - Mask(j).fov_offsets(3)) < 1e-5
                 tmp2_data(:,:,:,i) = tmp_data(:,:,:,i).*repmat(Mask(j).value,[1 1 E 1]);          
            end
        end
    end
    tmp_data = tmp2_data;
    clear tmp2_data
end

tmp_data = permute(tmp_data, [1 2 4 3]);
data_in_vector = reshape(tmp_data, [size(reco.data,1)*size(reco.data,2)*size(reco.data,4), size(tmp_data,4)]);
fit_result_a = NaN(size(data_in_vector,1),1);
tic
switch method{:}
    case 'Mono expo'    
        parfor voxel_nbr = 1:size(data_in_vector,1)
            
            tempydata=data_in_vector(voxel_nbr,:);
            if tempydata(1)>maxim
                [fit_result_a(voxel_nbr),~,~]=fit_exp(echotimes_used,tempydata,1);
            end
        end
        
    case 'Multi expo'
        NbComp = 100;
        fit_result_b = NaN(size(data_in_vector,1),1);
        T2MeraP = zeros([size(data_in_vector,1) numel(echotimes_used)]);
        T2MeraS = zeros([size(data_in_vector,1) NbComp]);
        T2MeraT = zeros([size(data_in_vector,1) NbComp]);
        parfor voxel_nbr=1:size(data_in_vector,1)
            tmp_voxel_data=data_in_vector(voxel_nbr,:);
            out= MERA_1D(tmp_voxel_data',echotimes_used,0,'none',-1,1,NbComp);
            T2MeraS(voxel_nbr,:) = out.S;
            T2MeraT(voxel_nbr,:) = out.T;
        end
        [~, indice] = max(T2MeraS, [], 2);
        for i=1:numel(indice)
            fit_result_a(i) = T2MeraT(i,indice(i));
        end
        fit_result_b = sum(T2MeraT.*T2MeraS,2)./(sum(T2MeraS,2));
        
    case 'EPG'
        [ix] = csvread('Indexx.csv');
        [iy] = csvread('Indexy.csv');
        % t = [1.260000e-002	2.520000e-002	3.780000e-002	5.040000e-002	6.300000e-002	7.560000e-002	8.820000e-002	1.008000e-001	1.134000e-001	1.260000e-001	1.386000e-001	1.512000e-001];
        % figure,plot(t,y)
        
        esp = [6.3e-3 diff(echotimes_used'*0.001)];
        flipangle = pi * ones([numel(esp) 1]);
        etl = numel(flipangle);
        T1 = 4000e-3;
        maxim=max(data_in_vector(:,1)) * 5 / 100;
        warning('off'); 
        B = NaN(size(data_in_vector,1),3);
        tmpnum = 1:round((size(data_in_vector,1)/1000)):size(data_in_vector,1);
        for aa=1:size(data_in_vector,1)
            y = data_in_vector(aa,:);
            
            if y(1)>maxim && ~isnan(y(1))
%                 if sum(aa == tmpnum) ==1
%                     aa
%                 end
                xdat.esp = esp;
                xdat.flipangle = flipangle;
                xdat.etl = numel(flipangle);
                xdat.T1 = T1;
                xdat.y = y;
                ydat = y;
                
                T2 = 0.120;
                B1 = 0.6;
                x0 = [T2 B1 ydat(1)];
                b_num = numel(x0);
                pars       = struct(...
                    'N',     numel(ydat),...       % number of Y data points
                    'b_num', b_num,...             % number of parameters of the model
                    'A',     zeros(b_num),...      % variance covariance matrix of partial derivatives
                    'G',     zeros(b_num,1),...    % gradient vector to update paramters
                    'D',     zeros(b_num,1),...    % vector of parameter std's
                    'q',     [],...                % ??? (some sort of normalized covariance matrix)
                    'B',     x0(:),...         % current parameter guess
                    'INC',   1+eye(b_num)*0.001,...% perturbation of B to estimate partial derivatives
                    'L',     0.1*eye(b_num),...    % lambda
                    'dL',    10,...                % factorial increment or decrement of lambda
                    'CONV',  0.001,...             % convergence criterion
                    'T1',    30,...                % limit of number of iterations (and output status)
                    'T2',    0,...                 % iteration counter
                    'F',     @evalfunction,...                 % model function
                    'E',     [],...                % error vector for current parameter estimate
                    'SO',    [],...                % sum of squares error of current parameter guess
                    'verbose',0);                  % flag indicating whether progress messages should be printed
                
                
                
                [B(aa,:), ~, STATUS]=levenbergmarquardt(pars,xdat,ydat,x0);
                %%% code to plot both fit (with and without taking in
                %%% account the stimulated echoes
                %                 s =  evalfunction(xdat,B(aa,:));
                %                 sini =  evalfunction(xdat,x0);
                %                 [t2,cte,~]=fit_exp(echotimes_used,xdat.y',1);
                %                 x2=cte *  exp(-echotimes_used/t2);
                %                 figure(10);
                %                 plot(abs(s)),hold on,plot(x2,'x'),plot(xdat.y,'.')
                %                 legend({[num2str(B(aa,1)*1000) '-echostim'] num2str(t2)}, 'Location','NorthEast');
                %                 hold off
            end
        end
         fit_result_a =B(:,1)*1000;
        %         B1map = nan([128 128]);
        %         T2map = nan([128 128]);
        %         %%%%%%%%
%         tmp = B;
%         for iii=12088:128*128
%         tmp(iii,1:3) = [0 0 0];
%         end
%         tmp_a=reshape(tmp,[128,128,3]);
%         T2map_a.reco.data=permute(tmp_a, [1 2 4 3]);
%         figure,imagesc(tmp_a(:,:,1));
%          figure,imagesc(tmp_a(:,:,2));
%          figure,imagesc(tmp_a(:,:,3));
        
        %%%%%%%
             
        
%         
%         
%         for a=1:numel(ix)
%             T2map(ix(a),iy(a)) = B(a,1);
%             B1map(ix(a),iy(a)) = B(a,2);
%         end
%         figure,imagesc(T2map)
%         figure,imagesc(B1map)
%         
%         
%         T2guess = T2FitN(data_in_vector,echotime,'dim',2);
%         
%         sprintf('T2 exp = %.4f +/- %.4f',mean(T2guess(:,2)),std(T2guess(:,2)))
%         sprintf('T2 EPG = %.4f +/- %.4f',mean(B(:,1)),std(B(:,1)))
%         sprintf('B1 EPG = %.4f +/- %.4f',mean(B(:,2)),std(B(:,2)))
end
toc

% Complete the T2star_map.reco structure
T2map_a.reco = reco;
% update the structure 
T2map_a.reco.echo_label ='';
T2map_a.reco.echotime = [];
T2map_a.reco.fov_offsets = reco.fov_offsets(:,1,:);
T2map_a.reco.fov_orientation = reco.fov_orientation(:,1,:);
T2map_a.reco.fov_phase_orientation = T2map_a.reco.fov_phase_orientation(1,:);
T2map_a.reco.label = T2map_a.reco.label (1,:);
T2map_a.reco.no_echoes = 1;
T2map_a.reco.phaselabel = T2map_a.reco.phaselabel(1,:);
T2map_a.reco.scaling_factor = T2map_a.reco.scaling_factor(1,:);
T2map_a.reco.scaling_offset =  T2map_a.reco.scaling_offset(1,:);
T2map_a.reco.unit = 'ms';
T2map_a.reco.mask_to_use='';
T2map_a.reco.texte='T2map';
T2map_a.reco.unit='ms';
T2map_a.reco.date=date;

ParamConfig=sprintf('##$QuantifMethod=levenbergmarquardt\n##$Fit threshold=%s\n#$Echo trash below (ms)=%s\n##$Echo trash after (ms)=%s\n##$Raw scan used=%s\n##$Method used=%s\n',...
    add_parameters{:}{1},...
    add_parameters{:}{2},...
    add_parameters{:}{3},...
    add_parameters{:}{4},...
    MSE_map_filename);
T2map_a.reco.paramQuantif = ParamConfig;

tmp_a=reshape(fit_result_a,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
T2map_a.reco.data=permute(tmp_a, [1 2 4 3]);
T2map_a.reco.globalmax=max(T2map_a.reco.data(:));
T2map_a.reco.globalmin=min(T2map_a.reco.data(:));

T2map_a.scan_number = 104;
T2map_a.clip = [0 200 1];
T2map_a.reco_number = 1;
T2map_a.reco=orderfields(T2map_a.reco);
switch method{:}
    case 'Mono expo'
        T2map_b = [];
    case 'Multi expo'
        T2map_b = T2map_a;
        T2map_b.reco = rmfield(T2map_b.reco, 'data');
        tmp_b=reshape(fit_result_b,[size(reco.data,1),size(reco.data,2),size(reco.data,4)]);
        T2map_b.reco.data=permute(tmp_b, [1 2 4 3]);
        T2map_b.reco.globalmax=max(T2map_b.reco.data(:));
        T2map_b.reco.globalmin=min(T2map_b.reco.data(:));
        
        T2map_b.reco=orderfields(T2map_b.reco);
    case 'EPG'
        T2map_b = [];    
end




