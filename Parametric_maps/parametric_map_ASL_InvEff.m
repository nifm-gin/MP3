function ASL_InvEff = parametric_map_ASL_InvEff(ASL_filename)
% function ASL_InvEff = parametric_map_ASL_InvEff(ASL_filename, add_parameters)
% generate a map of inversion efficienct from a complex ASL_GEFC scan
% coded by C. Debacker

% load the ASL file
fid=fopen(ASL_filename ,'r');
if fid>0
    fclose(fid);
    ASL = load(ASL_filename);
    ASL = ASL.uvascim.image;
else
    warning_text = sprintf('##$ Can not calculate the ASL_InvEff map because there is\n##$ Somthing wrong with the data\n##$ASL=%s\n##$',...
        ASL_filename);
    msgbox(warning_text, 'ASL_InvEff map warning') ;
    ASL_InvEff = [];
    return
end


if(isempty(regexpi(ASL.acq.ppl_name,'\w*casl\w*')) ||...
        (isempty(regexpi(ASL.acq.ppl_name,'\w*gefc\w*'))&&isempty(regexpi(ASL.acq.ppl_name,'\w*FcFlash\w*'))) || ...
        isreal(ASL.reco.data))
    msgbox(sprintf('%s \n \t - aquisition:\t %s \n \t - labeling:\t %s \n \t - reco:\t %s','Scan must combine:', 'GEFC','CASL or pCASL','complex'),'Wrong scan!', 'ASL_InvEff map warning');
    ASL_InvEff = [];
    return
else
    ASL_InvEff=ASL;
   
    % calcul de l'efficacite d'inversion PV6.0 (parties reelles et
    % imaginaires sur images separees)
    
    
    DataControl = ASL.reco.data(:,:,2,:,:);
    DataLabel = ASL.reco.data(:,:,1,:,:);
    alpha = abs( (DataControl - DataLabel) ./ (2 * DataControl) )*100;
    
    
    % Mise en forme donnees
    ASL_InvEff.reco.data = [];
    ASL_InvEff.reco.data(:,:,1,:,:)=alpha;
    
    ParamConfig=sprintf('##$QuantifMethod=ASL_InvEff\n##$ASL=%s\n##$ASL scan info\n%s\n##END=',...
        ASL_filename,[ASL.reco.iminfos{:}]);
    
    ASL_InvEff.reco.paramQuantif = ParamConfig;
    ASL_InvEff.reco.texte = 'ASL_InvEff';
    ASL_InvEff.reco.date = date;
    ASL_InvEff.reco.no_echoes = 1;
    ASL_InvEff.reco.globalmin=min(ASL_InvEff.reco.data(:));
    ASL_InvEff.reco.globalmax=max(ASL_InvEff.reco.data(:));
    ASL_InvEff.reco.fov_offsets = zeros([3,ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
    ASL_InvEff.reco.fov_orientation = zeros([9,ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
    ASL_InvEff.reco.label = {''};
    ASL_InvEff.reco.phaselabel  = {''};
    ASL_InvEff.reco.fov_phase_orientation  = zeros([ASL_InvEff.reco.no_echoes,ASL_InvEff.reco.no_slices,ASL_InvEff.reco.no_expts]);
    ASL_InvEff.reco.scaling_factor = [];
    ASL_InvEff.reco.scaling_offset = [];
    for m_expt=1:ASL_InvEff.reco.no_expts,
        for m_slice=1:ASL_InvEff.reco.no_slices,
            for m_echo=1:ASL_InvEff.reco.no_echoes
                ASL_InvEff.reco.fov_offsets(:,m_echo,m_slice,m_expt) = ASL.reco.fov_offsets(:,1,m_slice,m_expt);
                ASL_InvEff.reco.fov_orientation(:,m_echo,m_slice,m_expt) = ASL.reco.fov_orientation(:,1,m_slice,m_expt);
                ASL_InvEff.reco.label(m_echo,m_slice,m_expt) = ASL.reco.label(1,m_slice,m_expt);
                ASL_InvEff.reco.phaselabel(m_echo,m_slice,m_expt) = ASL.reco.phaselabel(1,m_slice,m_expt);
                ASL_InvEff.reco.fov_phase_orientation(m_echo,m_slice,m_expt) = ASL.reco.fov_phase_orientation(1,m_slice,m_expt);
                ASL_InvEff.reco.scaling_factor(m_echo,m_slice,m_expt) = 1;
                ASL_InvEff.reco.scaling_offset(m_echo,m_slice,m_expt) = 0;
            end
        end
    end
    if isfield(ASL_InvEff.reco, 'fov')
        ASL_InvEff.reco.displayedecho=ASL_InvEff.reco.fov;
    end
    if isfield(ASL_InvEff.reco, 'reco_number')
        ASL_InvEff.reco.reco_number=ASL_InvEff.reco.reco_number;
    else
        ASL_InvEff.reco.reco_number=1;
    end
    if isfield(ASL_InvEff.reco, 'scan_number')
        ASL_InvEff.reco.scan_number=ASL_InvEff.reco.scan_number;
    end
    
    ASL_InvEff.reco.displayedecho=ASL_InvEff.reco.no_echoes;
    ASL_InvEff.reco.displayedslice=ASL_InvEff.reco.no_slices;
    ASL_InvEff.reco.displayedexpt=1;
    ASL_InvEff.reco=orderfields(ASL_InvEff.reco);
    
    ASL_InvEff.clip=[0 100 1];
    
end
