folder_list = dir;
folder_list = folder_list(3:end);
pathname = pwd;
for j = 1:numel(folder_list)
    if folder_list(j).isdir
        par_file_listing = dir(strcat(pathname,filesep,folder_list(j).name,filesep, '*.PAR'));
        if isempty(par_file_listing)
            continue
        end
        subfold_pathname = [pathname filesep,folder_list(j).name];
        for i = 1:numel(par_file_listing) 
            %     if strcmp(par_file_listing(i).name, 'H09_v2_morpho_111107_MSME_TR_COURT_SENSE_9_1.PAR')
            if ~isempty(strfind(par_file_listing(i).name, 'vaso'))
                continue
            end
            DirPat= [subfold_pathname filesep par_file_listing(i).name];
            
            options.outputformat = 1;
            options.angulation = 1;%def=1
            options.subaan = 1;
            options.rescale = 1;
            options.pathpar = [subfold_pathname filesep];
            options.usefullprefix = 0;
            options.usealtfolder = 0;
            options.altfolder = '';
            
            [~, options.parname, ~]= fileparts(par_file_listing(i).name);
            
            options.prefix = '';
            options.parname = [options.parname, '.PAR'];
            % convert_r2a_bis({options.parname},options);
            % convert_r2a({options.parname},options);
            GFB_convert_r2a({options.parname},options);
            
            % Move par/Rec file into the same directory as the nii(s) file(s)
            [pathstr, name, ext]=fileparts(options.parname);
            parrec = GFB_convert_read_parrec([subfold_pathname, filesep,  name, '.par']);
            fNameRoot = [parrec.hdr.scn.protocol_name];
            % clean up file name somewhat
            fNameRoot = strrep(fNameRoot,'WIP','');
            fNameRoot = strrep(fNameRoot,'SENSE','');
            fNameRoot = strtrim(fNameRoot);
            fNameRoot = strrep(fNameRoot,' ','_');
            %move the par/rec files into the nii files
            movefile([subfold_pathname, filesep,  name, '.PAR'], [subfold_pathname, filesep,  fNameRoot, filesep,  name, '.PAR']);
            movefile([subfold_pathname, filesep,  name, '.REC'], [subfold_pathname, filesep,  fNameRoot, filesep,  name, '.REC']);
            % rename the nii files in order to have the same name as the
            % Par/Rec
            %         movefile([subfold_pathname, filesep,  fNameRoot,filesep, fNameRoot, '.nii'], [subfold_pathname, filesep,  fNameRoot, filesep,  name, '.nii']);
            %     end
        end
    end
end
disp('Conversion Done!');