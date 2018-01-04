% customize the macro (for this specific file"


[macro_pathstr,macro_name,~]  = fileparts(which('parametric_map_histology_vascular_staining_ImageJ'));

tiff_fileliste = dir('*.jpg');

for i=1:size(tiff_fileliste,1)
    [file_path, file_name, file_ext] =fileparts([pwd, filesep, tiff_fileliste(i).name]);
    file_path = strrep(file_path, filesep, [filesep filesep]);
    macro_fid = fopen(fullfile(macro_pathstr, [macro_name '_macro.txt']));
    
    if (macro_fid>0)
        macro=fread(macro_fid, '*char');
        fclose(macro_fid);
        
        macro=macro(:)';
        macro=strrep(macro, 'file_path', [file_path filesep filesep]);
        macro=strrep(macro, 'file_name', file_name);
        macro=strrep(macro, 'file_ext', file_ext);
        custum_macroformatlab = fopen(fullfile(file_path, [file_name '_macro_tmp.txt']), 'wt');
        fwrite(custum_macroformatlab, macro(:));
        fclose(custum_macroformatlab);
    end
end


%open the imageJ GUI
% imagej_GUI = Miji();
% IJ=ij.IJ();
% execture the Imagej's macro to quantify the data
% IJ.runMacroFile(fullfile(macro_pathstr, [macro_name '_macro_tmp.txt']));
% 
% % close the imageJ GUI
% imagej_GUI.exit
% 
% % delete the custum macro
% delete(fullfile(macro_pathstr, [macro_name '_macro_tmp.txt']))
%clear variables