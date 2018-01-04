function outfiles = convert_ParUvascim2nii_BL(varargin)


filename = varargin{1};
if size(varargin) ==1
    dim = 3;
else
    dim = varargin{2};
end
options.outputformat = 1;
options.angulation = 1;%def=1
options.subaan = 0;
options.rescale = 1;
options.pathpar = [pwd filesep];
options.usefullprefix = 0;
options.usealtfolder = 0;
options.altfolder = '';

options.parname = filename;

options.prefix = '';
options.parname = [options.parname, '.PAR'];
if isempty(dim)
    options.dim = 3;
else
options.dim = dim;
end

outfiles = GFB_convert_ParUvascim2nii_BL({options.parname},options);

outfiles = outfiles{:};




