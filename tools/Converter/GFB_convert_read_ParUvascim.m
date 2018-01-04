function parrec = GFB_convert_read_ParUvascim(varargin)
% function parrec = GFB_convert_read_parrec(parfilename, uvascim_filename)
%
% Using the PAR header file defined by parfilename and REC data format file
% in uvascim_filename, load the header data and image data. 
%
% (C)opyright 2005, Bennett Landman, bennett@bme.jhu.edu
% Revision History:
% Created: 2/11/2005
% 3/26/05 Forced correction for floating point images
% 10/23/2006 Simplified error handling 
% 11/5/2006 Added de-interleaving for PAR/REC so that save/write is
% compatable and simple
% 26/8/2010 remove option to load only selected scans (IT, JW)

if(nargin==2)
    parfilename=varargin{1};
    uvascim_filename=varargin{2};
end
if(nargin==1)
    [pth,fl] = fileparts(varargin{1});
    varargin{1} = [pth filesep fl];
    parfilename=[varargin{1} '.par'];
    if(exist(parfilename,'file') || exist([parfilename 'v2'],'file'))        
        uvascim_filename=[varargin{1} '.mat'];
        parrec = GFB_convert_read_ParUvascim(parfilename,uvascim_filename);
        parrec.caps = 0;        
        return;
    else
        parfilename=[varargin{1} '.PAR'];
        uvascim_filename=[varargin{1} '.REC'];
        parrec = GFB_convert_read_ParUvascim(parfilename,uvascim_filename);
        parrec.caps = 1;        
        return;
    end
end
parrec.hdr = GFB_convert_read_par(parfilename);

[parrec.scans parrec.hdr] = GFB_convert_read_Uvascim2(uvascim_filename, parrec.hdr,1);