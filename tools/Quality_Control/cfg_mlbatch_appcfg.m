function [cfg, def] = cfg_mlbatch_appcfg(varargin)
% 'System' - MATLABBATCH cfg_util initialisation
% This MATLABBATCH initialisation file can be used to load application
%              'System'
% into cfg_util. This can be done manually by running this file from
% MATLAB command line or automatically when cfg_util is initialised.
% The directory containing this file and the configuration file
%              'cfg_System_def'
% must be in MATLAB's path variable.
% Created at 2018-10-17 11:41:16.

if ~isdeployed
    % Get path to this file and add it to MATLAB path.
    % If the configuration file is stored in another place, the path must be adjusted here.
    p = fileparts(mfilename('fullpath'));
    addpath(genpath(p));
end
% run configuration main & def function, return output
cfg = cfg_QC;
def = [];
