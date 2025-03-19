% Finger extension movement
% MRI - perprocessing - need freesurfer
% creator: WxyZ - Yu Zheng
% Date: 20250314

%%
clc; clear; close all;
ft_defaults;
ft_version;
[ftv, ftpath] = ft_version;

%%
opt = [];
opt.filepath = 'xx\';
opt.mripath = 'xx\';
opt.matdir = 'xx\';

%% read mri - dicom
mri = ft_read_mri([opt.mripath 'xxx']);

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
mri_realigned = ft_volumerealign(cfg,mri);
% Assign the nasion (pressing "n"), left ("l") and right ("r") with the crosshairs on the ear markers. Then finish with "q".

cfg  = [];
cfg.resolution = 1;
cfg.dim        = [256 256 256];
mri_reslice = ft_volumereslice(cfg,mri_realigned);

ft_determine_coordsys(mri_reslice,'interactive','no');

%% prepare mgz for Freesurfer
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
mri2freesurfer = ft_volumerealign(cfg, mri_reslice);
ft_determine_coordsys(mri2freesurfer,'interactive','no');

mri2freesurfer = ft_convert_coordsys(mri_reslice,'acpc');
transform_vox2acpc = mri2freesurfer.transform;
transform_vox2ctf = mri_reslice.transform;

cfg          = [];
cfg.filename = fullfile(opt.mripath,sprintf('%s_2fs.mgz',opt.sub));
cfg.filetype = 'mgz';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri2freesurfer);


%% prepare headmodel
cfg = [];
cfg.output = 'brain';
segmentedmri = ft_volumesegment(cfg, mri_reslice);

cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, segmentedmri);

%% Create leadfield 
cfg = [];
cfg.channel    = label; % ensure that rejected sensors are not present
cfg.grad       = sens_avg;
cfg.headmodel  = headmodel;
cfg.sourcemodel = sourcemodel;  % workbench results
cfg.reducerank = 2; % default for MEG is 2, for EEG is 3
leadfield = ft_prepare_leadfield(cfg);








