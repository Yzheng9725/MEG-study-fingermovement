% Finger extension movement
% MEG - perprocessing
% creator: WxyZ - Yu Zheng
% Date: 20250312

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

MEGDATA_CLEAN_definetrial = [];
HEADPOS_definetrial = [];

for runidx = 1:8
    opt.MEGfilename = 'xx';

    cfg              = [];
    cfg.dataset      = opt.MEGfilename;
    cfg.channel      = {'MLC', 'MLF', 'MLO', 'MLP', 'MLT',...
        'MRC', 'MRF', 'MRO', 'MRP', 'MRT',...
        'MZC', 'MZF', 'MZO', 'MZP', 'MZT'}; 
    cfg.hpfilter     = 'yes';
    cfg.hpfreq       = 0.1;
    cfg.hpfiltord     = 4;
    cfg.lpfilter     = 'yes';
    cfg.lpfreq       = 90;
    cfg.trl          = trl;  % Movement initation - detected by EMG - see func:wxyz_emgdetect
    MEGdata_definetrial  = ft_preprocessing(cfg);

    %% 
    cfg              = [];
    cfg.dataset      = opt.MEGfilename;
    cfg.channel      = {'HLC0011','HLC0012','HLC0013', ...
        'HLC0021','HLC0022','HLC0023', ...
        'HLC0031','HLC0032','HLC0033'};
    cfg.trl          = trl;
    headpos_definetrial = ft_preprocessing(cfg);

    if isempty(HEADPOS_definetrial)
        HEADPOS_definetrial = ft_appenddata([], headpos_definetrial);
    else
        HEADPOS_definetrial = ft_appenddata([], HEADPOS_definetrial, headpos_definetrial);
    end

    %% ICA
    % downsample to save time
    cfg                             = [];
    cfg.resamplefs                  = 300;
    cfg.detrend                     = 'no';
    data2ica                        = ft_resampledata(cfg, MEGdata_definetrial);

    % ICA
    cfg                             = [];
    cfg.method                      = 'runica';
    cfg.numcomponent                = 20;
    comp_ICA                        = ft_componentanalysis(cfg, data2ica);

    % show results
    cfg = [];
    cfg.layout   = 'CTF275.lay';
    comp2rej = ft_icabrowser(cfg,comp_ICA);

    % Project and clean
    cfg                             = [];
    cfg.unmixing                    = comp_ICA.unmixing;
    cfg.topolabel                   = comp_ICA.topolabel;
    comp                            = ft_componentanalysis(cfg, MEGdata_definetrial);

    cfg                             = [];
    cfg.component                   = find(comp2rej);
    MEGdata_CLEAN_definetrial       = ft_rejectcomponent(cfg, comp, MEGdata_definetrial);

    %% Visual to check ICA results
    cfg          = [];
    cfg.method   = 'trial';
    dummy        = ft_databrowser(cfg, MEGdata_definetrial);

    clear dummy

  
    %%
    if isempty(MEGDATA_CLEAN_definetrial)
        MEGDATA_CLEAN_definetrial = ft_appenddata([], MEGdata_CLEAN_definetrial);
    else
        MEGDATA_CLEAN_definetrial = ft_appenddata([], MEGDATA_CLEAN_definetrial, MEGdata_CLEAN_definetrial);
    end

     %% Headpos
    cfg         = [];
    cfg.trl     = trl;
    headpos_definetrial = ft_redefinetrial(cfg, headpos);

    if isempty(HEADPOS_definetrial)
        HEADPOS_definetrial = ft_appenddata([], headpos_definetrial);
    else
        HEADPOS_definetrial = ft_appenddata([], HEADPOS_definetrial, headpos_definetrial);
    end

end

                                                                                               
%%
MEGDATA2task = [];

for t = 1:numel(opt.task)
    cfg = [];
    cfg.trials = MEGDATA_CLEAN_definetrial.trialinfo == t;
    MEGdata2task = ft_selectdata(cfg, MEGDATA_CLEAN_definetrial);

    if isempty(MEGDATA2task)
        MEGDATA2task = ft_appenddata([], MEGdata2task);
    else
        MEGDATA2task = ft_appenddata([], MEGDATA2task, MEGdata2task);
    end
end

for t = 1:numel(opt.task)
    cfg = [];
    cfg.trials = MEGDATA2task.trialinfo == t;
    eval(['MEGDATA_' opt.task{t} '= ft_selectdata(cfg, MEGDATA2task);']);
end


%%
headpos2lf = HEADPOS_definetrial;

% cal average head movement / task
cfg = [];
cfg.dataset = opt.hdrfile;
cfg.method = 'avgoverrpt';
cfg.numclusters = 1;
headcrt4lf = ft_headmovement_wxyz(cfg, headpos2lf);

sens_avg = headcrt4lf.grad;


