% Finger extension movement
% MEG - Timelock & Timefrequency & ITPC & GFP % ERD/ERS
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

%%
for t = 1:numel(opt.task)
    if ~exist(['MEGDATA_' opt.task{t}],'var')
        load([opt.matdir '\MEGDATA_' opt.task{t} '.mat'])
    end
    eval(['data2process = MEGDATA_' opt.task{t} ';'])

    %     cfg = [];
    %     cfg.demean        = 'yes';
    %     cfg.baselinewindow = [-1.5 -1.0];
    %     data2process = ft_preprocessing(cfg,data2process);

    %% Timelock
    cfg = [];
    AVG = ft_timelockanalysis(cfg, data2process);

    eval(['AVG_' opt.task{t} '= AVG;']);

    %% Time Freq
    data2freq = data2process;
    data2freq.trial = cellfun(@(x) ft_preproc_padding(x,'zero',2*data2freq.fsample,2*data2freq.fsample), data2freq.trial,'UniformOutput',false);
    data2freq.time = cellfun(@(x) -5:1/data2freq.fsample:5-1/data2freq.fsample, data2freq.time,'UniformOutput',false);

    cfg                             = [];
    cfg.channel                     = 'meg';
    cfg.method                      = 'wavelet';
    cfg.output                      = 'pow';
    cfg.foi                         = 0.5:0.5:90;
    cfg.toi                         = -2.5:0.1:2.5;
    cfg.width                       = 4;  
    TFWAVE_AVG = ft_freqanalysis(cfg, data2freq);

    eval(['TFWAVE_AVG_' opt.task{t} '= TFWAVE_AVG;']);

    %%  ITPC
    cfg = [];
    cfg.method = 'wavelet';
    cfg.output = 'fourier';
    cfg.toi = -2.0:0.02:2;
    cfg.foi = 0.5:0.5:90;
    cfg.width  = 4;
    freq = ft_freqanalysis(cfg, data2freq);

    itc = [];
    itc.label = freq.label;
    itc.freq = freq.freq;
    itc.time = freq.time;
    itc.dimord = 'chan_freq_time';

    F = freq.fourierspctrm;
    N = size(F,1);

    itc.itpc = F./abs(F);
    itc.itpc = sum(itc.itpc,1);
    itc.itpc = abs(itc.itpc)/N;
    itc.itpc = squeeze(itc.itpc);

    % TFR plot
    % itc2plot = itc;
    % itc2plot.cumtapcnt = freq.cumtapcnt;
    % itc2plot.powspctrm = itc.itpc;

    eval(['ITPC_' opt.task{t} '= itc;']);

    %% GFP
    % per task
    eval(['AVG2gfp = AVG_' opt.task{t} ';']);

    cfg = [];
    cfg.method = 'amplitude';
    AVG_gfp = ft_globalmeanfield(cfg, AVG2gfp);

    eval(['GFP_' opt.task{t} '= AVG_gfp;']);

    % all task avg
    %     AVG2gfp = AVG_Thumb;
    %     AVG2gfp.avg = (AVG_Thumb.avg+AVG_Index.avg+AVG_Middle.avg+AVG_Little.avg)/4;
    %
    %     cfg = [];
    %     cfg.method = 'amplitude';
    %     AVG_gmfp = ft_globalmeanfield(cfg, AVG2gfp);

    %% ERD/ERS
    FreqName = {'Alpha','Beta','High-Gamma'};
    FreqBand = {[8 13],[13 30],[60 90]};

    for r = 1:numel(FreqName)
        opt.freqname = FreqName{r};  
        opt.freqband = FreqBand{r};  

        cfg = [];
        cfg.resamplefs = 300;
        data2erds = ft_resampledata(cfg, data2process);

        ERDS = wxyz_erdsanalysis(data2erds, opt.freqband, baseline);
        eval(['ERDS_' opt.freqname '_' opt.task{t} '= ERDS;']);
    end

end


%%  check
baseline = [-1.5 -1.0];

% Timelock
cfg                 = [];
cfg.baseline        = baseline;
cfg.layout          = 'CTF275_helmet.mat';
cfg.colorbar        = 'no';
cfg.masknans        = 'no';
cfg.showlabels      = 'yes';
cfg.xlim            = [-2.0 2.0];
cfg.ylim            = [-2.5e-13 2.5e-13];
cfg.showoutline     = 'yes'; 
cfg.showscale       = 'no'; 
cfg.linewidth       = 1;
ft_multiplotER(cfg, AVG);

% TimeFreq
cfg                 = [];
cfg.baseline        = baseline;
cfg.baselinetype    = 'absolute';
cfg.xlim            = [-2.0 2.0];
cfg.ylim            = [0.5 90];
cfg.zlim            = [-1e-24 1e-24];
cfg.layout          = 'CTF275_helmet.mat';
cfg.colorbar        = 'no';
cfg.masknans        = 'no';
cfg.showlabels      = 'no';
cfg.showoutline     = 'yes'; 
cfg.showscale       = 'no'; 
cfg.figure          = gcf;

for t = 1:length(opt.task)
    figure
    eval(['ft_multiplotTFR(cfg, TFWAVE_AVG_' opt.task{t} ');']);
end






