% Finger extension movement
% MEG - Sourcelevel - Low-frequency - event-related SAM
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

load([opt.mripath 'mri.mat'])
load('xx\headmodel.mat')
load('xx\sourcemodel.mat')
load('xx\leadfield_avg.mat')
load('xx\sens_avg.mat')

baseline = [-1.5 -1.0];
epoch = [-0.5 2.0];

%% create vertex ROI
% atlas
load([opt.mripath 'brainatlas.mat'])
atlas_wb_L = brainatlasL;
atlas_wb_R = brainatlasR;

atlasLabel = [21 23];        % label:Postcentral/Precentral   aparc
Atlasannot_Lshpere = [atlas_wb_L.parcellation];      % brain area label in atlas, only left

% select ROI pos index
ROIidx = [];
for i = 1:numel(atlasLabel)
    idx = find(Atlasannot_Lshpere == atlasLabel(i));
    ROIidx = [ROIidx;idx];
    clearvars idx
end

%% change leadfield.inside 2 keep ROI
leadfield.inside = false(size(leadfield.inside));
leadfield.inside(ROIidx) = true;
sourcemodel.inside = leadfield.inside;

for t = 1:length(opt.task)
    %% load MEGDATA_task
    eval(['data2source = MEGDATA_' opt.task{t} ';'])
    data2source.grad = sens; % replace grad -> grad_alltrials_avg

    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    % cfg.demean        = 'yes';
    % cfg.baselinewindow = baseline;
    data2sam = ft_preprocessing(cfg, data2source);

    %% SAM - 20 ms time window
    % use baseline to cal noise_avg/noise_avg.cov
    cfg = [];
    cfg.toilim = baseline;
    dataNoise_sam = ft_redefinetrial(cfg, data2sam);

    cfg = [];
    cfg.covariance = 'yes';
    avgNosie_sam = ft_timelockanalysis(cfg, dataNoise_sam);

    cfg = [];
    cfg.toilim = epoch;
    data_sam = ft_redefinetrial(cfg, data2sam);

    cfg = [];
    cfg.covariance = 'yes';
    avg_sam = ft_timelockanalysis(cfg, data_sam);

    Timeseg = [];
    Timestep = 0.02;
    Timeseg = -0.5:Timestep:1-Timestep;  % Low-frequency task activation

    for c = 1:length(Timeseg)
        Time2cal = knnsearch(avg_sam.time', Timeseg(c));

        % cal event-related SAM, using computed leadfield
        Source_sam_avg = ft_inverse_sam(leadfield, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
            'reducerank',2, 'toi', [Time2cal Time2cal+Timestep*data2sam.fsample], 'noisecov',avgNosie_sam.cov);

        [MaxPseudoZ, MaxPosidx] = max(Source_sam_avg.pseudoZ);
        if ~ismember(MaxPosidx,ROIidx)
            [MaxPseudoZ, ROIPosidx] = max(Source_sam_avg.pseudoZ(ROIidx));
            MaxPosidx = ROIidx(ROIPosidx);
        end

        % select max pseudoZ position -> source localization; save pos & ori & mom & posidx(whole sourcemodel.pos)
        Source_sam_max(c).time = Timeseg(c);
        Source_sam_max(c).posidx = MaxPosidx;
        Source_sam_max(c).maxpseudoZ = MaxPseudoZ;

    end

end


%% find pZ peak window
for t = 1 : length(opt.task)
    v2find = nan(1,numel(Source_sam_max));
    for i = 1:numel(Source_sam_max)
        v2find(i) = Source_sam_max(i).maxpseudoZ;
    end

    timeidx = {[],[],[],[]};  % 4 Peak timeidx of group-avg GFP results in [-0.5 1.0]

    for c = 1:4
        % find peak
        [o, p] =findpeaks(v2find(timeidx{c}));
        if isempty(o)
            [peakv(c), pidx(c)] = max(v2find(timeidx{c}));
            peakidx(c) = timeidx{c}(1) + pidx(c) -1;
        else
            [peakv(c), ~] = max(findpeaks(v2find(timeidx{c})));
            peakidx(c) = find(v2find == peakv(c));
        end

        Time4peak(t,c) = Source_sam_max(peakidx(c)).time;
    end
end


%% refine the peak periods 
for t = 1:numel(opt.task)
    % load MEGDATA_task
    eval(['data2source = MEGDATA_' opt.task{t} ';'])

    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    % cfg.demean        = 'yes';
    % cfg.baselinewindow = baseline;
    data2sam = ft_preprocessing(cfg, data2source);

    % use baseline to cal noise_avg/noise_avg.cov
    cfg = [];
    cfg.toilim = baseline;
    dataNoise_sam = ft_redefinetrial(cfg, data2sam);

    cfg = [];
    cfg.covariance = 'yes';
    avgNosie_sam = ft_timelockanalysis(cfg, dataNoise_sam);

    cfg = [];
    cfg.toilim = epoch;
    data_sam = ft_redefinetrial(cfg, data2sam);

    cfg = [];
    cfg.covariance = 'yes';
    avg_sam = ft_timelockanalysis(cfg, data_sam);

    Time2cal = [];
    for i = 1:numel(Time4peak(t,:))
        Time2cal(i,:) = Time4peak(t,i)-0.015:0.005:Time4peak(t,i)+0.015;
    end
    Time2cal = reshape(Time2cal,1,[]);

    %%
    Source_sam_max = [];
    for c = 1:numel(Time2cal)
        Time2cal = knnsearch(avg_sam.time', Time2cal(c));

        % cal event-related SAM, using computed leadfield
        Source_sam_avg = ft_inverse_sam(leadfield_ROI, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
            'reducerank',2, 'toi', [Time2cal Time2cal+Timestep*data2sam.fsample], 'noisecov',avgNosie_sam.cov);

        [MaxPseudoZ, MaxPosidx] = max(Source_sam_avg.pseudoZ);

        %%
        Source_sam_max(c).time = Time2cal(c);
        Source_sam_max(c).posidx = MaxPosidx;
        Source_sam_max(c).maxpseudoZ = Source_sam_avg.pseudoZ(MaxPosidx);

    end

    %% find peak periods
    timenum = numel(Time2cal)/numel(Time4peak);
    for c = 1:numel(Time4peak)
        v2find = [];
        for i = 1:timenum
            v2find(i) = Source_sam_max((c-1)*timenum + i).maxpseudoZ;
        end
        [~, idx2save] = max(v2find);
        Time2max(c) =  Source_sam_max((c-1)*timenum + idx2save).time;
    end

end





