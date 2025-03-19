% Finger extension movement
% MEG - ITPC & GFP & ERD/ERS - Group-level statistic
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
opt.matdir = 'xx\';

%% ITPC stat
% need ITPC mat for all tasks and all individual
allsub = [];
allbase = [];

for ss = 1:16
    % need ITPC mat for all tasks and all individual - task_itpc

    base_itpc = task_itpc;
    for c = 1:numel(base_itpc.label)
        for f = 1:numel(base_itpc.freq)
            base_itpc.itpc(c,f,:) = ones(size(base_itpc.time))*mean(base_itpc.itpc(c,f,26:51));  % baseline period
        end
    end

    allsub{ss} = task_itpc;
    allbase{ss} = base_itpc;
end

% statistic
cfg = [];
cfg.latency = [-1.5 2];
cfg.method = 'montecarlo';
cfg.parameter = 'itpc';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 5;
cfg_neigh.method = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg_neigh, itc);
cfg.tail = 1;  % only positive
cfg.clustertail = 1;
cfg.alpha = 0.05;
cfg.numrandomization = 2000;

Nsubj = 16;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1, Nsubj) ones(1, Nsubj)*2];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

ITPC_stat = ft_freqstatistics(cfg, allsub{:}, allbase{:});



%% GFP stat
% need GFP mat for all tasks and all individual
baseline_mean = mean(GFP_alltask(:,1201:1800),2);  % baseline: -1.5 ~ -1.0 s

allsub = [];
allbase = [];
for ss = 1:16  % 16 participants
    gfp_task = AVG2gfp;
    gfp_task.avg(1,:) = GFP_alltask(ss,:);

    cfg = [];
    cfg.channel = AVG2gfp.label(1);
    gfp_task = ft_selectdata(cfg, gfp_task);

    gmfp_base = gfp_task;
    gmfp_base.avg(1,:) = ones(size(gmfp_base.avg(1,:)))*baseline_mean(ss);

    allsub{ss} = gfp_task;
    allbase{ss} = gmfp_base;
end

cfg = [];
cfg.latency = [-1.5 2.0];  % baseline + epoch periods

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'fdr';
cfg.clusterstatistic = 'maxsum';
cfg.neighbours = [];
cfg.tail = 1;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 5000;

Nsubj = 16;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1, Nsubj) ones(1, Nsubj)*2];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

GFP_stat = ft_timelockstatistics(cfg, allsub{:}, allbase{:});

% the significant GFP time period were used as low-frequency task activation period


%%  ERD/ERS stat
% need ERD/ERS mat for all tasks and all individual

% baseline period
baseline_mean = mean(ERDS_alltask(:,:,knnsearch(ERDS_alltaskavg.time',-1.5):knnsearch(ERDS_alltaskavg.time',-1.0)),3);

allsub = [];
allbase = [];
for i = 1:16
    erds_task = ERDS_alltaskavg;
    erds_task.avg(:,:) = squeeze(ERDS_alltask(i,:,:));

    erds_base = erds_task;
    erds_base.avg(:,:) = ones(size(erds_base.avg(:,:)))*baseline_mean(i);

    allsub{i} = erds_task;
    allbase{i} = erds_base;
end


cfg = [];
cfg.latency = [-1.5 2];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg_neigh.method = 'distance';
cfg.neighbours = ft_prepare_neighbours(cfg_neigh, ERDS_alltaskavg);
cfg.tail = 0;  % positive + negtive
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 2000;

Nsubj = 16;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1, Nsubj) ones(1, Nsubj)*2];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

ERDS_stat = ft_timelockstatistics(cfg, allsub{:}, allbase{:});

% the significant ERD/ERS time periods were used in sourcelevel peak time selection (Ana3_2)





