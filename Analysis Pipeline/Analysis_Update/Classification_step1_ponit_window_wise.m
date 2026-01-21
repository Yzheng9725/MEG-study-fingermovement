% point & window-wise classification

%% A-window-based
Source_Amap = cell(1, numel(opt.task));

for t = 1:4
    A_task = ActIndex(t,:);
    A2trl = cat(3, A_task{:});
    Source_Amap{t}.cfg = [];
    for l = 1:length(ActIndex{t}(1,:))
        Source_Amap{t}.label{l} = ['vertex' num2str(l)];
    end
    Source_Amap{t}.dimord = 'rpt_chan_time';
    for trl = 1:length(ActIndex{t}(:,1))
        Source_Amap{t}.time{trl} = Timeseg;
        Source_Amap{t}.trial{trl} = squeeze(A2trl(trl,:,:));
    end
end

nC1 = numel(Source_Amap{1}.trial);
nC2 = numel(Source_Amap{2}.trial);
nC3 = numel(Source_Amap{3}.trial);
nC4 = numel(Source_Amap{4}.trial);

cfg = [] ;
cfg.method          = 'mvpa';
cfg.avgovertime     = 'no';
cfg.design          = [ones(nC1,1); 2*ones(nC2,1); 3*ones(nC3,1); 4*ones(nC4,1)];
cfg.features        = 'chan';
cfg.mvpa            = [];
cfg.mvpa.classifier = 'libsvm';
cfg.mvpa.metric     = 'accuracy'; 
cfg.mvpa.k          = 5;
cfg.mvpa.repeat     = 1;
cfg.mvpa.preprocess = 'zscore'; 
ClassifyResult = ft_timelockstatistics(cfg, Source_Amap{1}, Source_Amap{2},...
     Source_Amap{3}, Source_Amap{4});


%% Source-singal-based
nC1 = numel(Sourcedata_task{1}.trial);
nC2 = numel(Sourcedata_task{2}.trial);
nC3 = numel(Sourcedata_task{3}.trial);
nC4 = numel(Sourcedata_task{4}.trial);

cfg = [] ;
cfg.method          = 'mvpa';
cfg.avgovertime     = 'no';
cfg.design          = [ones(nC1,1); 2*ones(nC2,1); 3*ones(nC3,1); 4*ones(nC4,1)];
cfg.features        = 'chan';
cfg.mvpa            = [];
cfg.mvpa.classifier = 'libsvm';
cfg.mvpa.metric     = 'accuracy'; 
cfg.mvpa.k          = 5;
cfg.mvpa.repeat     = 1;
cfg.mvpa.preprocess = 'zscore'; 
ClassifyResult = ft_timelockstatistics(cfg, Sourcedata_task{1}, Sourcedata_task{2},...
    Sourcedata_task{3}, Sourcedata_task{4});


