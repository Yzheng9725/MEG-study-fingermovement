% classification using features from multi-window 

%% data preprocessed
opt.freqband = [0.5 4];  % [4 8],[8 13],[13 30],[30 60],[60 90]
opt.freqname = 'Delta';  % 'Theta','Alpha','Beta','Gamma','High-gamma'
opt.task = {'Thumb','Index','Middle','Little'};

cfg                 = [];
cfg.lpfilter        = 'yes';
cfg.lpfreq          = opt.freqband(2);
cfg.hpfilter        = 'yes';
cfg.hpfreq          = opt.freqband(1);
cfg.demean          = 'yes';
cfg.baselinewindow  = [-2 -1.5];
data2ana = cell(1, numel(opt.task));
for t = 1:numel(opt.task)
    data2ana{t}     = ft_preprocessing(cfg, MEGDATA_task{t});
end
% append data
data2ana_all = ft_appenddata([],data2ana{1},data2ana{2},data2ana{3},data2ana{4});

%% feature extraction & classification
clfopt = [];
clfopt.trialnum = 100;
clfopt.kfold = 5;
clfopt.Label_all = [ones(1,clfopt.trialnum)*1 ones(1,clfopt.trialnum)*2 ones(1,clfopt.trialnum)*3 ones(1,clfopt.trialnum)*4]';
clfopt.ClassIdx_all = crossvalind('Kfold',clfopt.Label_all,clfopt.kfold);
clfopt.TaskLabel = {'Thumb', 'Index', 'Middle','Little'};

for k = 1:5
    posROIidx = [];
    trl2train = find(clfopt.ClassIdx_all ~= k);

    maxAI = nan(numel(opt.task), length(Timeseg)); maxidx = nan(numel(opt.task), length(Timeseg));
    minAI = nan(numel(opt.task), length(Timeseg)); minidx = nan(numel(opt.task), length(Timeseg));
    posidx_max = nan(numel(opt.task), length(Timeseg));
    posidx_min = nan(numel(opt.task), length(Timeseg));
    ActIndex_classify = cell(1, numel(opt.task));

    for t = 1:4
        trl2train_task = trl2train((t-1)*80+1:80*t) - (t-1)*100;
        for c = 1:length(Window_decodable)  % results in A-window-based window-wise classification
            ActIndex_classify{t}(c,:) = mean(ActIndex{t,c}(trl2train_task,:));
        end

        % find max/min-A location
        [maxAI(t,:), maxidx(t,:)] = max(ActIndex_classify{t},[],2);
        [minAI(t,:), minidx(t,:)] = min(ActIndex_classify{t},[],2);
        posidx_max(t,:) = ROIidx(maxidx(t,:));
        posidx_min(t,:) = ROIidx(minidx(t,:));
    end

    %%  Source reconstruction
    Timewin = Actseg(1):0.05:Actseg(2);
    timewin2clf = knnsearch(Timewin', Actseg(1)):knnsearch(Timewin', Actseg(2));
    posidx4clf_max = posidx_max(:,timewin2clf);
    posidx4clf_min = posidx_min(:,timewin2clf);
    posidx = reshape([posidx4clf_max posidx4clf_min]', 1, []);  % virtual channel set locations
    Sourcedata_lcmv = cell(1, numel(opt.task));

    for t = 1:numel(opt.task)
        bmfilter = commonfilter_ROI;
        Sourcedata_lcmv{t} = [];
        Sourcedata_lcmv{t}.time(1:length(data2ana{t}.trial)) = data2ana{t}.time(1);
        for ch = 1:numel(posidx)
            Sourcedata_lcmv{t}.label(ch,1) = {['Vertex ' num2str(ch)]};
        end

        for trl = 1:length(MEGDATA_task{t}.trial)
            Sourcedata_lcmv{t}.trial{trl} = bmfilter * data2ana{t}.trial{trl};
        end

        Sourcedata_lcmv{t}.trialinfo = t * ones(size(data2ana{t}.trialinfo));
    end

    % select data & append data
    for t = 1:numel(opt.task)
        if t == 1
            Sourcedata_classify{1} = Sourcedata_lcmv{t};
        else
            Sourcedata_classify{1} = ft_appenddata([], Sourcedata_lcmv{1}, Sourcedata_lcmv{t});
        end
    end

    %% train/test set
    cfg = [];
    cfg.latency = Time_decodable;  % results in source reconstruction signal-based point-wise classification
    cfg.trials = find(clfopt.ClassIdx_all ~= k);
    train2feature = ft_selectdata(cfg, Sourcedata_classify{clfnum});

    cfg.trials = find(clfopt.ClassIdx_all == k);
    test2feature = ft_selectdata(cfg, Sourcedata_classify{clfnum});

    % resample data
    cfg = [];
    cfg.resamplefs = 30;
    train2feature = ft_resampledata(cfg, train2feature);
    test2feature = ft_resampledata(cfg, test2feature);

    DATA2train = train2feature.trial;
    DATA2test = test2feature.trial;

    trainset = cellfun(@(x) reshape(x', 1, []), DATA2train,'UniformOutput',false);
    testset = cellfun(@(x) reshape(x', 1, []), DATA2test,'UniformOutput',false);
    trainset = cell2mat(trainset');
    testset = cell2mat(testset');
end