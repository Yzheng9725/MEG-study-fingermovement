% classification - different feature extraction
% PCA & A_epoch & A_window
% classification using SVM

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

%% source reconstruction - ROI
Sourcedata = [];
Sourcedata.time(1:length(data2ana_all.trial)) = data2ana_all.time{1};
for ch = 1:numel(ROIidx)
    Sourcedata.label(ch,1) = {['Vertex ' num2str(ROIidx(ch))]};
end

for trl = 1:length(data2ana_all.trial)
    Sourcedata.trial{trl} = commonfilter_ROI * data2ana_all.trial{trl}; % use commonfilter
end

Sourcedata.trialinfo = data2ana_all.trialinfo;
Sourcedata.fsample = data2ana_all.fsample;

%% PCA
cfg = [];
cfg.method       = 'pca';
cfg.numcomponent = 5;  % >70% var
Source_pca = ft_componentanalysis(cfg, Sourcedata);

% >70% var
% dat = cell2mat(Sourcedata.trial);
% C1 = (dat*dat')./(size(dat,2)-1);
% [~,D] = eig(C1);                                    % eigenvalue decomposition (EVD)
% d = cat(2,(1:1:numel(Sourcedata.label))',diag(D));  % sort eigenvectors in descending order of eigenvalues
% d = sortrows(d, -2);
% explained_var = d(:,2) / sum(d(:,2));
% cum_explained = cumsum(explained_var);
% compnum2save = find(cum_explained >= 0.7, 1);

% select data for each task
for t = 1:numel(opt.task)
    trl2sel = find(Source_pca.trialinfo == t);
    cfg = [];
    % cfg.channel = Label2sel;
    cfg.trials = trl2sel(1:100);
    Sourcedata_task = ft_selectdata(cfg, Source_ROI_pca_70var_alltask);
end

%%  A-epoch
Actseg = [-1 2.0];
% cross validation
clfopt.trialnum = 100;
clfopt.kfold = 5;
clfopt.Label_all = [ones(1,clfopt.trialnum)*1 ones(1,clfopt.trialnum)*2 ones(1,clfopt.trialnum)*3 ones(1,clfopt.trialnum)*4]';
clfopt.ClassIdx_all = crossvalind('Kfold',clfopt.Label_all,clfopt.kfold);
clfopt.TaskLabel = {'Thumb', 'Index', 'Middle','Little'};

for k = 1:5
    trl2train = find(clfopt.ClassIdx_all ~= k);
    % clearvars AI_alltask_kt_avg
    ActIndex = cell(1, numel(opt.task));
    maxAI = nan(numel(opt.task), 1); maxidx = nan(numel(opt.task), 1);
    minAI = nan(numel(opt.task), 1); minidx = nan(numel(opt.task), 1);

    for t = 1:numel(opt.task)
        trl2train_task = trl2train((t-1)*80+1:80*t) - (t-1)*100;
        ActIndex{t} = mean(AI_alltask_kt.AI_alltask_kt{t}(trl2train_task,:));

        % find max&min-A location
        [maxAI(t,:), maxidx(t,:)] = max(ActIndex{t},[],2);
        [minAI(t,:), minidx(t,:)] = min(ActIndex{t},[],2);

        posidx_max(t,:) = ROIidx(maxidx(t,:));
        posidx_min(t,:) = ROIidx(minidx(t,:));
    end

    %% Source reconstruction
    posidx = [posidx_max; posidx_min]; 
    commonfilter = cell2mat(commonfilter_wb(posidx));
    Sourcedata_lcmv = cell(1, numel(opt.task));

    for t = 1:numel(opt.task)
        bmfilter = commonfilter;
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

    % append data
    Sourcedata_classify = [];
    for t = 1:numel(opt.task)
        if t == 1
            Sourcedata_classify = Sourcedata_lcmv{t};
        else
            Sourcedata_classify = ft_appenddata([], Sourcedata_classify, Sourcedata_lcmv{t});
        end
    end
end

%% A-window
% select max/min-A per window
Timewin = Actseg(1):0.05:Actseg(2);
timewin2clf = knnsearch(Timewin', Actseg(1)):knnsearch(Timewin', Actseg(2));
posidx4clf_max = posidx_max(:,timewin2clf);
posidx4clf_min = posidx_min(:,timewin2clf);

for  t = 1:4
    pos2plot = sourcemodel.pos(posidx4clf_max(t,:),:);
    pos2plot_min = sourcemodel_inflated.pos(posidx4clf_min(t,:),:);
end

% Source reconstruction
posidx = reshape([posidx4clf_max posidx4clf_min]', 1, []);  
commonfilter = cell2mat(commonfilter_wb(posidx));

for t = 1:numel(opt.task)
    bmfilter = commonfilter;
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

% append data
Sourcedata = [];
Sourcedata4clf = cell(1, numel(opt.task));
for t = 1:numel(opt.task)
    cfg = [];
    cfg.latency = Actseg;
    data2append = ft_selectdata(cfg, Sourcedata_lcmv{t});

    %%
    cfg = [];
    cfg.channel = data2append.label(1:2*numel(opt.task));
    Sourcedata4clf{t} = ft_selectdata(cfg, data2append);

    for trl = 1:numel(data2append.trial)
        for tw = 1:numel(timewin2clf)
            timwin2append = knnsearch(data2append.time{1}', Timewin(timewin2clf(tw)));
            for tch = 1:2*numel(opt.task)
                if tw == numel(timewin2clf)
                    data4timwin{trl}(tch,timwin2append:timwin2append+15) =...
                        data2append.trial{trl}((tch-1)*numel(timewin2clf)+tw, timwin2append:timwin2append+15);
                else
                    data4timwin{trl}(tch,timwin2append:timwin2append+14) =...
                        data2append.trial{trl}((tch-1)*numel(timewin2clf)+tw, timwin2append:timwin2append+14);
                end
            end
        end
    end

    Sourcedata4clf{t}.trial = data4timwin;
   
    if t == 1
        Sourcedata_classify = Sourcedata4clf{t};
    else
        Sourcedata_classify = ft_appenddata([], Sourcedata_classify, Sourcedata4clf{t});
    end
end


