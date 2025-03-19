% Finger extension movement
% MEG - Digit map
% creator: WxyZ - Yu Zheng
% Date: 20250313

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
% need clfopt

baseline = [-1.5 -1.0];
epoch = [-0.5 2.0];

%% create vertex ROI
%%%%%  NOTE: need index of vertice in Active-ROI! %%%%%

%% Source space Projection - erSAM
for k = 1:5
    %% cal baseline 2 cal NOISE cov
    for t = 1: numel(opt.task)
        %% Select trials 2 calculate common COV
        ClassIdx4task = clfopt.ClassIdx_all((t-1)*clfopt.trialnum+1 : t*clfopt.trialnum);

        %% load raw trial data
        eval(['data2ana = MEGDATA_' opt.task{t} ';'])

        cfg = [];
        cfg.trials = find(ClassIdx4task ~= k);
        data2base = ft_selectdata(cfg, data2ana);

        if t == 1
            DATA2BASE = ft_appenddata([],data2base);
        else
            DATA2BASE = ft_appenddata([],DATA2BASE,data2base);
        end
    end

    %% cal cov
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    % cfg.demean        = 'yes';
    % cfg.baselinewindow = baseline;
    DATA2BASE = ft_preprocessing(cfg, DATA2BASE);

    % use baseline to cal noise_avg/noise_avg.cov
    cfg = [];
    cfg.toilim = baseline;
    dataNoise_sam = ft_redefinetrial(cfg, DATA2BASE);

    cfg = [];
    cfg.covariance = 'yes';
    avgNoise_sam = ft_timelockanalysis(cfg, dataNoise_sam);

    %% cal SAM
    Timeseg = [];
    Timestep = 0.02;
    Timeseg = [];  % time window of group-level Active-ROI exceeding chance level

    for t = 1: numel(opt.task)
        leadfield_ROI = leadfield;
        sourcemodel_ROI = sourcemodel;
        leadfield_ROI.inside = false(size(leadfield.inside));
        leadfield_ROI.inside(ROIidx) = true;
        sourcemodel_ROI.inside = leadfield_ROI.inside;

        %% Select trials 2 calculate common COV
        ClassIdx4task = clfopt.ClassIdx_all((t-1)*clfopt.trialnum+1 : t*clfopt.trialnum);

        %% load raw trial data
        eval(['data2ana = MEGDATA_' opt.task{t} ';'])
     
        cfg = [];
        cfg.trials = find(ClassIdx4task ~= k);
        data2task = ft_selectdata(cfg, data2ana);

        %% cal cov
        cfg = [];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 8;
        % cfg.demean        = 'yes';
        % cfg.baselinewindow = baseline;
        data2sam = ft_preprocessing(cfg, data2task);

        cfg = [];
        cfg.toilim = epoch;
        dataAvg_sam = ft_redefinetrial(cfg, data2sam);

        cfg = [];
        cfg.covariance = 'yes';
        avg_sam = ft_timelockanalysis(cfg, dataAvg_sam);

        %% cal Weights  - TOI: slip windows -500 - 1000ms  by 20ms step
        filter_sam = nan(numel(Timeseg),sum(sourcemodel_ROI.inside),numel(label));
        ori_sam{t} = nan(numel(Timeseg),sum(sourcemodel_ROI.inside),3);
        ori_correct{t} = nan(numel(Timeseg),sum(sourcemodel_ROI.inside),3);
        filter4recon{t} = nan(numel(Timeseg),sum(sourcemodel_ROI.inside),numel(label));

        for c = 1:numel(Timeseg)
            Time2cal = knnsearch(avg_sam.time', Timeseg(c));

            % cal event-related SAM, using computed leadfield
            Source_sam_avg = ft_inverse_sam(leadfield_ROI, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
                'reducerank',2, 'toi', [Time2cal Time2cal+Timestep*data2sam.fsample],  'noisecov', avgNoise_sam.cov);   %

            % Whole brain + ROI inside: space filter & pseudoZ & ori  2  SAVE
            filter_sam(c, :, :) = cell2mat(Source_sam_avg.filter');
            ori_sam{t}(c, :, :) = cell2mat(Source_sam_avg.ori)';

            %% adjust ori
            vertices = sourcemodel.pos;
            tri = sourcemodel.tri;
            ori4cal = squeeze(ori_sam{t}(c, :, :));

            for i = 1:numel(ROIidx)
                [f,~] = find(tri == ROIidx(i));
                avgNormal = [0 0 0];
                for j = 1:length(f)
                    face = tri(f(j),:);
                    v1 = vertices(face(2),:) - vertices(face(1),:);
                    v2 = vertices(face(3),:) - vertices(face(1),:);
                    faceNormal = cross(v1, v2);
                    avgNormal = avgNormal + faceNormal / norm(faceNormal);
                end
                avgNormal = avgNormal / length(f);
                avgNormal = avgNormal / norm(avgNormal);
                avgNormal = -1 * avgNormal;

                dotv = dot(avgNormal,ori4cal(i, :));
                if dotv >= 0
                    ori_correct{t}(c, i, :) = ori_sam{t}(c, i, :);
                    filter4recon{t}(c, i, :) = filter_sam(c, i, :);
                else
                    ori_correct{t}(c, i, :) = -1*ori_sam{t}(c, i, :);
                    filter4recon{t}(c, i, :) = -1*filter_sam(c, i, :);
                end
            end

        end
    end



    %% Data Recon
    for t = 1:numel(opt.task)
        %% load raw trial data
        eval(['data2ana = MEGDATA_' opt.task{t} ';'])

        if t == 1
            DATA2ANA = ft_appenddata([], data2ana);
        else
            DATA2ANA = ft_appenddata([], DATA2ANA, data2ana);
        end

    end

    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    % cfg.demean        = 'yes';
    % cfg.baselinewindow = baseline;
    data2recon = ft_preprocessing(cfg, DATA2ANA);

    %%
    Acc = [];  % zeros(numel(Timeseg),numel(ROIidx))

    for c = 1:numel(Timeseg)
        %% estimate ROI-prob MAP per window
        ReconData_train = [];
        ReconData_test = [];

        for t = 1:numel(opt.task)
            disp(['-------- ' opt.task{t} ' --------']);
            ClassIdx4task = clfopt.ClassIdx_all((t-1)*clfopt.trialnum+1 : t*clfopt.trialnum);

            %% 划分train/test set
            cfg = [];
            cfg.trials = find(data2recon.trialinfo == t);
            data2recon_task = ft_selectdata(cfg, data2recon);

            cfg = [];
            cfg.toilim = [Timeseg(c) Timeseg(c)+Timestep];
            cfg.trials = find(ClassIdx4task ~= k);
            data2recon_train = ft_selectdata(cfg, data2recon_task);

            cfg.trials = find(ClassIdx4task == k);
            data2recon_test = ft_selectdata(cfg, data2recon_task);

            ReconData_train{t} = nan(numel(data2recon_train.trial),numel(ROIidx),numel(data2recon_train.time{1}));
            ReconData_test{t} = nan(numel(data2recon_test.trial),numel(ROIidx),numel(data2recon_test.time{1}));

            Filter4cal = squeeze(filter4recon{t}(c, :, :));
            Filter4cal(isnan(Filter4cal(:,1)),:) = [];

            for i = 1:numel(data2recon_train.trial)
                ReconData_train{t}(i,:,:) = Filter4cal * data2recon_train.trial{i};
            end

            for i = 1:numel(data2recon_test.trial)
                ReconData_test{t}(i,:,:) = Filter4cal * data2recon_test.trial{i};
            end
        end


        %% Classification
        trainlabel = clfopt.Label_all(clfopt.ClassIdx_all ~= k);
        testlabel = clfopt.Label_all(clfopt.ClassIdx_all == k);

        Vol2clf = Map_task_wholebrain(c,ROIidx);
        Idx2clf = find(Vol2clf ~= 0);


        for n = 1:numel(Idx2clf)
            trainset = [];
            testset =[];

            for i = 1:numel(ReconData_train)
                trainset = [trainset; squeeze(ReconData_train{i}(:,Idx2clf(n),:))];
                testset = [testset; squeeze(ReconData_test{i}(:,Idx2clf(n),:))];
            end

            %% normalization
            [~,PStrian] = mapminmax(trainset');
            PStrian.ymin = -1;
            PStrian.ymax = 1;
            [trainset1,PStrian] = mapminmax(trainset,PStrian);

            [~,PStest] = mapminmax(testset');
            PStest.ymin = -1;
            PStest.ymax = 1;
            [testset1,PStest] = mapminmax(testset,PStest);


            %%
            model = svmtrain(trainlabel, trainset1,'-t 0 -b 1');
            [predict_label,accuracy_test, prob_estimates(c,Idx2clf(n),:,:)] = svmpredict(testlabel,testset1,model,'-b 1');

            Acc.mean(c,Idx2clf(n)) = accuracy_test(1);
            Acc.predictlabel{c,Idx2clf(n)} = predict_label;

            Acc.Thumb(c,Idx2clf(n)) = sum(predict_label(testlabel == 1) == testlabel(testlabel == 1)) / numel(testlabel(testlabel == 1));
            Acc.Index(c,Idx2clf(n)) = sum(predict_label(testlabel == 2) == testlabel(testlabel == 2)) / numel(testlabel(testlabel == 2));
            Acc.Middle(c,Idx2clf(n)) = sum(predict_label(testlabel == 3) == testlabel(testlabel == 3)) / numel(testlabel(testlabel == 3));
            Acc.Little(c,Idx2clf(n)) = sum(predict_label(testlabel == 4) == testlabel(testlabel == 4)) / numel(testlabel(testlabel == 4));

        end

    end

end


%% call mean Acc 5-fold per task
Acc_alltrial.mean = zeros(numel(Timeseg),numel(ROIidx));
Acc_alltrial.Thumb = zeros(numel(Timeseg),numel(ROIidx));
Acc_alltrial.Index = zeros(numel(Timeseg),numel(ROIidx));
Acc_alltrial.Middle = zeros(numel(Timeseg),numel(ROIidx));
Acc_alltrial.Little = zeros(numel(Timeseg),numel(ROIidx));

for k = 1:5
    Acc_alltrial.mean = Acc_alltrial.mean + Acc.mean;
    Acc_alltrial.Thumb = Acc_alltrial.Thumb + Acc.Thumb;
    Acc_alltrial.Index = Acc_alltrial.Index+ Acc.Index;
    Acc_alltrial.Middle = Acc_alltrial.Middle + Acc.Middle;
    Acc_alltrial.Little = Acc_alltrial.Little + Acc.Little;
end

Acc_alltrial.mean = Acc_alltrial.mean/5;
Acc_alltrial.Thumb = Acc_alltrial.Thumb/5;
Acc_alltrial.Index = Acc_alltrial.Index/5;
Acc_alltrial.Middle = Acc_alltrial.Middle/5;
Acc_alltrial.Little = Acc_alltrial.Little/5;

%% cal Task label for each vertex in ROI
for c = 1:numel(Timeseg)
    for t = 1:4
        eval(['Acc2com(t,:) = Acc_alltrial.' opt.task{t} '(c,:);']);
    end

    for i = 1:length(Acc2com(1,:))
        if Acc_alltrial.mean(c,i) > 25
            if find(Acc2com(:,i) ~= 0)>=1
                [~, tasklabel(c,i)] = max(Acc2com(:,i));
            else
                tasklabel(c,i) = 0;
            end
        else
            tasklabel(c,i) = 0;
        end
    end
end

Map_alltask_ROI = tasklabel;
Map_alltask_wholebrain = nan(numel(Timeseg),numel(sourcemodel.thickness));
Map_alltask_wholebrain(:,ROIidx) = Map_alltask_ROI;

%% group-level ACC & Label 
for ss = 1:16
    Acc_mean_allsub(ss,:,:) = zeros(numel(Timeseg), size(sourcemodel.pos,1));
    Acc_task_allsub(ss,:,1,:) = zeros(numel(Timeseg), size(sourcemodel.pos,1));

    Acc_mean_allsub(ss,:,ROIidx) = Acc_alltrial.mean;
    Acc_task_allsub(ss,:,1,ROIidx) = Acc_alltrial.Thumb;
    Acc_task_allsub(ss,:,2,ROIidx) = Acc_alltrial.Index;
    Acc_task_allsub(ss,:,3,ROIidx) = Acc_alltrial.Middle;
    Acc_task_allsub(ss,:,4,ROIidx) = Acc_alltrial.Little;

    Map_alltask_allsub(ss,:,:) = Map_alltask_wholebrain;
end


%%  Nonparameters Statistic methods - fdr method
% need fieldtrip source structure & MNI template sourcemodel
Source2stat.dimord = 'subj_time_pos';
Base2stat = Source2stat;  % chance level: 25%
Base2stat.dimord = 'subj_time_pos';

Acc2inside = Acc_mean_allsub;
Acc2inside(:,:,setdiff(1:15684,ROIidx)) = 0;

Source2stat.time = Timeseg;
Base2stat.time = Timeseg;
for ss = 1:16
    Source2stat.acc(ss,:,:) = squeeze(Acc2inside(ss,:,ROIidx));
    Base2stat.acc(ss,:,:) = zeros(size(Source2stat.acc(ss,:,:)));
    Base2stat.acc(ss,:,:) = 25;
end

cfg = [];
cfg.parameter = 'acc';
cfg.latency = 'all';
cfg.tri = sourcemodel.tri;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'fdr';
cfg.tail = 1;
cfg.alpha = 0.05;
cfg.numrandomization = 5000;

Nsubj = 16;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1, Nsubj) ones(1, Nsubj)*2];

cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;

Source_stat_nps = ft_sourcestatistics(cfg, Source2stat, Base2stat);


%% generate Digit-map  -  DEVOTE
for ss = 1:16
    for t = 1:4
        if ss == 1
            Map_alltask_group{t} = zeros(numel(Timeseg),numel(sourcemodel.thickness));
        end

        map_ind = [];
        map_ind = Map_alltask_wholebrain_perindiv;
        map_ind(map_ind ~= t) = 0;
        map_ind(map_ind == t) = 1;

        Map_alltask_group{t} = Map_alltask_group{t} + map_ind;
    end
end

for c = 1:numel(Timeseg)
    for t = 1:4
        Map2com(t,:) = Map_alltask_group{t}(c,:);
    end

    for i = 1:length(Map2com(1,:))
        if sum(Map2com(:,i)) ~= 0
            [~,Map_4taskin1brain(c,i)] = max(Map2com(:,i));
        else
            Map_4taskin1brain(c,i) = 0;
        end
    end

end





