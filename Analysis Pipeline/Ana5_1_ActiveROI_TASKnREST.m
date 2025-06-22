% Finger extension movement
% MEG - Active ROI
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


%%  %% Source space Projection - SAM
for t = 1: numel(opt.task)
    %% Cross Val
    for k = 1:5
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

        % use baseline to cal noise_avg/noise_avg.cov
        cfg = [];
        cfg.toilim = baseline;
        dataNoise_sam = ft_redefinetrial(cfg, data2sam);

        cfg = [];
        cfg.covariance = 'yes';
        avgNoise_sam = ft_timelockanalysis(cfg, dataNoise_sam);

        % Epoch
        cfg = [];
        cfg.toilim = epoch;
        dataAvg_sam = ft_redefinetrial(cfg, data2sam);

        cfg = [];
        cfg.covariance = 'yes';
        avg_sam = ft_timelockanalysis(cfg, dataAvg_sam);


        %% cal Weights - REST - 20ms in baseline - ERSAM
        Timestep = 0.02;
        Time2cal =  knnsearch(avgNoise_sam.time', -1);

        ori_correct_rest = nan(numel(sourcemodel.thickness),3);
        filter4recon_rest = nan(numel(sourcemodel.thickness),numel(label));

        Source_sam_avg = ft_inverse_sam(leadfield_ROI, sens, headmodel, avgNoise_sam.avg, avgNoise_sam.cov, ...
            'reducerank',2,  'toi', [Time2cal-Timestep*data2sam.fsample  Time2cal], 'noisecov', avgNoise_sam.cov);  % 'toi', [Time2cal-Timestep*data2sam.fsample  Time2cal],

        filter_sam(sourcemodel_ROI.inside, :) = cell2mat(Source_sam_avg.filter');
        ori_sam(sourcemodel_ROI.inside, :) = cell2mat(Source_sam_avg.ori)';

        % adjust ori
        vertices = sourcemodel.pos;
        tri = sourcemodel.tri;
        ori4cal = ori_sam;

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

            dotv = dot(avgNormal,ori4cal(ROIidx(i), :));

            if dotv >= 0
                ori_correct_rest(ROIidx(i), :) = ori_sam(ROIidx(i), :);
                filter4recon_rest(ROIidx(i), :) = filter_sam(ROIidx(i), :);
            else
                ori_correct_rest(ROIidx(i), :) = -1*ori_sam(ROIidx(i), :);
                filter4recon_rest(ROIidx(i), :) = -1*filter_sam(ROIidx(i), :);
            end
        end

        clearvars filter_sam ori_sam Source_sam_avg

        %% cal Weights  - TOI: slip windows -500 - 1000ms  by 20ms step - ERSAM
        Timeseg = [];
        Timeseg = -0.5:Timestep:1-Timestep;

        filter_sam = nan(numel(Timeseg),numel(sourcemodel.thickness),numel(label));
        ori_sam = nan(numel(Timeseg),numel(sourcemodel.thickness),3);
        ori_correct = nan(numel(Timeseg),numel(sourcemodel.thickness),3);
        filter4recon = nan(numel(Timeseg),numel(sourcemodel.thickness),numel(label));

        for c = 1:numel(Timeseg)
            Time2cal = knnsearch(avg_sam.time', Timeseg(c));

            %%
            % cal event-related SAM, using computed leadfield
            Source_sam_avg = ft_inverse_sam(leadfield_ROI, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
                'reducerank',2, 'toi', [Time2cal Time2cal+Timestep*data2sam.fsample],  'noisecov',avgNoise_sam.cov);   %

            filter_sam(c, sourcemodel_ROI.inside, :) = cell2mat(Source_sam_avg.filter');
            ori_sam(c, sourcemodel_ROI.inside, :) = cell2mat(Source_sam_avg.ori)';

            % adjust ori
            vertices = sourcemodel.pos;
            tri = sourcemodel.tri;
            ori4cal = squeeze(ori_sam(c, :, :));

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

                dotv = dot(avgNormal,ori4cal(ROIidx(i), :));

                if dotv >= 0
                    ori_correct(c, ROIidx(i), :) = ori_sam(c, ROIidx(i), :);
                    filter4recon(c, ROIidx(i), :) = filter_sam(c, ROIidx(i), :);
                else
                    ori_correct(c, ROIidx(i), :) = -1*ori_sam(c, ROIidx(i), :);
                    filter4recon(c, ROIidx(i), :) = -1*filter_sam(c, ROIidx(i), :);
                end
            end

        end

        %% Data Recon
        %% estimate ROI-prob MAP for each window
        data2recon = data2sam;

        ReconData_train = [];
        ReconData_test = [];

        cfg = [];
        cfg.toilim = [-1.0-Timestep -1.0];
        cfg.trials = find(ClassIdx4task ~= k);
        data2recon_train = ft_selectdata(cfg, data2recon);

        cfg.trials = find(ClassIdx4task == k);
        data2recon_test = ft_selectdata(cfg, data2recon);

        ReconData_train{2} = nan(numel(data2recon_train.trial),numel(ROIidx),numel(data2recon_train.time{1}));
        ReconData_test{2} = nan(numel(data2recon_test.trial),numel(ROIidx),numel(data2recon_test.time{1}));

        % rest
        Filter4cal = squeeze(filter4recon_rest(:, :));
        Filter4cal(isnan(Filter4cal(:,1)),:) = [];

        for i = 1:numel(data2recon_train.trial)
            ReconData_train{2}(i,:,:) = Filter4cal * data2recon_train.trial{i};
        end

        for i = 1:numel(data2recon_test.trial)
            ReconData_test{2}(i,:,:) = Filter4cal * data2recon_test.trial{i};
        end

        % task
        for c = 1:numel(Timeseg)
            Filter4cal = squeeze(filter4recon(c, :, :));
            Filter4cal(isnan(Filter4cal(:,1)),:) = [];

            % train/test set
            cfg = [];
            cfg.toilim = [Timeseg(c) Timeseg(c)+Timestep];
            cfg.trials = find(ClassIdx4task ~= k);
            data2recon_train = ft_selectdata(cfg, data2recon);

            cfg.trials = find(ClassIdx4task == k);
            data2recon_test = ft_selectdata(cfg, data2recon);

            %%
            ReconData_train{1} = nan(numel(data2recon_train.trial),numel(ROIidx),numel(data2recon_train.time{1}));
            ReconData_test{1} = nan(numel(data2recon_test.trial),numel(ROIidx),numel(data2recon_test.time{1}));

            for i = 1:numel(data2recon_train.trial)
                ReconData_train{1}(i,:,:) = Filter4cal * data2recon_train.trial{i};
            end

            for i = 1:numel(data2recon_test.trial)
                ReconData_test{1}(i,:,:) = Filter4cal * data2recon_test.trial{i};
            end

            %% Classification
            trainlabel = [ones(size(ReconData_train{1},1),1);zeros(size(ReconData_train{1},1),1)];
            testlabel = [ones(size(ReconData_test{1},1),1);zeros(size(ReconData_test{1},1),1)];

            for n = 1:numel(ROIidx)
                trainset = [];
                testset =[];

                for i = 1:numel(ReconData_train)
                    trainset = [trainset; squeeze(ReconData_train{i}(:,n,:))];
                    testset = [testset; squeeze(ReconData_test{i}(:,n,:))];
                end

                %% normalization
                [~,PStrain] = mapminmax(trainset');
                PStrain.ymin = -1;
                PStrain.ymax = 1;
                [trainset1,PStrain] = mapminmax(trainset',PStrain);
                trainset1 = trainset1';
                testset1 = mapminmax('apply',testset',PStrain);
                testset1 = testset1';

                %% SVM
                model = svmtrain(trainlabel, trainset1, '-t 0 -b 1');
                [predict_tslabel, accuracy_test, tss] = svmpredict(testlabel, testset1, model);
                Acc(c,n) = accuracy_test(1);

            end
        end
        %% cal sum(ACCperfold)
        if k == 1
            AccAlltrial = zeros(size(Acc));
        end
        AccAlltrial = AccAlltrial+Acc;
    end
    %% cal mean(ACC)
    AccAlltrial = AccAlltrial/5;

end



%% Group-lavel
for ss = 1:16   % 16 participants
    AccAlltask(ss,:,:) = zeros(75, size(sourcemodel.pos,1));
    for t = 1:numel(opt.task)
        AccAlltask(ss,:,ROIidx) = squeeze(AccAlltask(ss,:,ROIidx))+AccAlltrial;
    end

    AccAlltask(ss,:,:) = AccAlltask(ss,:,:)/4;
end

%%  Nonparameters Statistic methods - fdr method
% need fieldtrip source structure & MNI template sourcemodel
Source2stat.dimord = 'subj_time_pos';
Base2stat = Source2stat;  % chance level: 50%
Base2stat.dimord = 'subj_time_pos';

Acc2inside = AccAlltask;
Acc2inside(:,:,setdiff(1:15684,ROIidx)) = 0;

Source2stat.time = Timeseg;
Base2stat.time = Timeseg;
for ss = 1:16
    Source2stat.acc(ss,:,:) = squeeze(Acc2inside(ss,:,ROIidx));
    Base2stat.acc(ss,:,:) = zeros(size(Source2stat.acc(ss,:,:)));
    Base2stat.acc(ss,:,:) = 50;
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


