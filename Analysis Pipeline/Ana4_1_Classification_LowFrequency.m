% Finger extension movement
% MEG - Classification - Low-frequency - add temporal features 
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


%% Classify option set
clfopt.trialnum = 100;
clfopt.kfold = 5;
clfopt.Label_all = [ones(1,clfopt.trialnum)*1 ones(1,clfopt.trialnum)*2 ones(1,clfopt.trialnum)*3 ones(1,clfopt.trialnum)*4]';
clfopt.ClassIdx_all = crossvalind('Kfold',clfopt.Label_all,clfopt.kfold);
clfopt.TaskLabel = {'Thumb', 'Index', 'Middle','Little'};
clfopt.Time4peak = Time2max;


%% Cross Validation
for k = 1:clfopt.kfold
    %% change leadfield.inside 2 keep ROI
    leadfield4ROI  = leadfield;
    leadfield4ROI.inside = false(size(leadfield.inside));
    leadfield4ROI.inside(ROIidx) = true;
    sourcemodel4ROI = sourcemodel;
    sourcemodel4ROI.inside = leadfield4ROI.inside;

    Timestep = 0.02;

    %% cal source pos & filter per task per toi
    for t = 1:numel(opt.task)
        disp(['-------- ' opt.task{t} ' --------']);
        ClassIdx4task = clfopt.ClassIdx_all(find(clfopt.Label_all == t));

        %% load raw trial data
        eval(['data2source = MEGDATA_' opt.task{t} ';'])

        cfg = [];
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 8;
        % cfg.demean        = 'yes';
        % cfg.baselinewindow = baseline;
        data2sam = ft_preprocessing(cfg, data2source);

        %% cal ER-SAM  - using 4 folds
        % use baseline to cal noise_avg/noise_avg.cov
        cfg = [];
        cfg.toilim = baseline;
        cfg.trials = find(ClassIdx4task ~= k);
        dataNoise_sam = ft_redefinetrial(cfg, data2sam);

        cfg = [];
        cfg.covariance = 'yes';
        avgNosie_sam = ft_timelockanalysis(cfg, dataNoise_sam);

        cfg = [];
        cfg.toilim = epoch;
        cfg.trials = find(ClassIdx4task ~= k);
        data_sam = ft_redefinetrial(cfg, data2sam);

        cfg = [];
        cfg.covariance = 'yes';
        avg_sam = ft_timelockanalysis(cfg, data_sam);

        %% condition "four peak periods" - 20 ms time window
        Time4peak = clfopt.Time4peak(t,:);
        
        for c = 1:numel(Time4peak)
            Time2cal = knnsearch(avg_sam.time', Time4peak(c));
            Source_sam_avg = ft_inverse_sam(leadfield4ROI, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
                'reducerank',2, 'toi', [Time2cal Time2cal+Timestep*data_sam.fsample], 'noisecov', avgNosie_sam.cov);

            [MaxPseudoZ, MaxPosidx] = max(Source_sam_avg.pseudoZ);

            Pos4Clf{t}(c,:) = sourcemodel.pos(MaxPosidx,:);
            Posidx4Clf(t,c) = MaxPosidx;
            Ori4clf{t,c} = nan(numel(Source_sam_avg.ori),3);
            for i = 1:numel(ROIidx)
                Ori4clf{t,c}(ROIidx(i),:) = Source_sam_avg.ori{ROIidx(i)};
            end
        end

        %% condition "Active" - [-0.5 1.0]
        % Time2cal = [knnsearch(avg_sam.time', -0.5) knnsearch(avg_sam.time', 1.0)];
        % Source_sam_avg = ft_inverse_sam(leadfield4ROI, sens, headmodel, avg_sam.avg, avg_sam.cov, ...
        %     'reducerank',2, 'toi', Time2cal, 'noisecov', avgNosie_sam.cov);
        %
        % [MaxPseudoZ, MaxPosidx] = max(Source_sam_avg.pseudoZ);
        %
        % Pos4Clf{t}(1,:) = sourcemodel.pos(MaxPosidx,:);
        % Posidx4Clf(t,1) = MaxPosidx;
        % Ori4clf{t,1} = nan(numel(Source_sam_avg.ori),3);
        % for i = 1:numel(ROIidx)
        %     Ori4clf{t,1}(ROIidx(i),:) = Source_sam_avg.ori{ROIidx(i)};
        % end

        %% define train / test set
        cfg = [];
        cfg.trials = find(ClassIdx4task ~= k);
        data2train = ft_selectdata(cfg, data2source);

        cfg = [];
        cfg.trials = find(ClassIdx4task == k);
        data2test = ft_selectdata(cfg, data2source);

        if t == 1
            Data2cov = ft_appenddata([], data2train);
            Data2test = ft_appenddata([], data2test);
        else
            Data2cov = ft_appenddata([], Data2cov, data2train);
            Data2test = ft_appenddata([], Data2test,data2test);
        end

    end


    %% Source reconstruction - using all folds
    %% cal cov
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    % cfg.demean        = 'yes';
    % cfg.baselinewindow = baseline;
    data2avg_recon = ft_preprocessing(cfg, Data2cov);
    data2test_recon = ft_preprocessing(cfg, Data2test);

    cfg = [];
    cfg.covariance = 'yes';
    cfg.covariancewindow = epoch;
    cfg.latency = epoch;
    avg2cov = ft_timelockanalysis(cfg, data2avg_recon);


    %% create ReconData structure
    TimeName = {'t1','t2','t3','t4'};
    % TimeName = 'Active';
    Filter4cal = [];
    ReconData = [];
    ReconData.pos = [];
    ReconData.ori = [];
    ReconData.label = [];
    ReconData.time = data2avg_recon.time;
    Labelidx = 1;

    for c = 1:numel(TimeName)
        for t = 1:numel(opt.task)
            ReconData.label{Labelidx,1} = [opt.task{t} '_' TimeName{c}];  % 对应filter顺序：task+time
            Labelidx = Labelidx + 1;
            ReconData.ori = [ReconData.ori; Ori4clf{t,c}(Posidx4Clf(t,c),:)];
            ReconData.pos = [ReconData.pos; Pos4Clf{t}(c,:)];
        end
    end

    Recon2test = ReconData;
    Recon2test.time = data2test_recon.time;

    Filter4cal = [];

    %% cal space(inverse) filter
    vertices = sourcemodel.pos;
    tri = sourcemodel.tri;

    for c = 1:numel(TimeName)
        for t = 1:numel(opt.task)
            %% adjust source orientation on the cortex
            [f,~] = find(tri == Posidx4Clf(t,c));
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

            dotv = dot(avgNormal,Ori4clf{t,c}(Posidx4Clf(t,c),:));

            if dotv >= 0
                ori4cal{t,c} = Ori4clf{t,c}(Posidx4Clf(t,c),:);
            else
                ori4cal{t,c} = -1*Ori4clf{t,c}(Posidx4Clf(t,c),:);
            end

            %% recon - cov inverse
            Cact = avg2cov.cov;
            [U,S,V] = svd(Cact,0);
            s = diag(S);
            [m, n] = size(Cact);
            tolerance = 10 * max(m,n) * eps;
            kappa     = sum(s./s(1) > tolerance);
            lambda    = 0;
            % lambda = 0.05 * mean(s);  % 5% lambda
            S = diag(1 ./ (s(1:kappa) + lambda));
            Cinv = V(:,1:kappa)*S*U(:,1:kappa)';

            %% recon - cal filt
            gain = leadfield.leadfield{Posidx4Clf(t,c)} * ori4cal{t,c}';
            trgain_Cinv = gain' * Cinv;
            filt  = trgain_Cinv / (trgain_Cinv * gain);
            Filter4cal = [Filter4cal; filt];

        end
    end

    %% recon
    % Train - recon
    for i = 1:numel(data2avg_recon.trial)
        ReconData.trial{i} = Filter4cal * data2avg_recon.trial{i};
    end

    ReconData.trainlabel = clfopt.Label_all(clfopt.ClassIdx_all ~= k);
    Recon2train = ReconData;

    % Test - recon
    for i = 1:numel(data2test_recon.trial)
        Recon2test.trial{i} = Filter4cal * data2test_recon.trial{i};
    end

    Recon2test.testlabel = clfopt.Label_all(clfopt.ClassIdx_all == k);



    %% SVM
    TimeSel = {'t1','t2','t3','t4'};
    TimeIdx = {1,2,3,4,[1 2],[1 3],[1 4],[2 3],[2 4],[3 4],[1 2 3],[1 2 4],[1 3 4],[2 3 4],[1 2 3 4]};

    % TimeSel = 'Active';  % for condition "Active"
    % TimeIdx = 1;
    
    for m = 1:numel(TimeIdx)
        Chan4Sel = [];
        TimeIdx2cal = TimeIdx{m};

        n = 1;
        for t = 1:numel(opt.task)
            for i = 1:numel(TimeIdx2cal)
                Chan4Sel{n} = [opt.task{t} '_' TimeSel{TimeIdx2cal(i)}];
                n = n+1;
            end
        end

        % train/test data set
        DATA2train = [];
        DATA2test = [];
        trainset = [];
        testset =[];

        cfg = [];
        cfg.channel = Chan4Sel;
        cfg.latency = [-0.5 1.0];
        train2feature = ft_selectdata(cfg, Recon2train);
        test2feature = ft_selectdata(cfg, Recon2test);

        cfg = [];
        cfg.resamplefs = 16;
        train2feature = ft_resampledata(cfg, train2feature);
        test2feature = ft_resampledata(cfg, test2feature);

        DATA2train = train2feature.trial;
        DATA2test = test2feature.trial;

        for i = 1:numel(DATA2train)
            trainset(i,:) = reshape(DATA2train{i},1,[]);
        end

        for i = 1:numel(DATA2test)
            testset(i,:) = reshape(DATA2test{i},1,[]);
        end

        trainlabel = clfopt.Label_all(clfopt.ClassIdx_all ~= k);
        testlabel = clfopt.Label_all(clfopt.ClassIdx_all == k);

        % normalization
        [~,PStrain] = mapminmax(trainset');
        PStrain.ymin = -1;
        PStrain.ymax = 1;
        [trainset1,PStrain] = mapminmax(trainset',PStrain);
        trainset1 = trainset1';
        testset1 = mapminmax('apply',testset',PStrain);
        testset1 = testset1';

        % OVA
        numClasses = numel(opt.task);
        models = cell(numClasses, 1);

        % train
        for i = 1:numClasses
            labels = trainlabel;
            labels(labels ~= i) = -1;
            labels(labels == i) = 1;
            models{i} = svmtrain(labels, trainset1, '-t 0 -b 1');  % linear
        end

        % predict
        scores = zeros(numel(testlabel), numClasses);

        for i = 1:numClasses
            tlabels = testlabel;
            tlabels(tlabels ~= i) = -1;
            tlabels(tlabels == i) = 1;
            [predict_label{k,i}, accuracy_test{k,i}, s] = svmpredict(tlabels, testset1, models{i}, '-b 1');
            scores(:, i) = s(:, max(models{i}.Label));
        end

        [~, predictions] = max(scores, [], 2);

        Acc_ova.test(k) = sum(predictions == testlabel) / numel(testlabel);
        Acc_ova.predictlabel{k} = predictions ;

        Acc_ova.Thumb(k) = sum(predictions(testlabel == 1) == testlabel(testlabel == 1)) / numel(testlabel(testlabel == 1));
        Acc_ova.Index(k) = sum(predictions(testlabel == 2) == testlabel(testlabel == 2)) / numel(testlabel(testlabel == 2));
        Acc_ova.Middle(k) = sum(predictions(testlabel == 3) == testlabel(testlabel == 3)) / numel(testlabel(testlabel == 3));
        Acc_ova.Little(k) = sum(predictions(testlabel == 4) == testlabel(testlabel == 4)) / numel(testlabel(testlabel == 4));

    end
end
