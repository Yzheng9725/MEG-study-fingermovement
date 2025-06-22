% Finger extension movement
% MEG - Classification - Alpha/Beta/High-Gamma
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

FreqName = {'Alpha','Beta','High-Gamma'};
FreqBand = {[8 13],[13 30],[60 90]};

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

%%
for r = 1:numel(FreqName)
    opt.freqname = FreqName{r};
    opt.freqband = FreqBand{r};


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

        %% cal source pos & filter per task per toi
        for t = 1:numel(opt.task)
            disp(['-------- ' opt.task{t} ' --------']);
            ClassIdx4task = clfopt.ClassIdx_all(find(clfopt.Label_all == t));

            for c = 1:size(clfopt.Time4peak,2)
                ERDseg = [clfopt.Time4peak(t,1) clfopt.Time4peak(t,1)+0.5];

                %% load raw trial data
                eval(['data2source = MEGDATA_' opt.task{t} ';'])

                % cfg = [];
                % cfg.demean        = 'yes';
                % cfg.baselinewindow = baseline;
                % data2source = ft_preprocessing(cfg, data2source);

                cfg = [];
                cfg.trials = find(ClassIdx4task ~= k);
                data2source = ft_selectdata(cfg, data2source);

                %% Baseline
                cfg = [];
                cfg.toilim = baseline;
                dataPre = ft_redefinetrial(cfg, data2source);

                %% ERD
                cfg = [];
                cfg.toilim =  ERDSseg;
                dataPost = ft_redefinetrial(cfg, data2source);

                dataAll = ft_appenddata([], dataPre, dataPost);

                cfg             = [];
                cfg.method      = 'mtmfft';
                cfg.output      = 'powandcsd';
                cfg.taper       = 'hanning';
                cfg.channel     = 'meg';
                cfg.tapsmofrq   = 1;
                cfg.foilim      = opt.freqband;
                cfg.pad         = 'nextpow2';
                cfg.keeptrials  = 'yes';
                freqAll         = ft_freqanalysis(cfg, dataAll);
                freqPre        = ft_freqanalysis(cfg, dataPre);
                freqPost          = ft_freqanalysis(cfg, dataPost);

                cfg             = [];
                cfg.method      = 'dics';
                cfg.channel     = 'meg';
                cfg.sourcemodel = leadfield_ROI;
                cfg.headmodel   = headmodel;
                cfg.dics.projectnoise = 'yes';
                cfg.dics.lambda = '5%';
                cfg.dics.keepfilter = 'yes';
                cfg.dics.realfilter = 'yes';
                cfg.dics.fixedori  ='yes';
                SourceAll      = ft_sourceanalysis(cfg, freqAll);
                clearvars freqAll

                cfg.sourcemodel.filter = SourceAll.avg.filter;
                SourcePre      = ft_sourceanalysis(cfg, freqPre);
                SourcePost      = ft_sourceanalysis(cfg, freqPost);
                clearvars freqPost

                SourceDiff = SourcePost;
                SourceDiff.avg.pow = (SourcePost.avg.pow-SourcePre.avg.pow)./ (2*SourcePre.avg.pow);

                [MPseudoT, MPosidx] = min(SourceDiff.avg.pow);
                if ~ismember(MPosidx,ROIidx)
                    [MPseudoT, ROIPosidx] = min(SourceDiff.avg.pow(ROIidx));
                    MPosidx = ROIidx(ROIPosidx);
                end

                Pos4Clf{t}(c,:) = sourcemodel.pos(MPosidx,:);
                Posidx4Clf(t,c) = MPosidx;
                Ori4clf{t,c} = nan(numel(SourceAll.avg.ori),3);
                for i = 1:numel(ROIidx)
                    Ori4clf{t,c}(ROIidx(i),:) = SourceAll.avg.ori{ROIidx(i)};
                end
            end

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
        cfg.lpfreq = opt.freqband(2);
        cfg.hpfilter = 'yes';
        cfg.hpfreq = opt.freqband(1);
        % cfg.demean        = 'yes';
        % cfg.baselinewindow = baseline;
        data2train_recon = ft_preprocessing(cfg, Data2cov);
        data2test_recon = ft_preprocessing(cfg, Data2test);

        cfg = [];
        cfg.latency = actseg;
        data4cov = ft_selectdata(cfg, data2train_recon);

        Cact = cellfun(@(x) (x-mean(x,2))*(x-mean(x,2))', data4cov.trial, 'UniformOutput', false);
        Cact = mean(cat(3, Cact{:}), 3);

        [U,S,V] = svd(Cact,0);
        s = diag(S);
        [m, n] = size(Cact);
        tolerance = 10 * max(m,n) * eps;
        kappa     = sum(s./s(1) > tolerance);
        % lambda    = 0;
        lambda = 0.05 * mean(s);  % 5% lambda
        S = diag(1 ./ (s(1:kappa) + lambda));
        Cinv = V(:,1:kappa)*S*U(:,1:kappa)';

        %% 要预先生成ReconData结构
        ERDSName = {'ERD','ERS'};
        % ERDSName = 'ERS';  % for freqband "High-Gamma"

        Filter4cal = [];
        ReconData = [];
        ReconData.pos = [];
        ReconData.ori = [];
        ReconData.label = [];
        ReconData.time = data2train_recon.time;
        Labelidx = 1;
        Filter4cal = [];

        %% adjust source orientation on the cortex
        vertices = sourcemodel.pos;
        tri = sourcemodel.tri;

        for l = 1:numel(ERDSName)
            for t = 1:numel(opt.task)
                [f,~] = find(tri == Posidx4Clf(t,l));
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

                dotv = dot(avgNormal,Ori4clf{t,l}(Posidx4Clf(t,l),:));

                if dotv >= 0
                    ori4cal{t,l} = Ori4clf{t,l}(Posidx4Clf(t,l),:);
                else
                    ori4cal{t,l} = -1*Ori4clf{t,l}(Posidx4Clf(t,l),:);
                end
            end
        end

        %% recondata cal filt
        for c = 1:numel(ERDSName)
            for t = 1:numel(opt.task)
                ReconData.label{Labelidx,1} = [opt.task{t} '_' ERDSName{c}];
                Labelidx = Labelidx + 1;
                ReconData.ori = [ReconData.ori; Ori4clf{t,c}(Posidx4Clf(t,c),:)];
                ReconData.pos = [ReconData.pos; Pos4Clf{t}(c,:)];

                %% recon - cal filt
                ori4cal{t,c} = Ori4clf{t,c}(Posidx4Clf(t,c),:);
                gain = leadfield.leadfield{Posidx4Clf(t,c)} * ori4cal{t,c}';
                trgain_Cinv = gain' * Cinv;
                filt  = trgain_Cinv / (trgain_Cinv * gain);
                Filter4cal = [Filter4cal; filt];
            end
        end


        %% recon
        Recon2test = ReconData;
        Recon2test.time = data2test_recon.time;

        Recon2train = ReconData;

        % Train - recon
        for i = 1:numel(data2train_recon.trial)
            Recon2train.trial{i} = Filter4cal * data2train_recon.trial{i};
        end

        % Test - recon
        for i = 1:numel(data2test_recon.trial)
            Recon2test.trial{i} = Filter4cal * data2test_recon.trial{i};
        end

        cfg = [];
        cfg.resamplefs = 200;
        Recon2test = ft_resampledata(cfg, Recon2test);
        Recon2train = ft_resampledata(cfg, Recon2train);

        %%
        cfg = [];
        Train2TF = ft_preprocessing(cfg,Recon2train);
        Test2TF = ft_preprocessing(cfg,Recon2test);

        % Time Freq - keeptrials
        cfg                             = [];
        cfg.method                      = 'wavelet';
        cfg.output                      = 'pow';
        cfg.foi                         = opt.freqband(1):0.5:opt.freqband(2);
        cfg.toi                         = actseg(1):0.1:actseg(2);
        cfg.width                       = 7;
        cfg.keeptrials                  = 'yes'; % calculate all trials and keep each results
        TFWAVE2train = ft_freqanalysis(cfg, Train2TF);
        TFWAVE2test = ft_freqanalysis(cfg, Test2TF);



        %% SVM
        TimeSel = {'ERD','ERS'};
        Chan4Sel = [];

        n = 1;
        for l = 1:numel(opt.task)
            for i = 1:numel(TimeSel)
                Chan4Sel{n} = [opt.task{l} '_' TimeSel{i}];
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
        cfg.latency = actseg;
        train2feature = ft_selectdata(cfg, TFWAVE2train);
        test2feature = ft_selectdata(cfg, TFWAVE2test);

        DATA2train = squeeze(mean(train2feature.powspctrm(:,:,:,:),3));
        DATA2test = squeeze(mean(test2feature.powspctrm(:,:,:,:),3));

        for i = 1:size(DATA2train,1)
            trainset(i,:) = reshape(squeeze(DATA2train(i,:,:))',1,[]);
        end

        for i = 1:size(DATA2test,1)
            testset(i,:) = reshape(squeeze(DATA2test(i,:,:))',1,[]);
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

    

