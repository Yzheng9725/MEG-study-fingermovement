% EMG Classification
% 4 classes

%% flit data to freqband
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 20;
cfg.lpfilter = 'yes';
cfg.lpfreq = 450;
cfg.bsfilter = 'yes';
cfg.bsfreq = [49.5 50.5; 99 101; 147 153; 196 204];
data2ana = cell(1, numel(opt.task));
for t = 1:numel(opt.task)
    data2ana{t} = ft_preprocessing(cfg, EMGdata_task{t});
end

%% append data
Sensordata = [];
for t = 1:numel(opt.task)
    if t == 1
        Sensordata = data2ana{t};
    else
        Sensordata = ft_appenddata([], Sensordata, data2ana{t});
    end
end

%% classification
% cross val
clfopt.trialnum = 100;
clfopt.kfold = 5;
clfopt.Label_all = [ones(1,clfopt.trialnum)*1 ones(1,clfopt.trialnum)*2 ones(1,clfopt.trialnum)*3 ones(1,clfopt.trialnum)*4]';
clfopt.ClassIdx_all = crossvalind('Kfold',clfopt.Label_all,clfopt.kfold);
clfopt.TaskLabel = {'Thumb', 'Index', 'Middle','Little'};

% 5-kfold - train/test set
for k = 1:5
    cfg = [];
    cfg.trials = find(clfopt.ClassIdx_all ~= k);
    train2feature = ft_selectdata(cfg, Sensordata);

    cfg.trials = find(clfopt.ClassIdx_all == k);
    test2feature = ft_selectdata(cfg, Sensordata);

    DATA2train = train2feature.trial;
    DATA2test = test2feature.trial;

    for i = 1:numel(DATA2train)
        trainset(i,:) = reshape(DATA2train{i}',1,[]);
        trainset(i,:) = abs(trainset(i,:));
        trainset(i,:)  = ft_preproc_lowpassfilter(trainset(i,:), 1000, 50, 6);
    end

    for i = 1:numel(DATA2test)
        testset(i,:) = reshape(DATA2test{i}',1,[]);
        testset(i,:) = abs(testset(i,:));
        testset(i,:)  = ft_preproc_lowpassfilter(testset(i,:), 1000, 50, 6);
    end

    trainlabel = train2feature.trialinfo;
    testlabel = test2feature.trialinfo;
end

%% SVM
% normalization
[~,PStrain] = mapminmax(trainset');
PStrain.ymin = -1;
PStrain.ymax = 1;
[trainset1,PStrain] = mapminmax(trainset',PStrain);
trainset1 = trainset1';
[testset1,PStest] = mapminmax('apply', testset', PStrain);
testset1 = testset1';

% OVA - SVM
numClasses = numel(opt.task);
models = cell(numClasses, 1);

% train model
for i = 1:numClasses
    labels = trainlabel;
    labels(labels ~= i) = -1;
    labels(labels == i) = 1;
    models{i} = svmtrain(labels, trainset1, '-t 0');
end

% predict
scores = zeros(numel(testlabel), numClasses);
for i = 1:numClasses
    tlabels = testlabel;
    tlabels(tlabels ~= i) = -1;
    tlabels(tlabels == i) = 1;
    [predict_label, accuracy_test, s] = svmpredict(tlabels, testset1, models{i}, '-b 1');
    scores(:, i) = s(:, max(models{i}.Label));
end

[~, predictions] = max(scores, [], 2);

Acc_ova.test = sum(predictions == testlabel) / numel(testlabel);
Acc_ova.predictlabel = predictions ;

Acc_ova.Thumb = sum(predictions(testlabel == 1) == testlabel(testlabel == 1)) / numel(testlabel(testlabel == 1));
Acc_ova.Index = sum(predictions(testlabel == 2) == testlabel(testlabel == 2)) / numel(testlabel(testlabel == 2));
Acc_ova.Middle = sum(predictions(testlabel == 3) == testlabel(testlabel == 3)) / numel(testlabel(testlabel == 3));
Acc_ova.Little= sum(predictions(testlabel == 4) == testlabel(testlabel == 4)) / numel(testlabel(testlabel == 4));