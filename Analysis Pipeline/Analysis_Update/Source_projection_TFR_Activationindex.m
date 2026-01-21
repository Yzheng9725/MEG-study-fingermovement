% Source projection & TFR

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

%% Common filter
data_alltask = ft_appenddata([], data2ana{1}, data2ana{2}, data2ana{3}, data2ana{4});

cfg = [];
cfg.covariance = 'yes';
data4cov = ft_timelockanalysis(cfg, data_alltask);

cfg = [];
cfg.method          = 'lcmv';
cfg.headmodel       = headmodel;
cfg.sourcemodel     = leadfield;
cfg.unit            = sourcemodel.unit;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori   = 'yes';
source_lcmv         = ft_sourceanalysis(cfg, data4cov);
commonfilter_wb = source_lcmv.avg.filter;

%% Source reconstruction - lcmv
% select vertices in ROI
ROI2use = [ ];  % Precentral & postcentral index in parcellationLabel
ROIidx = find(brainatlas.parcellation == ROI2use);
commonfilter_ROI = cell2mat(commonfilter_wb(ROIidx));

Sourcedata = [];
Sourcedata.time(1:length(data2ana_all.trial)) = data2ana_all.time{1};
for trl = 1:length(data2ana_all.trial)
    Sourcedata.trial{trl} = commonfilter_ROI * data2ana_all.trial{trl};
end
Sourcedata.trialinfo = data2ana_all.trialinfo;

%% Time-frequency Response
data2process = ft_appenddata([], Sourcedata_lcmv{1}, Sourcedata_lcmv{2}, ...
    Sourcedata_lcmv{3}, Sourcedata_lcmv{4});
data2freq = data2process;
% padding
data2freq.trial = cellfun(@(x) ft_preproc_padding(x,'zero',(7-2.5)*data2freq.fsample,(7-2.5)*data2freq.fsample), ...
    data2freq.trial,'UniformOutput',false);
data2freq.time = cellfun(@(x) -7:1/data2freq.fsample:7-1/data2freq.fsample, data2freq.time,'UniformOutput',false);

cfg                             = [];
cfg.method                      = 'wavelet';
cfg.output                      = 'pow';
cfg.foi                         = [0.5:0.5:8 9:1:90];
cfg.toi                         = -2.5:0.1:2.5;
cfg.width                       = 4;  
cfg.keeptrials                  = 'yes'; 
TFWAVE_task                     = ft_freqanalysis(cfg, data2freq);

%% Activation Index - epoch
basewindow = [-2 -1.5];
Actseg = [-1 2];

ActIndex = cell(1, numel(opt.task));
for t = 1:numel(opt.task)
    for trl = 1:numel(data2ana{t}.trial)
        % cov - baseline
        Time4base = knnsearch(data2ana{t}.time{1}', basewindow(1)):knnsearch(data2ana{t}.time{1}', basewindow(2));
        data2base = data2ana{t}.trial{trl}(:, Time4base);
        data2base  = ft_preproc_baselinecorrect(data2base);
        covb = (data2base * data2base') ./ length(Time4base);

        % cov - active
        Time4act = knnsearch(data2ana{t}.time{1}', Actseg(1)):knnsearch(data2ana{t}.time{1}', Actseg(2));
        data2act = data2ana{t}.trial{trl}(:, Time4act);
        data2act  = ft_preproc_baselinecorrect(data2act);
        cova = (data2act * data2act') ./ length(Time4act);

        %% cal A
        for i = 1:numel(ROIidx)
            ActIndex{t}(trl,i) = (Commonfilter(i,:)*cova*Commonfilter(i,:)' - Commonfilter(i,:)*covb*Commonfilter(i,:)')...
                /(Commonfilter(i,:)*cova*Commonfilter(i,:)' + Commonfilter(i,:)*covb*Commonfilter(i,:)');
        end
    end
end

%% Activation Index - window
Timewin = Actseg(1):0.05:Actseg(2);
ActIndex = cell(numel(opt.task), numel(Timewin));

for t = 1:numel(opt.task)
    for trl = 1:numel(data2ana{t}.trial)
        % cov - baseline
        Time4base = knnsearch(data2ana{t}.time{1}', basewindow(1)):knnsearch(data2ana{t}.time{1}', basewindow(2));
        data2base = data2ana{t}.trial{trl}(:, Time4base);
        data2base  = ft_preproc_baselinecorrect(data2base);
        covb = (data2base * data2base') ./ length(Time4base);

        %% cal A
        for c = 1:numel(Timewin)
            % cov - active
            time2cal = [Timewin(c) Timewin(c)+0.05];
            Time4act = knnsearch(data2ana{t}.time{1}', time2cal(1)):knnsearch(data2ana{t}.time{1}', time2cal(2));
            data2act = data2ana{t}.trial{trl}(:, Time4act);
            data2act = ft_preproc_baselinecorrect(data2act);
            cova = (data2act * data2act') ./ length(Time4act);

            for i = 1:numel(ROIidx)
                ActIndex{t,c}(trl,i) = (Commonfilter(i,:)*cova*Commonfilter(i,:)' - Commonfilter(i,:)*covb*Commonfilter(i,:)')...
                    /(Commonfilter(i,:)*cova*Commonfilter(i,:)' + Commonfilter(i,:)*covb*Commonfilter(i,:)');
            end
        end
    end
end
