% Finger extension movement
% MEG - Sourcelevel - Alpha/Beta/High-Gamma - DICS
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

    %% DICS-Source - 500 ms sliping window, 100 ms step
    for t = 1:numel(opt.task)

        Timeseg = [];
        Timewindowlength = 0.5;
        Timeseg = -0.5:0.1:2-Timewindowlength;  
        
        cfg = [];
        cfg.toilim = baseline;
        dataPre = ft_redefinetrial(cfg, data2ana);

        for c = 1:length(Timeseg)
            Time2cal = [Timeseg(c) Timeseg(c)+Timewindowlength];

            cfg = [];
            cfg.toilim =  Time2cal;
            dataPost = ft_redefinetrial(cfg, data2ana);

            dataAll = ft_appenddata([], dataPre, dataPost);

            %% cal fft
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

            %% cal inverse filter
            cfg             = [];
            cfg.method      = 'dics';
            cfg.channel     = 'meg';
            cfg.sourcemodel = leadfield;
            cfg.headmodel   = headmodel;
            cfg.dics.projectnoise = 'yes';
            cfg.dics.lambda = '5%';
            cfg.dics.keepfilter = 'yes';
            cfg.dics.realfilter = 'yes';
            cfg.dics.fixedori  ='yes';
            SourceAll      = ft_sourceanalysis(cfg, freqAll);
            clearvars freqAll

            %% cal source 
            cfg.sourcemodel.filter = SourceAll.avg.filter;
            SourcePre      = ft_sourceanalysis(cfg, freqPre);
            SourcePost      = ft_sourceanalysis(cfg, freqPost);
            clearvars freqPost

            SourceDiff = SourcePost;
            SourceDiff.avg.pow = (SourcePost.avg.pow-SourcePre.avg.pow) ./ (2*SourcePre.avg.pow);

            % ERS - Max pT in ROIidx
            [MaxPseudoT, MaxPosidx] = max(SourceDiff.avg.pow);
            if ~ismember(MaxPosidx,ROIidx)
                [MaxPseudoT, ROIPosidx] = max(SourceDiff.avg.pow(ROIidx));
                MaxPosidx = ROIidx(ROIPosidx);
            end

            Source_dics_ERDS(c).time = Time2cal;
            Source_dics_ERDS(c).ERSpos = sourcemodel.pos(MaxPosidx,:);
            Source_dics_ERDS(c).ERSposidx = MaxPosidx;
            Source_dics_ERDS(c).ERSori = SourceAll.avg.ori{1,MaxPosidx};
            Source_dics_ERDS(c).maxpseudoT = MaxPseudoT;
            Source_dics_ERDS(c).pseudoT = SourceDiff.avg.pow;

            %%
            % ERD - Min pT in ROIidx - only Alpha/Beta
            if r == 1 || r == 2
                [MinPseudoT, MinPosidx] = min(SourceDiff.avg.pow);
                if ~ismember(MinPosidx,ROIidx)
                    [MinPseudoT, ROIPosidx] = min(SourceDiff.avg.pow(ROIidx));
                    MinPosidx = ROIidx(ROIPosidx);
                end

                Source_dics_ERDS(c).ERDpos = sourcemodel.pos(MinPosidx,:);
                Source_dics_ERDS(c).ERDposidx = MinPosidx;
                Source_dics_ERDS(c).ERDori = SourceAll.avg.ori{1,MinPosidx};
                Source_dics_ERDS(c).minpseudoT = MinPseudoT;
            end

        end
    end
end


%% find pZ peak window
% find in group ERD/Stimeseg
ERDtime2sel = [knnsearch(Timeseg', ERDtimeseg(1)) knnsearch(Timeseg', ERDtimeseg(2)-0.5)];  % 0.5 - 500 ms time window length
ERStime2sel = [knnsearch(Timeseg', ERStimeseg(1)) knnsearch(Timeseg', ERStimeseg(2))];

for t = 1 : length(opt.task)
    v2find_ERS(:,:) = nan(1,numel(Source_dics_ERDS));
    v2find_ERD(:,:) = nan(1,numel(Source_dics_ERDS));  % if has  ERD

    for i = 1:numel(Source_dics_ERDS)
        v2find_ERS(:,i) = Source_dics_ERDS(i).maxpseudoT;
        v2find_ERD(:,i) = Source_dics_ERDS(i).minpseudoT;  % if has  ERD
    end

    [~,peakidx(1)] = min(v2find_ERD(ERDtime2sel(1):ERDtime2sel(2)));  % if has  ERD
    [~,peakidx(2)] = max(v2find_ERS(ERStime2sel(1):ERStime2sel(2))); 
    peakidx(2) = peakidx(2) + ERStime2sel(1)-1;

    Time2max(t,:) = [Source_dics_ERDS(peakidx).time(1)];
end















