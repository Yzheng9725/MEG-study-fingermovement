%%
% Finger extension movement
% Fig.6
% creator: WxyZ
% Date: 20250314

%%
ft_defaults;

SubIdx = 1:16;
baseline = [-1.5 -1.0];

Timeseg = [];
Timestep = 0.02;
Timeseg = [-0.5:Timestep:1-Timestep] + 0.01;

%%
opt = [];
opt.mripath = '.\mat\S00\forMEG\';
opt.wbpath = '.\mat\S00\workbench\';

opt.task = {'Thumb','Index','Middle','Little'};

%%
load([opt.mripath 'sourcemodel.mat'])
load([opt.mripath 'brainatlas.mat'])
atlas_wb_L = brainatlasL;
atlas_wb_R = brainatlasR;

atlasLabel = [21 23];
Atlasannot_Lshpere = [atlas_wb_L.parcellation]; 

% select pre/post-all ROI pos index
ROIidx = [];
for i = 1:numel(atlasLabel)
    idx = find(Atlasannot_Lshpere == atlasLabel(i));
    ROIidx = [ROIidx;idx];  
    clearvars idx
end

ROIidx = sort(ROIidx,'ascend');  % 升序

load([opt.mripath 'sourcemodel.mat'])
load([opt.mripath 'sourcemodel_inflated.mat'])

%% plot Vertices number - ACT-ROI Fig.6(a)
load('.\mat\Digit-map\AllTask_Acc_TASKnREST_Group.mat')
load('.\mat\Digit-map\Accstat_TASKnREST_fdr.mat')

ACTROImask = Source_stat_nps.mask;

for c = 1:75
    vernum(c) = numel(find(Source_stat_nps.mask(c,:)));
end

h = figure('Units','centimeters','Position',[10 10 7 4], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','tight','Padding','tight');
drawnow;

nexttile([1 1])
time2plot  = -0.1:0.1:0.5;

for c = 1:numel(time2plot)
    rectangle('Position',[time2plot(c) 0.1 0.02 799.9],'FaceColor',[1 0 0 0.2],...
        'EdgeColor','none')
    hold on
end

plot(Timeseg, vernum,'Linestyle','-','LineWidth',1.2,'Color','k')
hold on
plot([0 0], [0 1000],'Color',[0.5 0.5 0.5],'LineStyle','--')

xlim([-0.5 1.0])
ylim([0 800])

set(gca,'XTick',[-0.5 0 0.5 1.0],'FontSize', 10)
set(gca,'YTick',[0 400 800],'FontSize', 10)
box on

ylabel('Number of vertices','VerticalAlignment','bottom');
xlabel('Time (s)','VerticalAlignment','cap');
set(gca,'FontName', 'Calibri', 'FontWeight', 'bold');  % Times New Roman



%%  ACT-ROI plot  - Fig.6(b)
idx2plot = [21:5:51];

h = figure('Units','centimeters','Position',[10 10 14 6.8], 'IntegerHandle', 'off');
tiledlayout(2, 4,  'TileSpacing','tight','Padding','tight');  % ,'TileIndexing','columnmajor'

cmap = colormap(gray(1000));
normsulc = -1* sourcemodel_inflated.sulc;
normsulc = (normsulc - min(normsulc)) / (max(normsulc-min(normsulc)));
num_colors = 850;
color_indices = round(normsulc * (num_colors - 1)) +1;
cmapped = cmap(150+color_indices,:);

for i = idx2plot
    acc2plot = nan(size(sourcemodel.thickness));
    acc2plot(ROIidx) = squeeze(mean(AccAlltask(SubIdx,i,ROIidx)))' .* Source_stat_nps.mask(i,:);
    acc2plot(acc2plot == 0) = NaN;
    acc2plot = acc2plot/100;

    %%
    nexttile;
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', cmapped);
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', acc2plot);  %lighting gouraud; material dull;light
    view(-110,40)
    clim([0.5 0.7])
    colormap('hot')

    title([num2str(Timeseg(i)) ' s'],'Position',[-10 -15 120])
    set(gca,'FontSize', 10, 'FontName', 'Calibri')
    set(gca, 'PlotBoxAspectRatio', [1 1 1]);
    drawnow
    hold off
end

cb1 = colorbar;
set(cb1,'Position', [0.87 0.10 0.03 0.3])
set(cb1,'Box','off','FontSize',10,'Ticks',[0.5 0.7])
cb1.Label.String = 'ACC';
cb1.Label.VerticalAlignment = 'baseline';

set(gca,'FontName', 'Calibri' ,'FontSize', 10, 'FontWeight', 'bold');  
sgtitle('Group Active-ROIs','FontName', 'Calibri' ,'FontSize', 12, 'FontWeight', 'bold');




%% DIGIT MAP - Fig.6(c)
load('.\mat\Digit-map\Accstat_TASKnREST_fdr.mat')
ACTROImask = Source_stat_nps.mask;

load('.\mat\Digit-map\AllTask_Acc_4TASKs_Group.mat')
load('.\mat\Digit-map\DigitMap_Acc_4TASKs_Group.mat')
load('.\mat\Digit-map\Accstat_4TASK_fdr.mat')

Time2map = Source_stat_nps.time;

%% cal Group - 4 tasks devote - ACT-ROI fdr
Map_alltask_group = [];

for ss = 1:16
    for t = 1:4
        if ss == 1
            Map_alltask_group{t} = zeros(numel(Time2map),numel(sourcemodel.thickness));
        end
        
        map_ind = [];
        map_ind = zeros(numel(Time2map), numel(sourcemodel.thickness));
        for c = 1:numel(Time2map)
            map_ind(c,ROIidx) = squeeze(Map_alltask_allsub(ss,c+20,ROIidx))' .* Source_stat_nps.mask(c,:);  % 去除未通过检验的vertex
        end

        map_ind(map_ind ~= t) = 0;
        map_ind(map_ind == t) = 1;

        Map_alltask_group{t} = Map_alltask_group{t} + map_ind;
 
    end

end

%  DEVOTE- ACTIVE ROI - 4 TASKS IN 1 BRAIN - 4TASKS - Group
for c = 1:numel(Time2map)
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

%% Group Digit-map - show
h = figure('Units','centimeters','Position',[10 10 14 16], 'IntegerHandle', 'off');
tiledlayout(4, 4,  'TileSpacing','tight','Padding','tight');  % ,'TileIndexing','columnmajor'
Cidx = 5:20;

cm = [nan nan nan; 59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;
% colormap for tasks

cmap = colormap(gray(1000));
% colormap for surface
normsulc = -1* sourcemodel_inflated.sulc;
normsulc = (normsulc - min(normsulc)) / (max(normsulc-min(normsulc)));
num_colors = 850;
color_indices = round(normsulc * (num_colors - 1)) +1;
cmapped = cmap(150+color_indices,:);

for c = Cidx
    map2plot = zeros(size(sourcemodel.thickness));
    map2plot = Map_4taskin1brain(c,:)+1; 

    nexttile

    %% plot sulc/gri  
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', cmapped);
    
    %% plot map
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', cm(map2plot,:)); 
    view(-110,40)

    title([num2str(Time2map(c)) ' s'],'Position',[-5 -15 120]);
    set(gca,'FontSize', 10, 'FontName', 'Calibri')
    set(gca, 'PlotBoxAspectRatio', [1 1 1]); 
    drawnow

end

colormap(cm(2:end,:))
sgtitle('Group Digit-MAPs', 'FontName', 'Calibri','FontSize', 12, 'FontWeight', 'Bold')




%% Group Digit-map - show - 0.05s 
cm = [nan nan nan; 59 135 199; 242 120 115; 255 211 115; 54 151 88]/255;
h = figure('Units','centimeters','Position',[10 10 6 4], 'IntegerHandle', 'off');
tiledlayout(1, 5,  'TileSpacing','compact','Padding','tight');  

for c = 8
    map2plot = zeros(size(sourcemodel.thickness));
    map2plot = Map_4taskin1brain(c,:)+1;  

    nexttile([1 4])
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', cmapped);
    hold on

    %% plot map
    ft_plot_mesh(sourcemodel_inflated,'vertexcolor', cm(map2plot,:));  % ,'facealpha',0.6,  lighting gouraud; material dull;light
    view(-110,40)

    set(gca, 'PlotBoxAspectRatio', [1 1 1]); % 所有图窗大小一致
    drawnow
end
hold on
 
colormap(cm(2:end,:))
cb1 = colorbar;
set(cb1,'Position', [0.7 0.35 0.05 0.4])
set(cb1,'YTick', 0.125:0.25:0.875)
set(cb1,'YTickLabel',{'Thumb' ,'Index','Middle','Little'})   % ,'ALL'
set(gca,'FontSize', 10, 'FontName', 'Calibri','FontWeight', 'Bold') % , 'FontWeight', 'bold'




%% cal vertices number for task  - Fig.6(d)
digitnum = [];
for t = 1:4
    for c = 1:35
        digitnum(t,c) = numel(find(Map_4taskin1brain(c,:) == t));
    end
end

digitnumsum = sum(digitnum, 2);
digitnum_thres = mean(digitnum, 1);

% joy plot
data2plot = digitnum';

lgLable= {'Thumb','Index','Middle','Little'};
colors=[0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

h = figure('Units','centimeters','Position',[10 10 7 4.5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing','tight','Padding','tight');
drawnow;
nexttile([1 1])

p = 120;
x = -0.09:0.02:0.59;
for i = size(data2plot, 2):-1:1
    f = data2plot(:,i)';
    fShifted = f + i * p;
    pHandle = plot(x, fShifted,'color',colors(i,:),'LineStyle','-','LineWidth', 1.2,'HandleVisibility', 'off');
    hold on;
    plot(x, ones(size(x))*(max(digitnum_thres)*0.25 + i*p),'LineStyle','--', 'color', [0.5 0.5 0.5]);
    yline(i * p, '-', 'LineWidth', 0.5, 'Color',[0.7 0.7 0.7],'HandleVisibility', 'off')
    Xfill = [x, fliplr(x)];
    Yfill = [fShifted, ones(1, length(x)) * i * p];
    fill(Xfill, Yfill, pHandle.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
yTick = (1:size(data2plot, 2)) * p;
box off
xlim([-0.1 0.6])
ylim([80 780])
set(gca, 'YTick', yTick, 'YTickLabel', lgLable, 'TickDir', 'out','YColor','k');
set(gca,'XTick',[-0.1 0 0.6])

title('Number of vertices','FontName', 'Calibri','FontSize', 12, 'FontWeight', 'bold');
xlabel('Time (s)','VerticalAlignment','Middle');
set(gca,'FontName', 'Calibri','FontSize', 10, 'FontWeight', 'bold');  % Times New Roman











