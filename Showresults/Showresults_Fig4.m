%%
% Finger extension movement
% Fig.4
% creator: WxyZ - Yu Zheng
% Date: 20250314

%%
ft_defaults;

SubIdx = 1:16;
cm = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0.466 0.674 0.188];

load('.\mat\Classification\MeanAcc_allsub_Lowfreqband_4features.mat')
MeanAcc_LFB = MeanAcc;
load('.\mat\Classification\MeanAcc_allsub_Alpha_Beta.mat')
MeanAcc_AnB = MeanAcc;
load('.\mat\Classification\MeanAcc_allsub_Gamma.mat')
MeanAcc_G = MeanAcc;

clear MeanAcc

MeanAcc_Allband = [MeanAcc_LFB MeanAcc_AnB MeanAcc_G];

%% plot Acc   
subname = arrayfun(@(x) strcat('S', num2str(x, '%02d')), 1:numel(SubIdx), 'UniformOutput', false); subname = cat(2, subname, {'AVG'});
Freqname = {'LowFreq','Alpha','Beta','Gamma'};
acc = MeanAcc_Allband(SubIdx,:)*100; accstd = cat(1, zeros(size(acc)), std(acc, 0, 1)); 
acc2show = cat(1, acc, mean(acc, 1));

% ACC
fig = figure('Units','centimeters','Position',[10 10 14 6], 'IntegerHandle', 'off');
nexttile;
[~, hlgd] = wxyz_barplot(acc2show, accstd, flip([255 255 255; 255 210 217; 253 164 169; 242 119 115]/255), 'legend', Freqname, 'xtick', subname,...
    'ylabel', 'ACC (%)', 'show0err', false,'errBarWidth',0.6,'BarWidth',1);   % ,'errBarCapSize',3,'edgecolor',[0 0 0]
ylim([0 100]); set(gca, 'YTick', 0:20:100);
hold on; plot([size(acc2show, 1)-0.5 size(acc2show, 1)-0.5], get(gca, 'YLim'), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.6);
plot([0 size(acc2show, 1)+0.5], [28.5 28.5], '--', 'Color', 'r', 'LineWidth', 0.6)
set(hlgd, 'NumColumns', numel(Freqname), 'Location', 'northoutside', 'String', hlgd.String);
wxyz_decorFigure(fig, 'FontSize', 10, 'FontName', 'Calibri','FontScale',1.0, 'BoxWidth', 0.8,'Box','on','FontWeight','bold');
set(gca,'FontName', 'Calibri','FontSize',10); 
box on


%% add features
%% Errorbar
load('.\mat\Classification\MeanAcc_allsub_Lowfreqband_addfeatures.mat')
load('.\mat\Classification\MeanAcc_allsub_Lowfreqband_ACTseg.mat')
MeanAcc_ACT = MeanAcc;
Acc2plot = importdata(fullfile('.\mat\Classification\', 'Acc_addfeatures.mat'));
stderr = std(Acc2plot,0,1);

x = 1:16;

% ACC
h = figure('Units','centimeters','Position',[10 10 14 5.5], 'IntegerHandle', 'off');
tiledlayout(1,1, 'TileSpacing','tight','Padding','loose');  % ,'TileIndexing','columnmajor'
nexttile
bar(x,mean(Acc2plot),0.6,'FaceColor',wxyz_color_gradient(1, 1),'EdgeColor','k','LineWidth',0.4)
hold on

errorbar(x,mean(Acc2plot),stderr,'k','LineStyle','-',...
    'Marker','^','MarkerSize',6,'LineWidth',1.0,'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k'); 
ylim([50 90])
xlim([0.4 16.6])
ylabel('ACC (%)') 
set(gca,'XTick', 1:16, 'XTickLabel',{'Active','t1','t2','t3','t4','t1+t2','t1+t3','t1+t4', ...
    't2+t3', 't2+t4', 't3+t4','t1+t2+t3', ...
    't1+t2+t4', 't1+t3+t4', 't2+t3+t4','All'}) 

set(gca,'FontSize', 10, 'FontName', 'Calibri', 'FontWeight', 'bold') 
box on
set(findobj(h, 'Type', 'Axes'), 'LineWidth', 0.8); % axes box line



%% 5 time points comb
Acc2plot = importdata(fullfile('.\mat\Classification\', 'Acc_addfeatures.mat'));
Acc2plot_5cmb = [Acc2plot(:, 1) mean(Acc2plot(:, 2:5), 2) mean(Acc2plot(:, 6:11), 2) mean(Acc2plot(:, 12:15), 2) Acc2plot(:, 16)];
stderr = std(Acc2plot_5cmb, 0, 1);

x = 1:5;

% ACC
h = figure('Units','centimeters','Position',[10 10 4 5.5], 'IntegerHandle', 'off');
tiledlayout(1,1, 'TileSpacing','tight','Padding','loose'); 
nexttile
bar(x,mean(Acc2plot_5cmb),0.6,'FaceColor',wxyz_color_gradient(1, 1),'EdgeColor','k','LineWidth',0.4)
hold on

errorbar(x,mean(Acc2plot_5cmb),stderr,'k','LineStyle','-',...
    'Marker','^','MarkerSize',6,'LineWidth',1.0,'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k'); 
ylim([50 90])
xlim([0.4 5.6])
ylabel('ACC (%)')  
set(gca,'XTick', 1:16, 'XTickLabel',{'Active','1 time points','1 time points','3 time points','All'}) 

set(gca,'FontSize', 10, 'FontName', 'Calibri', 'FontWeight', 'bold') 
box on
set(findobj(h, 'Type', 'Axes'), 'LineWidth', 0.8); 


%%  confusion matrix
load('.\mat\Classification\confusionMatrix.mat')
confusionMatrix = confusionMatrix/16;

xname = {'Thumb','Index','Middle','Little'};
yname = {'Thumb','Index','Middle','Little'};

h = figure('Units','centimeters','Position',[10 10 9 7], 'IntegerHandle', 'off');
tiledlayout(1,1, 'TileSpacing','tight','Padding','compact');  

nexttile
f1 = heatmap(xname,yname,confusionMatrix);
cmap = [255 255 255; 255 210 217; 253 164 169; 242 70 70]/255;
ci = linspace(1, 256, 4);
cq = 1:256;
cmap = [interp1(ci, cmap(:,1), cq, 'linear')',...
            interp1(ci, cmap(:,2), cq, 'linear')',...
            interp1(ci, cmap(:,3), cq, 'linear')'];
colormap(cmap);
f1.Title = 'Confusion matrix';
f1.YLabel = 'True';
f1.XLabel = 'Predict';
set(gca,'FontSize', 10, 'FontName', 'Calibri') 


