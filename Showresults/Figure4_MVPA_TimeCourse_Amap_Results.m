%% This figure is shown in Figure 4.
%%
clc; clear; close all;
ft_defaults;

%%
opt = [];
opt.filepath = fileparts(pwd);
opt.sub = 'S00';
opt.paradigm = 'Extension_4tasks';
opt.motion = 'Extension';
opt.mripath = fullfile('E:', 'MRI', opt.sub, 'forMEG');
opt.task = {'Thumb','Index','Middle','Little'};
opt.ana = 'Source_Space';
opt.subana = 'MVPA-LeftSMC';
opt.groupmatdir = fullfile(opt.filepath, 'mat', opt.ana, opt.subana);
opt.groupresultdir = fullfile(opt.filepath, 'results', opt.ana, opt.subana);

if ~exist(opt.groupresultdir, 'dir')
    mkdir(opt.groupresultdir);
end

FN = {'Allband','Delta','Theta','Alpha','Beta','Gamma','High-gamma','LowFreqBand'};
FB = {[0.5 90],[0.5 4],[4 8],[8 13],[13 30],[30 60],[60 90],[0.5 8]};

%% Amap MVPA ACC
ACC_group = importdata(fullfile(opt.groupmatdir, 'showresults_ACC_MVPA_AImap_16sub.mat'));
ACC_stat = importdata(fullfile(opt.groupmatdir, 'showresults_stat_ACC_MVPA_AImap_16sub.mat'));

%% Plot Amap MVPA ACC with significant time points
t = ACC_stat(1).time;
colors = [65, 82, 84; 73, 159, 162; 162, 212, 210; 245, 215, 221; 245, 175, 191; 123, 105, 110]/255;
fig = figure('Units', 'centimeters', 'Position', [10 10 15 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
nexttile; hold on;
plot([0 0], [10 60],'LineStyle','--','Color','k','LineWidth',0.5);
plot([-1 2], [25 25],'LineStyle','--','Color','k','LineWidth',0.5);
xlim([-1 2]); ylim([20 50]);
hlist = gobjects(1,6);
for f = 2:7 % delta/theta/alpha/beta/gamma/high-gamma
    acc_avg = mean(ACC_group{f})*100;
    acc_std = std(ACC_group{f})*100;
    acc_sem = acc_std / sqrt(size(ACC_group{f},1));
    [hl, hp] = boundedline(t, acc_avg, acc_sem, 'alpha');
    set(hl, 'Color', colors(f-1,:), 'LineWidth', 1.2);
    set(hp, 'FaceColor', colors(f-1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'EdgeAlpha', 0.3);
    plot_time_regions(t, ACC_stat(f).mask, 'color', colors(f-1,:), 'alpha', 1, 'ylim', [20+(f-2)*0.5 20+(f-1)*0.5]);
    hlist(f-1) = hl;
end
xlabel('Time (s)','FontWeight','bold');
ylabel('Accuracy (%)','FontWeight','bold');
legend(hlist, {'Delta','Theta','Alpha','Beta','Gamma','High-gamma'}, 'Location', 'northeast');
set(gca, 'Color', 'none', 'XTick', -1:1:2);
set(gcf, 'Name', strcat('MVPA-ACC-AImap-AllFreqBands-sem'));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% Time course MVPA ACC
ACC_group = importdata(fullfile(opt.groupmatdir, 'showresults_ACC_MVPA_Timecourse_16sub.mat'));
ACC_stat = importdata(fullfile(opt.groupmatdir, 'showresults_stat_ACC_MVPA_Timecourse_16sub.mat'));

%% Plot all together
t = ACC_stat(1).time;
colors = [65, 82, 84; 73, 159, 162; 162, 212, 210; 245, 215, 221; 245, 175, 191; 123, 105, 110]/255;
fig = figure('Units', 'centimeters', 'Position', [10 10 15 5], 'IntegerHandle', 'off');
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'tight');
nexttile; hold on;
plot([0 0], [10 60],'LineStyle','--','Color','k','LineWidth',0.5);
plot([-1 2], [25 25],'LineStyle','--','Color','k','LineWidth',0.5);
xlim([-1 2]); ylim([20 50]);
hlist = gobjects(1,6);
for f = 2:7 % delta/theta/alpha/beta/gamma/high-gamma
    acc_avg = mean(ACC_group{f})*100;
    acc_std = std(ACC_group{f})*100;
    acc_sem = acc_std / sqrt(size(ACC_group{f},1));
    [hl, hp] = boundedline(t, acc_avg, acc_sem, 'alpha');
    set(hl, 'Color', colors(f-1,:), 'LineWidth', 1.2);
    set(hp, 'FaceColor', colors(f-1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'EdgeAlpha', 0.3);
    plot_time_regions(t, ACC_stat(f).mask, 'color', colors(f-1,:), 'alpha', 1, 'ylim', [20+(f-2)*0.5 20+(f-1)*0.5]);
    hlist(f-1) = hl;
end
xlabel('Time (s)','FontWeight','bold');
ylabel('Accuracy (%)','FontWeight','bold');
legend(hlist, {'Delta','Theta','Alpha','Beta','Gamma','High-gamma'}, 'Location', 'northeast');
set(gca, 'Color', 'none', 'XTick', -1:1:2);
set(gcf, 'Name', strcat('MVPA-ACC-Timecourse-AllFreqBands-sem'));
print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.png')), '-dpng', '-r600', '-image');
% print(fig, fullfile(opt.groupresultdir, strcat(fig.Name, '.eps')), '-depsc2', '-r600', '-vector');

%% SUBFUNCTION
function plot_time_regions(t, sig, varargin)

p = inputParser;
addParameter(p,'method','patch');
addParameter(p,'color',[0.8 0.8 0.8]);
addParameter(p,'alpha',0.3);
addParameter(p,'width',[]);
addParameter(p,'ylim',ylim);
parse(p,varargin{:});

method = p.Results.method;
color  = p.Results.color;
alpha  = p.Results.alpha;
width  = p.Results.width;
y1     = p.Results.ylim(1);
y2     = p.Results.ylim(2);

sig = logical(sig(:).');
t = t(:).';
if isempty(width)
    dt = mean(diff(t));
else
    dt = width;
end
d = diff([false sig false]);
idx_start = find(d == 1);
idx_end   = find(d == -1) - 1;
hold on;
for k = 1:numel(idx_start)
    x1 = t(idx_start(k));
    x2 = t(idx_end(k));
    x1 = x1 - dt/4;
    x2 = x2 + dt/4;
    
    switch lower(method)
        case 'patch'
            patch([x1 x2 x2 x1], [y1 y1 y2 y2], color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
        case 'rectangle'
            rectangle('Position',[x1 y1 (x2-x1) (y2-y1)], 'FaceColor', [color alpha], 'EdgeColor', 'none');
    end
end
hold off;
end
