%%  
% Fig.2(d),(e)

%%
ft_defaults;

FN = {'Alpha','Beta','Gamma'};
FB = {[8 13],[13 30],[60 90]};

%%
for r = 1:numel(FN)
    opt.freqname = FN{r};
    opt.freqband = FB{r};

    load(['.\mat\SensorLevel\ERDS_stat_alltaskavg_' opt.freqname '.mat'])
    load(['.\mat\SensorLevel\ERDS_alltaskavg_' opt.freqname '.mat'])
    load(['.\mat\SensorLevel\ERDS_alltask_' opt.freqname '.mat'])

    ERDS_alltaskavg.avg = squeeze(mean(ERDS_alltask(:,:,:),1));
    eval(['ERDS_alltaskavg_' opt.freqname ' = ERDS_alltaskavg;']);


    % cal signaf time
    chanidx = find(strcmp(ERDS_alltaskavg.label, 'MLF67'));
    ersclusterv = ERDS_stat.posclusterslabelmat(chanidx,:);
    erstimeidx = find(ersclusterv == 1);
    erstime(r,:) =  [ERDS_stat.time(erstimeidx(1)) ERDS_stat.time(erstimeidx(end))];

    if r~=3
        erdclusterv = ERDS_stat.negclusterslabelmat(chanidx,:);
        erdtimeidx = find(erdclusterv == 1);
        erdtime(r,:) =  [ERDS_stat.time(erdtimeidx(1)) ERDS_stat.time(erdtimeidx(end))];
    end
end


%% plot
cm = [43 174 133; 18 107 174; 204 85 149]/255;
h = figure('Units','centimeters','Position',[10 10 9 4], 'IntegerHandle', 'off');
tiledlayout(1, 1,  'TileSpacing','tight','Padding','tight');  
nexttile

plot([0 0],[-150 150],'Color',[0.5 0.5 0.5],'LineWidth',0.8,'LineStyle','--')
hold on
plot([-2.5 2.5],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',0.8,'LineStyle','--')

% alpha & beta
for r = 1:2
    opt.freqname = FN{r};  %'Alpha','Beta',
    opt.freqband = FB{r}; 
    eval(['ERDS2plot = ERDS_alltaskavg_' opt.freqname ';']);

    timeidx(1) = find(ERDS2plot.time == erdtime(r,1));
    timeidx(2) = find(ERDS2plot.time == erdtime(r,2));
    timeidx(3) = find(ERDS2plot.time == erstime(r,1));

    t1 = ERDS2plot.time(1:timeidx(1));
    t2 = ERDS2plot.time(timeidx(1):timeidx(2));
    t3 = ERDS2plot.time(timeidx(2):timeidx(3));
    t4 = ERDS2plot.time(timeidx(3):end);

    %% plot
    plot(t1, ERDS2plot.avg(chanidx,1:timeidx(1)),'Color',[0.85 0.325 0.098],'LineWidth',1.5,'LineStyle',':');
    lin(r)= plot(t2, ERDS2plot.avg(chanidx,timeidx(1):timeidx(2)),'Color',cm(r,:),'LineWidth',1.5,'LineStyle','-');
    plot(t3, ERDS2plot.avg(chanidx,timeidx(2):timeidx(3)),'Color',cm(r,:),'LineWidth',1.5,'LineStyle',':');
    plot(t4, ERDS2plot.avg(chanidx,timeidx(3):end),'Color',cm(r,:),'LineWidth',1.5,'LineStyle','-');

end

% gamma
r = 3;
opt.freqname = FN{r}; 
opt.freqband = FB{r};
eval(['ERDS2plot = ERDS_alltaskavg_' opt.freqname ';']);

timeidx(1) = find(ERDS2plot.time == erstime(r,1));
timeidx(2) = find(ERDS2plot.time == erstime(r,2));

t1 = ERDS2plot.time(1:timeidx(1));
t2 = ERDS2plot.time(timeidx(1):timeidx(2));
t3 = ERDS2plot.time(timeidx(2):end);

plot(t1, ERDS2plot.avg(chanidx,1:timeidx(1)),'Color',cm(r,:),'LineWidth',1.5,'LineStyle',':');
lin(r) = plot(t2, ERDS2plot.avg(chanidx,timeidx(1):timeidx(2)),'Color',cm(r,:),'LineWidth',1.5,'LineStyle','-');
plot(t3, ERDS2plot.avg(chanidx,timeidx(2):end),'Color',cm(r,:),'LineWidth',1.5,'LineStyle',':');

plot([-0.2 -0.2], [-25 55], 'LineStyle','--', 'Color',[1 0 0 0.6],'LineWidth',0.8)
plot([1.0 1.0], [-25 55], 'LineStyle','--', 'Color',[1 0 0 0.6],'LineWidth',0.8)
plot([1.6 1.6], [-25 55], 'LineStyle','--', 'Color',[1 0 0 0.6],'LineWidth',0.8)
plot([0 0], [-25 55], 'LineStyle','--', 'Color',[1 0 0 0.6],'LineWidth',0.8)

xlim([-1.5 2])
ylim([-25 55])
set(gca,'XTick',[-1.5 0 2.0],'FontSize', 10)
set(gca,'YTick',[-25 0 55],'FontSize', 10)
box on
ylabel('ERD/ERS (%)','FontSize', 10,'VerticalAlignment','bottom');
xlabel('Time (s)','FontSize', 10);
set(gca,'FontName', 'Calibri', 'FontWeight', 'bold');

legend([lin(1) lin(2) lin(3)],'Alpha','Beta','Gamma','FontSize',10,'Box','off','Location','eastoutside'); 


%% cluster topo alpha/beta
baseline = [-1.5 -1];

for r = 1:3  % 1:3
    opt.freqname = FN{r};
    opt.freqband = FB{r};
    eval(['ERDS2plot = ERDS_alltaskavg_' opt.freqname ';']);
    load(['.\mat\SensorLevel\ERDS_stat_alltaskavg_' opt.freqname '.mat'])

    h = figure('Units','centimeters','Position',[10 10 2.5 4], 'IntegerHandle', 'off');
    drawnow;
    tiledlayout(2,1,'TileSpacing','tight','Padding','compact');
    cfg                 = [];
    cfg.baseline        = baseline;
    cfg.layout          = 'CTF275_helmet.mat';
    cfg.colorbar        = 'no';
    cfg.masknans        = 'no';
    cfg.showlabels      = 'no';
    cfg.zlim            = [-50 50];
    cfg.comment         = 'no'; % remove bottom left corner text
    cfg.marker          = 'off';
    cfg.highlight          = 'on';
    cfg.highlightsymbol    = '.';
    cfg.highlightcolor     = [0 0 0];
    cfg.highlightsize      = 5;
    cfg.style              = 'straight';
    cfg.gridscale          = 150;
    cfg.showoutline     = 'yes'; % show head/helmet outline
    cfg.showscale       = 'no'; % remove bottom right corner graph

    % ERS
    if r == 1
        cfg.xlim            = [1.6 1.6];
        cfg.highlightchannel   =  ERDS2plot.label(find(ERDS_stat.posclusterslabelmat(:,1081) == 1));  % 1081 901 601
    elseif r == 2
        cfg.xlim            = [1.0 1.0];
        cfg.highlightchannel   =  ERDS2plot.label(find(ERDS_stat.posclusterslabelmat(:,901) == 1));  % 1081 901 601
    elseif r == 3
        cfg.xlim            = [0 0];
        cfg.highlightchannel   =  ERDS2plot.label(find(ERDS_stat.posclusterslabelmat(:,601) == 1));  % 1081 901 601
    end


    nexttile
    cfg.figure          = gcf;
    ft_topoplotER(cfg, ERDS2plot);
    set(gca,'PlotBoxAspectRatio',[1 1 1])

    % ERD
    if r == 1||r ==2
        nexttile
        cfg.xlim            = [-0.2 -0.2];
        cfg.highlightchannel   =  ERDS2plot.label(find(ERDS_stat.negclusterslabelmat(:,540) == 1));
        cfg.figure          = gcf;
        ft_topoplotER(cfg, ERDS2plot);
        set(gca,'PlotBoxAspectRatio',[1 1 1])
    end

    colormap(wxyz_colormap(1))

    if r == 3
        c=colorbar;
        set(c,'Position',[0.54,0.12,0.08,0.35],'Ticks',[-50 0 50],'TickLabels',{'-50','0', '50'},'FontSize', 10)
        c.Label.String = 'ERD/ERS (%)';
        c.Label.VerticalAlignment = 'top';
    end

    set(gca,'FontName', 'Calibri','FontWeight', 'bold');  % Times New Roman

end






