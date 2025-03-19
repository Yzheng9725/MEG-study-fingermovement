function [trl] = wxyz_emgdetect_v5(opt, emg, fs, trlnum, channel)
% wxyz_emgdetect are used to detect trigger based on EMG measurement.
% Notice!
% 1. this function can only to process 1 channel per times.
% 2. Type of trials should be equal.
% Some code was referenced from Fieldtrip.
% see also https://www.fieldtriptoolbox.org/example/detect_the_muscle_activity_in_an_emg_channel_and_use_that_as_trial_definition/

% Threshold: μ+σ×k    μ & σ specifies mean & std; k(user set, row 32) set to 1.5 
% BandPass [20 450] - rectification - Lowpass 50 - Threshold

% set the threshold defaults
% mintrlgap   = ft_getopt(opt, 'mintrlgap', 3.0); % second
maxreact    = ft_getopt(opt, 'maxreact', 1.5); % second
minact      = ft_getopt(opt, 'minact', 0.2); % second
resflag     = ft_getopt(opt, 'resflag', false);
triallength = ft_getopt(opt, 'triallength', 3); % second

trl = zeros; 
figure;
tiledlayout(6,ceil(trlnum/6),'TileSpacing','compact','Padding','compact');

for n = 1:trlnum
    emg2det = emg.trial{1,n}(channel,:);

    % process
    emgflt  = ft_preproc_bandpassfilter(emg2det, fs, [20 450], 6);      % bandpassfilter
    emgrec = abs(emgflt);                                               % rect
    emgsmt  = ft_preproc_lowpassfilter(emgrec, fs, 50, 6);

    % Threshold
    threshold = mean(emgsmt) + 1.5*std(emgsmt);

    emgtrl = zeros(size(emgsmt));
    emgtrl(find(emgsmt > threshold)) = 1;
    emgtrldif  = diff(emgtrl, [], 2);

    % get raw onset and offset trig
    emgonraw   = find(emgtrldif(:) == 1);
    emgoffraw  = find(emgtrldif(:) ==-1);

    emgon = emgonraw; % used to update trigger, rather than emgonraw.
    emgoff = emgoffraw;

    [emgpow, num] = max(emgsmt);
 
    for i = 1:length(emgon)
        if emgon(i) < num-fs*maxreact
            emgtrldif(emgon(i)) = 0;
        else
            c1 = emgon(i) > num;
            c2 = emgon(i) < fs*0.2;  
            if c1 || c2 
                emgtrldif(emgon(i)) = 0;
                emgtrldif(emgoff(i)) = 0;
                emgon(i) = nan;
            end
        end
    end

    emgon(isnan(emgon)) = [];
    if isempty(emgon)
        trl(n,1) = 0;
    else
        trl(n,1) = emgon(1);
    end

    % show results
    if resflag
        % show processed emg signal and trigger.
        nexttile;
        hold on;
        h1 = plot(emg.time{n}(1,:), emgsmt, 'Color', 'k', 'LineWidth', 1);
        h2 = plot(emg.time{n}(2:end), (emgtrldif(:)==1)*20, 'Color', 'r', 'LineWidth', 1.2);
        plot([trl(n,1)/fs trl(n,1)/fs], [0 200*ones(1, length(trl(n,1)))] , '--', 'Color', 'r','LineWidth', 1);
        plot([emg.time{n}(1) emg.time{n}(end)], [threshold threshold], 'Color', 'g', 'LineWidth', 1.2);
        xlim([emg.time{n}(1) emg.time{n}(end)-1]);
        ylim([-2 250])
        xlabel('Time'); ylabel('Amplitude'); title('Processed emg signal and trigger');

        % plot EMG2check
        hold on
        yyaxis right
        plot(emg.time{n}(1,:),emg2det,'Color', 'b', 'LineWidth', 0.5)
        ylim([-1000 1000])
        set(gca, 'ytick',[])

        if n == trlnum
            legend([h1 h2], {'Processed emg signal', 'Trigger'});
        end
        set(gca, 'FontName', 'Times New Roman');
    end
end