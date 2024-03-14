function summary_plot_gain_tuning(results, plots)
%% plot the distribution of BF

%% let's plot the distribution of the best Frequence sampled
figure
hold on
h1 = plot(plots.wt.BestFs_counts/sum(plots.wt.BestFs_counts), '-ko', 'LineWidth', 1);
h2 = plot(plots.ko.BestFs_counts/sum(plots.ko.BestFs_counts), '-ro', 'LineWidth', 1);
set(h1, 'MarkerFaceColor', 'w')
set(h2, 'MarkerFaceColor', 'w')


legend([h1, h2], 'WT', 'KO')
xticks([1:2:9])
% xticklabels({'4k', '5.6k', '8k', '11.3k', '16k', '22.6k', '32k', '45.3k', '64k'})
xticklabels({'4k',  '8k',  '16k',  '32k', '64k'})

xlabel('Best Frequency')
ylabel('Proportion of Neurons')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
hold off

box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,400])

%% plot the Monotonic index
% figure
% [fig1, fig2] = ecdf_bar_plot([results.wt.MI], [results.ko.MI], 'Monotonic Index');
% set(get(fig1,'XLabel'), 'String', 'Monotonic Index');
% set(get(fig2,'YLabel'), 'String', 'Monotonic Index');
% ylim(fig2,[0,1])
figure
[fig1, fig2] = ecdf_bar_plot([plots.wt.MI], [plots.ko.MI], 'Monotonic Index');
set(get(fig1,'XLabel'), 'String', 'Monotonic Index');
set(get(fig2,'YLabel'), 'String', 'Monotonic Index');
ylim(fig2,[0,1])

%% plot the gains
figure
[fig1, fig2] = ecdf_bar_plot([results.wt.gains_value], [results.ko.gains_value], 'Gains');
set(get(fig1,'XLabel'), 'String', 'Gains (\Delta sp/s per 10 dB step)');
set(get(fig2,'YLabel'), 'String', 'Gains (\Delta sp/s per 10 dB step)');
xlim(fig1, [0,50])
ylim(fig2, [0,15])


%% plot the monotonic index with the gains
figure
h1 = scatter([results.wt.MI], [results.wt.gains_value], 'ko', 'filled');
h1.MarkerFaceAlpha = 0.3;
hold on
h2 = scatter([results.ko.MI], [results.ko.gains_value], 'ro', 'filled');
h2.MarkerFaceAlpha = 0.3
xlabel('Monotonic Index')
ylabel('Gains (\Delta sp/s per 10 dB step)')
box off
legend([h1, h2], {'WT', 'KO'})
set(gcf, 'Color', 'w')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
hold off
%% plot the gains with the BF
wt_gain = [results.wt.gains_value];
ko_gain = [results.ko.gains_value];
wt_BF = [results.wt.BestF];
ko_BF = [results.ko.BestF];

wt_BF(isnan(wt_gain)) =[];
ko_BF(isnan(ko_gain)) =[];

wt_gain(isnan(wt_gain)) =[];
ko_gain(isnan(ko_gain)) =[];

freqs = unique(wt_BF);

for i = 1:length(freqs)
    wt_indx = find(wt_BF == freqs(i));
    gain_wt{i} = wt_gain(wt_indx);
    wt_gain_avg(i) = mean(gain_wt{i});
    wt_gain_sem(i) = std(gain_wt{i})/sqrt(length(gain_wt{i}));
    
    ko_indx = find(ko_BF == freqs(i));
    gain_ko{i} = ko_gain(ko_indx);
    ko_gain_avg(i) = mean(gain_ko{i});
    if ko_gain_avg(i) == 0
        ko_gain_sem(i) = NaN;
    else
    ko_gain_sem(i) = std(gain_ko{i})/sqrt(length(gain_ko{i}))    ;
    
    end
end

figure;
h1 = errorbar(wt_gain_avg, wt_gain_sem, '-ko', 'LineWidth', 1);
hold on
h2 = errorbar(ko_gain_avg, ko_gain_sem, '-ro', 'LineWidth', 1);
set(h1, 'MarkerFaceColor', 'w')
set(h2, 'MarkerFaceColor', 'w')


legend([h1, h2], 'WT', 'KO')
xticks([1:2:9])
% xticklabels({'4k', '5.6k', '8k', '11.3k', '16k', '22.6k', '32k', '45.3k', '64k'})
xticklabels({'4k',  '8k',  '16k',  '32k', '64k'})

xlabel('Best Frequency')
ylabel('Gains (\Delta sp/s per 10 dB step)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
hold off

box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,400])

%% Plot the RLF (firing rates are normalized to the spontaneous FR)
figure;
line_errorbar_drc(plots.wt.rlf, plots.ko.rlf)
xticks([1:8])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlim([0,9])
xlabel('Sound Intensity(dB)')
ylabel('Firing Rate (Hz)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')

%% plot the Threshold
figure;
[fig1, fig2] = ecdf_bar_plot(plots.wt.threshold', plots.ko.threshold', 'Threshold');
set(get(fig1,'XLabel'), 'String', 'Threshold (dB) at BF');
set(get(fig2,'YLabel'), 'String', 'Threshold (dB) at BF');
ylim(fig2, [0,40])
set(gcf, 'Color', 'w')

%% plot the change of firing in a range
level = 0:10:70;
range_indx = find(level>=30 & level <=50);
genotype = {'wt', 'ko'};

for i = 1:length(genotype)
    temp_rlf = plots.(genotype{i}).rlf;
    plots.(genotype{i}).range_rlf = temp_rlf(:, range_indx);
    plots.(genotype{i}).range_gain      = plots.(genotype{i}).range_rlf(:, end) - plots.(genotype{i}).range_rlf(:, 1);
end

figure
[fig1, fig2] = ecdf_bar_plot(plots.wt.range_gain, plots.ko.range_gain, 'ranged gain');
set(get(fig1,'XLabel'), 'String', 'Gain (\Delta sp/s from 30 to 50 dB)');
set(get(fig2,'YLabel'), 'String', 'Gain (\Delta sp/s from 30 to 50 dB)');

%% plot the gain vs BF

% figure; 
% scatter(log2([results.wt.BestF]), [results.wt.gains_value])
% hold on
% scatter(log2([results.ko.BestF]), [results.ko.gains_value])

%% plot the d-prime
figure;
[fig1, fig2] = ecdf_bar_plot([results.wt.dprime], [results.ko.dprime], 'D-prime');
set(get(fig1,'XLabel'), 'String', 'D-prime');
set(get(fig2,'YLabel'), 'String', 'D-prime');
ylim(fig2, [0,8])
set(gcf, 'Color', 'w')

%% plot the summary tuning curve
figure
wt = plots.wt.tuning_avg_center;
ko = plots.ko.tuning_avg_center;
% wt = plots.wt.tuning_avg_center_n;
% ko = plots.ko.tuning_avg_center_n;
tuning_summary_plot(wt, ko)
title('Frequency tuning All neurons')

%only plot good_tuning
figure
wt_good = find([plots.wt.dprime]>4);
ko_good = find([plots.ko.dprime]>4);
tuning_summary_plot(wt(wt_good,:), ko(ko_good,:))
title('Frequency tuning D-prime > 4')
% set(gcf, 'Color', 'w')
% % export_fig('Frequency_tuning_good_n',  '-png', '-pdf')
% export_fig('Frequency_tuning_good',  '-png', '-pdf')

