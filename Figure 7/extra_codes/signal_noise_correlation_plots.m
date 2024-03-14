function signal_noise_correlation_plots(groupData, summaryData, freqs)
% this function is used to plot the data analyzed for signal-noise
% correlation
% load('summary_data_all.mat') to load the groupData, summaryData, freqs


% Plot 1: plot the distribution of BF at each session
figure
line_errorbar_drc(groupData.wt.bestFreq_ratio_session, groupData.ko.bestFreq_ratio_session)
% re-set the x-tick labels
x_ticks_indx = 1:1:length(freqs);
xticks(x_ticks_indx)
labels = {};
for i = 1:length(x_ticks_indx)
    labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
end
xticklabels(labels)
xlabel('Best Frequency (kHz)')
ylabel('Proportion')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])


% plot 2: plot the distribution of BF for all data
figure;
h1 = plot(groupData.wt.bestFreq_count/length(groupData.wt.bestFreq), '-ok', 'LineWidth',1, 'MarkerFaceColor', 'w');
hold on
h2 = plot(groupData.ko.bestFreq_count/length(groupData.ko.bestFreq), '-or', 'LineWidth',1, 'MarkerFaceColor', 'w');
hold on
x_ticks_indx = 1:1:length(freqs);
xticks(x_ticks_indx)
labels = {};
for i = 1:length(x_ticks_indx)
    labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
end
xticklabels(labels)
xlabel('Best Frequency (kHz)')
ylabel('Proportion')
legend([h1, h2], {'WT', 'KO'})
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])


% plot 3 and 4: plot the polulation signal correlation and noise correlation
%% assemble all data together
figure
title_label = 'Group Data';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.sig_cor70], [summaryData.ko.sig_cor70], title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.noise_cor70all], [summaryData.ko.noise_cor70all], title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');

% plot 5, 6, 7, 8: plot the population signal correlation and noise
% correlation based on the sign of the signal correlation (e.g., sig_cor>0 or sig_cor <0)
figure
wt_positive = groupData.wt.postive_indx;
ko_positive = groupData.ko.postive_indx;
title_label = 'Signal Correlation > 0';
[fig1, fig2] = ecdf_bar_plot(groupData.wt.all_sig_cor(wt_positive), groupData.ko.all_sig_cor(ko_positive), title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure
[fig1, fig2] = ecdf_bar_plot(groupData.wt.all_noise_cor(wt_positive), groupData.ko.all_noise_cor(ko_positive), title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');


figure
wt_negative = groupData.wt.negative_indx;
ko_negative = groupData.ko.negative_indx;
title_label = 'Signal Correlation < 0';
[fig1, fig2] = ecdf_bar_plot(groupData.wt.all_sig_cor(wt_negative), groupData.ko.all_sig_cor(ko_negative), title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure;
title_label = 'Signal Correlation < 0';
[fig1, fig2] = ecdf_bar_plot(groupData.wt.all_noise_cor(wt_negative), groupData.ko.all_noise_cor(ko_negative), title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');

% plot 9, 10 : plot the population signal correlation and noise
% correlation based on the sign of the signal correlation (e.g., sig_cor>0 or sig_cor <0)
figure;
title_label = 'Same BF';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.sameBest_sig_cor], [summaryData.ko.sameBest_sig_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure;
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.sameBest_noise_cor], [summaryData.ko.sameBest_noise_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');

% plot 11, 12 the correlation at different frequencies
%
figure
h_wt = line_sem_plot(1:length(freqs), groupData.wt.noise_cor, 'k');
hold on
h_ko = line_sem_plot(1:length(freqs), groupData.ko.noise_cor, 'r');
xticklabels(labels)
xlabel('Frequency (kHz)') 
ylabel('Noise Correlation Coefficient')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])
legend([h_wt, h_ko], {'WT', 'KO'})

figure
h_wt = line_sem_plot(1:9, groupData.wt.noise_cor(wt_positive,:),'k' )
hold on
h_wt_n = line_sem_plot(1:9, groupData.wt.noise_cor(wt_negative,:),'--k' )

h_ko = line_sem_plot(1:9, groupData.ko.noise_cor(ko_positive,:),'r' )
hold on
h_ko_n = line_sem_plot(1:9, groupData.ko.noise_cor(ko_negative,:),'--r' )
xticks(x_ticks_indx)
xticklabels(labels)
xlabel('Frequency (kHz)')
ylabel('Noise Correlation Coefficient')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])
legend([h_wt, h_wt_n, h_ko, h_ko_n], {'WT: positive', 'WT: negative', 'KO: positive', 'KO: negative'})
% 

%plot 13, 14, 15, 16 to see the correlation on the same shanks or different
%shanks
figure;
title_label = 'Same Shanks';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.sameShank_sig_cor], [summaryData.ko.sameShank_sig_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure;
title_label = 'Same Shanks';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.sameShank_noise_cor], [summaryData.ko.sameShank_noise_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');


figure;
title_label = 'Diff Shanks';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.diffShank_sig_cor], [summaryData.ko.diffShank_sig_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

figure;
title_label = 'Diff Shanks';
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.diffShank_noise_cor], [summaryData.ko.diffShank_noise_cor], title_label);
set(get(fig1,'XLabel'), 'String', 'Noise Correlation Coefficient');
set(get(fig2,'YLabel'), 'String', 'Noise Correlation Coefficient');
