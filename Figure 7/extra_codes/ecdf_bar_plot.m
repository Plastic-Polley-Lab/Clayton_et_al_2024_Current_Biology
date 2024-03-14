function [fig1, fig2] = ecdf_bar_plot(wt, ko, title_label)
% this function is used to plot the ecdf and bar of the wt and ko data

% plot the ecdf graph of the data
% figure(100);
sgtitle(title_label)
fig1 = subplot(1, 2, 1);
[wt_f, wt_x] = ecdf(wt);
[ko_f, ko_x] = ecdf(ko);
h_wt = plot(wt_x, wt_f, '-k', 'LineWidth',1);
hold on
h_ko = plot(ko_x, ko_f, '-r', 'LineWidth',1);
legend([h_wt, h_ko], {'WT', 'KO'})
xlabel('Change xlabel')
ylabel('Cumulative Probability')
box off
set(gcf, 'Color', 'w')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
hold off

% set(get(fig1,'XLabel'), 'String', 'Signal Correlation Coefficient');
% set(get(fig2,'YLabel'), 'String', 'Signal Correlation Coefficient');

% xlim([0,70])
% xlim([0,150])
% savefig('time_scale_exc_5msBin.fig')
% export_fig('time_scale_exc_5msBin',  '-png')
%

%plot the bar graph
% figure(101);
fig2 = subplot(1, 2, 2);
% CT=cbrewer('div', 'RdYlBu', 6);
wt_avg = mean(wt, 'omitnan');
wt_n = length(wt(~isnan(wt)));
wt_sem = std(wt, 'omitnan')/sqrt(wt_n);


ko_avg = mean(ko, 'omitnan');
ko_n = length(ko(~isnan(ko)));
ko_sem = std(ko, 'omitnan')/sqrt(ko_n);;


h1 = bar(1,wt_avg,0.7,'LineWidth',1);
h1.FaceColor = 'k';
h1.FaceAlpha = 0.3;

hold on
h2 = bar(2,ko_avg,0.7,'LineWidth',1);
h2.FaceColor = 'r';
h2.FaceAlpha = 0.7;

e = errorbar([wt_avg, ko_avg],[wt_sem, ko_sem],'k', 'LineWidth',1);
e.LineStyle = 'none';
e.CapSize = 15;


xticks([1,2])
xticklabels({'', 'WT', 'KO', ''})
ylabel('Change ylabel')
% ylabel('Time scale (ms)')
% ylabel('lags (ms) before 0 crossing')

set(gca,'XTick',0:3)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[500,200,800,400])
set(gcf, 'Color', 'w')
xlim([0,3])
hold off
% ylim([0.4,1])
% name = {'','No stimulation','Light stimulation',''}
% set(gca,'xticklabel',name)
wt_n = (wt - wt_avg)./std(wt, 'omitnan');
ko_n = (ko - ko_avg)./std(ko, 'omitnan');
h_ks1 = kstest(wt_n);
h_ks2 = kstest(ko_n);
if h_ks1 == 0 && h_ks2 ==0
    fprintf('Date are from normal distribution')
    [h, p,~,stats] = ttest2(wt,ko)
else
    fprintf('Data are not normal distributed')
    [h, p,stats] = ranksum(wt,ko)
    [h, p,stats] = kstest2(wt,ko)
end

