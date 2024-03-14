function tuning_summary_plot(wt, ko)
wt_m = mean(wt,1, 'omitnan');
wt_std = std(wt,0, 1, 'omitnan');
for i = 1:length(wt_m)
    wt_count(i) = length(find(~isnan(wt(:,i))));
end
wt_sem = wt_std./sqrt(wt_count);

ko_m = mean(ko,1, 'omitnan');
ko_std = std(ko,0, 1, 'omitnan');
for i = 1:length(ko_m)
    ko_count(i) = length(find(~isnan(ko(:,i))));
end
ko_sem = ko_std./sqrt(ko_count);


h1 = errorbar(wt_m, wt_sem, '-ko');
% h1 = errorbar(plot_data_m_n, plot_data_sem_n, '-ko');
set(h1, 'LineWidth', 1)
set(h1, 'MarkerFaceColor', 'w')
hold on
h2 = errorbar(ko_m, ko_sem, '-ro');
% % h2 = errorbar(plot_data_m_ko_n, plot_data_sem_ko_n, '-ro');
set(h2, 'LineWidth', 1)
set(h2, 'MarkerFaceColor', 'w')
legend([h1, h2], {'WT', 'KO'})
ylabel('Firing Rate (Hz)')
xticks([1:17])
xticklabels({'-4', '-3.5', '-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0', ...
           '0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4'})
xlabel('Relative to BF (oct)')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,400])
xlim([4, 14])