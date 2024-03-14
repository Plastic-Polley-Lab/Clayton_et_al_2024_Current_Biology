function line_errorbar_drc(data1, data2)
% h1 = plot(mean(summary.wt.CrossCoef,1), '-ko', 'MarkerFaceColor', 'k', 'LineWidth', 1)
[r, c]=find(isnan(data1)); % remove NaN
data1(r,:) = [];
h1 = errorbar(mean(data1,1), std(data1, 0, 1)./sqrt(size(data1, 1)), '-ko');
set(h1, 'LineWidth', 1)
set(h1, 'MarkerFaceColor', 'w')
hold on
[r, c]=find(isnan(data2)); % remove NaN
data2(r,:) = [];
h2 = errorbar(mean(data2,1), std(data2, 0, 1)./sqrt(size(data2, 1)), '-ro');
set(h2, 'LineWidth', 1)
set(h2, 'MarkerFaceColor', 'w')
legend([h1, h2], {'WT', 'KO'})
ylabel('CorrCoef')
xlim([0,8])
xticks([0:7])
xticklabels({'','20', '15', '10', '5', '0', '-5', '-10'})
xlabel('Signal to Noise Ratio (dB)')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,300])