function [h1, h2] = bar_plot(a)
% This code is only used for paired data
% Input:
%      a: n x 2 matrix, where each row represents a sample, a may contains
%      NaN
% figure;
CT=cbrewer('div', 'RdYlBu', 6);
a_avg = mean(a, 'omitnan');
a1_n = length(a(~isnan(a(:,1))));
a2_n = length(a(~isnan(a(:,2))));
a_std = std(a, 'omitnan');
a_sem = a_std./[sqrt(a1_n), sqrt(a2_n)];
h1 = bar(1,a_avg(1),0.7,'FaceColor',CT(2,:),'FaceAlpha',.5, 'LineWidth',1);
hold on
h2 = bar(2,a_avg(2),0.7,'FaceColor',CT(6,:),'FaceAlpha',.5, 'LineWidth',1);
e = errorbar(a_avg,a_sem,'k', 'LineWidth',1.5);
e.LineStyle = 'none';
e.CapSize = 15;
set(gca,'XTick',0:3)
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,100,300,400])
xlim([0,3])
% ylim([0.4,1])
% name = {'','No stimulation','Light stimulation',''}
% set(gca,'xticklabel',name)
[h, p,~,stats] = ttest2(a(1:a1_n,1),a(1:a2_n,2))
% title('Control group')
end