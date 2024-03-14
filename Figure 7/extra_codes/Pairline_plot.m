function Pairline_plot(a)
% This code is only used for paired data
% Input:
%      a: n x 2 matrix, where each row represents a sample, 
figure;
plot(a','-','Color','k', 'LineWidth',1)
hold on
CT=cbrewer('div', 'RdYlBu', 6);
bar(1,mean(a(:,1)),0.7,'FaceColor',CT(2,:),'FaceAlpha',.5, 'LineWidth',1);
hold on
bar(2,mean(a(:,2)),0.7,'FaceColor',CT(6,:),'FaceAlpha',.5, 'LineWidth',1);
scatter(ones(size(a(:,1))),a(:,1),36, CT(2,:),'filled')
hold on
scatter(2*ones(size(a(:,1))),a(:,2),36, CT(6,:),'filled')
hold on
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
[h, p,~,stats] = ttest(a(:,1),a(:,2))
% title('Control group')
end