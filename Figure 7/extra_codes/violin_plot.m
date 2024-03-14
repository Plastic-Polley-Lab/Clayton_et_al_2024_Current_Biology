function violin_plot(raw_table)
% This code is only used for paired data
% Input:
%      a: n x 2 matrix, where each row represents a sample, a may contains
%      NaN
figure;
CT=cbrewer('div', 'RdYlBu', 6);
fig = violinplot(raw_table);
set(gca,'XTick',0:3)
fig(1).ViolinColor = CT(2,:);
fig(2).ViolinColor = CT(6,:);
for i = 1:2
    fig(i).ViolinAlpha = 0.5;
    fig(i).ScatterPlot.MarkerFaceAlpha = 1;
    fig(i).ViolinPlot.LineWidth = 1;
    fig(i).ViolinPlot.EdgeColor = [0,0,0];
    fig(i).BoxPlot.EdgeColor = [0,0,0];
    fig(i).BoxPlot.FaceColor = [0,0,0];
    fig(i).BoxWidth = 0.03;
    fig(i).WhiskerPlot.LineWidth =1;
    fig(i).WhiskerPlot.Color =[0,0,0];
    fig(i).MedianPlot.MarkerEdgeColor =[0,0,0];
end

box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,100,300,400])
% ylim([0.4,1])
% name = {'','No stimulation','Light stimulation',''}
% set(gca,'xticklabel',name)
% title('Control group')
end