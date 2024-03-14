function across_sets_plot(data_ana, fiberid)
cyclesize = 12
switch cyclesize
    case 4
        idx = 1
    case 12
        idx = 5
end
location = {'NDB', 'SI/GP'}
color = cbrewer('div', 'RdYlBu',4);
switch fiberid
    case 1
        col = 4
    case 2
        col = 1
end
for i = 1:6
    subplot(1,2,fiberid)
    
    temp_value(i,:) = data_ana(i).REG.cross_signal(idx,:);
    s1 = scatter(1:4,data_ana(i).REG.cross_signal(idx,:), 'o', 'MarkerFaceColor', color(col,:), 'MarkerEdgeColor',color(col,:))
    s1.MarkerFaceAlpha = .2;
    hold on
    h = plot(1:4, data_ana(i).REG.cross_signal(idx,:),'color' ,color(col,:),'LineWidth', 1)
    h.Color(4) = 0.2
    xlim([0,5])
    ylim([-0.1,0.4])
    xticks([1:4])
    xticklabels({'set1', 'set2', 'set3', 'Control'})
end
hold on
plot(1:4, mean(temp_value),'-d', 'color','k', 'MarkerFaceColor','k', 'LineWidth',1)
ylabel('Peak Cross-correlation')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,600])
title([location{fiberid}, ' Cycle- ', num2str(cyclesize)])
