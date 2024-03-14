function plot_patternSeq(t, CycleLength)
figure

color =cbrewer('div', 'RdYlBu', CycleLength); % for nice color;
% color = lines(length(C));
% gscatter(1:length(t),t,group',color)

for i = 1:CycleLength
    scatter(i:CycleLength:size(t,2), t(1,i:CycleLength:end),'o','filled','MarkerEdgeColor', color(i,:), 'MarkerFaceColor', color(i,:))
    ylim([50,350])
    ylabel('Interval (ms)')
    hold on
end


xlabel('Sound sequence #')
set(gcf,'position',[100,200,1200,200])
