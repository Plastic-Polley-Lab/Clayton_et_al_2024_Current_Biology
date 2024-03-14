function raster_plot_simple(psth,window)
spike = psth;
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
hold on
for i = 1:length(spike.stimulus.delay)
    rectangle('Position',[spike.stimulus.delay(i),length(spike.raster)+2*i-1, spike.stimulus.width(i),2],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
end
hold on

for i = 1:length(spike.raster)
    if ~isempty(spike.raster(i).ts)
        scatter(spike.raster(i).ts, i*ones(size(spike.raster(i).ts)), 6, '.','k')
        for j = 1:length(spike.raster(i).ts)
            plot([spike.raster(i).ts(j), spike.raster(i).ts(j)], [i-1,i], 'k', 'LineWidth',1)
        end
    end
end

hold off
xlim([window(1),window(end)])
% ylim([0,length(spike.raster)])
axis off
box off
ylabel('Trial #')
