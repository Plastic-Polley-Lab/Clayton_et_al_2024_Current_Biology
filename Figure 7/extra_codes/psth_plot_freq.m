function psth_plot_freq(psth1, psth2, window)
% INPUT:
%       psth: from psth_implae
%       chan: which channel to plot the fra
%       now: max for two colors
% title(num2str(chan))
spike = psth1;
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
h1 = subplot(3,1,1);
set(h1, 'Position', [.15 0.7 .7 .25]);
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
hold on
rectangle('Position',[0,0,3000,length(spike.raster)],'FaceColor', CT(colorIdx(1),:), 'EdgeColor', CT(colorIdx(1),:))
% imagesc(spike.scmatrix),colormap(BW)
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
ylim([0,length(spike.raster)])
axis off
box off
ylabel('Trial #')
% title(['Channel ', num2str(chan)])

spike = psth2;
h2 = subplot(3,1,2);
set(h2, 'Position', [.15 0.45 .7 .25]);
rectangle('Position',[0,0, 3000,length(spike.raster)],'FaceColor', CT(colorIdx(2),:), 'EdgeColor', CT(colorIdx(2),:))
% imagesc(spike.scmatrix),colormap(BW)
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
ylim([0,length(spike.raster)])
axis off
box off
ylabel('Trial #')

h3 = subplot(3,1,3);
set(h3, 'Position', [.15 0.2 .7 .25]);
spike = psth1
% plot(mean(spike.scmatrix,1))
% hold on
psth_avg = mean(spike.scmatrix,1)/0.001;
psth_avg_smooth = smoothts(psth_avg,'b',10);
hold on
plot(psth_avg_smooth, 'color', CT(colorIdx(1),:),'LineWidth', 1) % box smooth
hold on

spike = psth2
% plot(mean(spike.scmatrix,1))
% hold on
psth_avg = mean(spike.scmatrix,1)/0.001;
psth_avg_smooth = smoothts(psth_avg,'b',10);
hold on
plot(psth_avg_smooth, 'color', CT(colorIdx(2),:), 'LineWidth', 1) % box smooth

ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
xlim([window(1),window(end)])
ylim([0,max(psth_avg_smooth)/0.8])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,500,400])
hold off
% CT=cbrewer('div', 'RdYlBu', 6); % for nice color
% e = errorbar(mean(data,2),err, '-s','MarkerSize',6,'MarkerFaceColor', CT(6,:), 'Color',CT(6,:), 'LineWidth',1);
