function drc_plot(psth, chan, window, k)
% INPUT:
%       psth: from psth_implae
%       chan: which channel to plot the fra
%       now: max for two colors
% title(num2str(chan))
SNR = 20:-5:-10;
spike = psth(chan);
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
subaxis(14,1,2*(k-1) +1 , 'sv', 0.01);
% set(h1, 'Position', [.15 0.8 - (k-1)/10 .7 .1]);
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
hold on
% for i = 1:length(spike.stimulus.delay)
     i = 1; % only plot drc stimulation
    rectangle('Position',[spike.stimulus.delay(i),0, spike.stimulus.width(i),length(spike.raster)],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
% end
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

subaxis(14,1,2*(k-1) +2 , 'sv', 0.01);
% set(h2, 'Position', [.15 0.6-(k-1)/10 .7 .1]);
% plot(mean(spike.scmatrix,1))
% hold on
psth_avg = mean(spike.scmatrix,1)/0.001;
psth_avg_smooth = smoothts(psth_avg,'b',10);
hold on
% for i = 1:length(spike.stimulus.delay)
    i = 1
    rectangle('Position',[spike.stimulus.delay(i),0, spike.stimulus.width(i), max(psth_avg_smooth)/0.8],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
% end
plot(psth_avg_smooth, 'k', 'LineWidth', 1) % box smooth
% ylim([0, 200])

if k == 4
    ylabel('Firing rate (Hz)')
    set(gca,'xtick',[])
    yyaxis right
    set(gca,'ytick',[])
    ylabel(['SNR ', num2str(SNR(k)), ' dB'],'FontWeight','bold','Color','k')
elseif k == 7
    yyaxis right
    set(gca,'ytick',[])
    ylabel(num2str(SNR(k)), 'FontWeight','bold','Color','k')
    xlabel('Time (ms)')

else
    set(gca,'xtick',[])
    yyaxis right
    set(gca,'ytick',[])
    ylabel(num2str(SNR(k)), 'FontWeight','bold','Color','k')
end
xlim([window(1),window(end)])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,500,400])
hold off
% CT=cbrewer('div', 'RdYlBu', 6); % for nice color
% e = errorbar(mean(data,2),err, '-s','MarkerSize',6,'MarkerFaceColor', CT(6,:), 'Color',CT(6,:), 'LineWidth',1);
