function psth_plot_tosca(psth, neuron_num, fig)
% INPUT:
%       psth: from spike2PSTH

    
spike = psth;
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
h1 = subplot(2,1,1);
set(h1, 'Position', [.15 0.55 .7 .35]);
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
% rectangle('Position',[100,length(spike.spikeraster)+2, 0.05 ,2],'FaceColor', CT(colorIdx(1),:), 'EdgeColor', CT(colorIdx(1),:))
hold on
if fig == 0
    imagesc(spike.timepoint, [], spike.scmatrix),colormap(BW)
else
    for i = 1:length(spike.spikeraster)
        if ~isempty(spike.spikeraster(i).times)
            scatter(spike.spikeraster(i).times, i*ones(size(spike.spikeraster(i).times)), 6, '.','k')
            hold on
            for j = 1:length(spike.spikeraster(i).times)
                plot([spike.spikeraster(i).times(j), spike.spikeraster(i).times(j)], [i-1,i], 'k', 'LineWidth',1)
            end
        end
    end
end
hold off
% xlim([window(1),window(end)])
% ylim([0,length(spike.raster)])
axis off
box off
ylabel('Trial #')
title(['Neuron ', num2str(neuron_num)])

h2 = subplot(2,1,2);
set(h2, 'Position', [.15 0.2 .7 .35]);
psth_avg_smooth = smoothts(spike.FR_avg,'b',5);
hold on
rectangle('Position',[0 ,0, 0.05,max(psth_avg_smooth)],'FaceColor', CT(colorIdx(1),:), 'EdgeColor', CT(colorIdx(1),:))

plot(spike.timepoint,psth_avg_smooth, 'k', 'LineWidth', 1) % box smooth


ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
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
