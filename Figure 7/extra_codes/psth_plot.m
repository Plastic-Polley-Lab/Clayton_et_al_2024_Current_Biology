function psth_plot(psth, neuron_num, window, fig, chan)
% INPUT:
%       psth: from psth_implae
%       chan: which channel to plot the fra
%       window: time range for plotting, support different binsize: for
%       instance: 1: 2:1000; 
%       fig: 1 for high resolution; 0 for low resolution plots
if nargin == 5 % used for Impale plotting from the main gui
    title(num2str(chan))
    spike = psth(chan);
else
    spike = psth;
    
end
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
h1 = subplot(2,1,1);
set(h1, 'Position', [.15 0.55 .7 .35]);
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
hold on
for i = 1:length(spike.stimulus.delay)
    rectangle('Position',[spike.stimulus.delay(i),length(spike.raster)+2*i-1, spike.stimulus.width(i),2],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
end
hold on
if fig == 0
    imagesc(spike.scmatrix),colormap(BW)
else
    for i = 1:length(spike.raster)
        if ~isempty(spike.raster(i).ts)
            scatter(spike.raster(i).ts, i*ones(size(spike.raster(i).ts)), 6, '.','k')
            for j = 1:length(spike.raster(i).ts)
                plot([spike.raster(i).ts(j), spike.raster(i).ts(j)], [i-1,i], 'k', 'LineWidth',1)
            end
        end
    end
end
hold off
xlim([window(1),window(end)])
% ylim([0,length(spike.raster)])
axis off
box off
ylabel('Trial #')
title(['Neuron ', num2str(neuron_num)])

h2 = subplot(2,1,2);
set(h2, 'Position', [.15 0.2 .7 .35]);
% plot(mean(spike.scmatrix,1))
% hold on
window_step = unique(diff(window));
if window_step == 1
    psth_avg = mean(spike.scmatrix,1)/0.001;
    psth_avg_smooth = smoothts(psth_avg,'g',10,4);
    time = window - window(1)+ window_step/2;
else
    for i = 1:length(window)
        spike_scmatrix(:,i) = sum(spike.scmatrix(:,window(i): window(i)+ window_step-1),2); %
    end
    psth_avg = mean(spike_scmatrix,1)/(window_step/1000);
    psth_avg_smooth = smoothts(psth_avg,'b',5);
    time = window- window(1) + window_step/2; 
end
hold on
% for i = 1:length(spike.stimulus.delay)
%     rectangle('Position',[spike.stimulus.delay(i),0, spike.stimulus.width(i), max(psth_avg_smooth)/0.8],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
% end
% if isfield(psth, 'timepoint')
%     time = psth.timepoint;
% end
plot(time,psth_avg_smooth, 'k', 'LineWidth', 1) % box smooth
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
