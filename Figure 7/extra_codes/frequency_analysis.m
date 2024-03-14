%% frequency pattern
neuron_num = 2;
pattern_idx = 1:2:length(spikedata(1).clusterData.psth.raster);
figure
psth_idx = pattern_idx;
psth = spikedata(neuron_num).clusterData.psth;
psth1.raster = psth.raster(psth_idx );
psth1.scmatrix = psth.scmatrix(psth_idx ,:);
psth1.stimulus = psth.stimulus;

pattern_idx = 2:2:length(spikedata(1).clusterData.psth.raster);
figure
psth_idx = pattern_idx;
psth = spikedata(neuron_num).clusterData.psth;
psth2.raster = psth.raster(psth_idx );
psth2.scmatrix = psth.scmatrix(psth_idx ,:);
psth2.stimulus = psth.stimulus;

% plot the data
psth_plot_freq(psth1,psth2,1:3000)
set(gcf,'position',[100,200,1000,400])
%% noise pattern
neuron_num = 2
spike = spikedata(neuron_num).clusterData.psth;
window = 1:100;
figure;
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
colorIdx = [2,6];
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
hold on
for i = 1:length(spike.stimulus.delay)
    rectangle('Position',[spike.stimulus.delay(i),0, spike.stimulus.width(i),length(spike.raster)],'FaceColor', CT(colorIdx(i),:), 'EdgeColor', CT(colorIdx(i),:))
end
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
axis on
ylabel('Trial #')
xlabel('Time (ms)')