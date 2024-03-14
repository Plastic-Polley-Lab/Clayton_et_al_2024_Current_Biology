



%% manually check the statbility of neuron firing

%%

latency_p2t = [spikedata.latency_p2t];

figure;
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
edges = 0.1:0.05:0.8;
h = histogram(latency_p2t,20, 'FaceColor','k');
hold on
plot([0.4, 0.4],[0,max(h.Values)],'--k','LineWidth', 1)
% xlim([0.2, 0.8])
% xticks([0.2:0.2:0.8])
ylabel('Number of units')
xlabel('Trough-peak interval (ms)')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,300])
hold off

%
indx = find(latency_p2t>0.4)
indx_thin = find(latency_p2t<0.4)
for i = 1:length(spikedata)
    waveform(:,i) = spikedata(i).waveform;
end

figure;
t = (1: 1: size(waveform(22:end,1),1))/fs * 1000;
plot(t,mean(waveform(22:end,indx),2)*1000,'Color', CT(6,:), 'LineWidth', 2)
hold on
plot(t,mean(waveform(22:end,indx_thin),2)*1000, 'Color', CT(2,:), 'LineWidth', 2)
ylabel('Amplitude (uV)')
xlabel('Time (ms)')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,300])
hold off


run_num = 1; % 9, 19, 69, 47, 55, 67
for i = 42:length(indx)
    psth(i) = spikedata(indx(i)).clusterData.run(run_num).psth;
    figure
    psth_plot(psth,i)
end

