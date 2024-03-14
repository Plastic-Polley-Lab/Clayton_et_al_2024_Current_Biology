function psth_plot_simple(spike, window)
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
plot(time,psth_avg_smooth, 'k', 'LineWidth', 1) % box smooth
ylabel('Firing rate (Hz)')
xlabel('Time (ms)')
xlim([window(1),window(end)])
ylim([0,max(psth_avg_smooth)/0.8])
box off
set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,200,500,400])
hold off