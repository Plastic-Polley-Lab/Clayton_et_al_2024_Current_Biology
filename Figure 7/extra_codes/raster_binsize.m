function [spike_scmatrix, psth_avg, psth_avg_smooth, time] = raster_binsize(spike, window, window_step)

% spike = summarydata.wt(1).clusterData.psth;
for i = 1:length(window)
    spike_scmatrix(:,i) = sum(spike.scmatrix(:,window(i): window(i)+ window_step-1),2); %
end
psth_avg = mean(spike_scmatrix,1)/(window_step/1000);
psth_avg_smooth = smoothts(psth_avg,'b',5);
time = window- window(1) + window_step/2;