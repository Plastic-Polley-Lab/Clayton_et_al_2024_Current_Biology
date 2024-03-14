% extract the latency from peak to trough
%% check waveform
figure
for i = 1:30
    plot(clusterData.run(1).spatiotemporalTemplate(:,i))
    hold on
end
%%
clear
file = dir('*.mat');
fs = 24414.0625
waveform =[];
latency_p2t = [];
for i = 1:length(file)
    load(file(i).name)
    template_indx = clusterData.run(1).templateMin;
    waveform(:,i) = clusterData.run(1).spatiotemporalTemplate(:,template_indx);
%     [~,peak] = max(waveform);
    [~,trough] = min(waveform(:,i));
    [~,peak] = max(waveform(trough:end,i));
    latency_p2t(i) = (peak-1)/fs * 1000;
end
%%
figure;
edges = 0.1:0.05:0.8;
histogram(latency_p2t,20)
indx = find(latency_p2t>0.4)
indx_thin = find(latency_p2t<0.4)
figure;
plot(mean(waveform(:,indx),2))
hold on
plot(mean(waveform(:,indx_thin),2))
