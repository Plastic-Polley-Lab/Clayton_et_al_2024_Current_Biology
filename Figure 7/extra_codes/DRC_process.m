function drc_result = DRC_process(neuron_num, spikedata, threshold, SNR) 
%% DRC no noise 
% neuron_num = 3; % 20, 25, 33, 35 
% psth = spikedata(neuron_num).clusterData.psth; 
%% DRC no noise 
%by default there are 60 trials 
% for i = 1:3 
%     figure 
%     psth1.raster = psth.raster(i:3:60); 
%     psth1.scmatrix = psth.scmatrix(i:3:60,:); 
%     psth1.stimulus = psth.stimulus; 
%     psth_plot(psth1,1, 1:1250) 
% end 
%% DRC with noise 
% neuron_num = 21; 
% SNR = 20:-10:-10; 
innerIndexes = spikedata(neuron_num).clusterData.stimData.innerIndexes; 
inner_var = spikedata(neuron_num).clusterData.stimData.inner_variables; 
drc1_indx = []; 
 
%% summary DRC responses for different SNR and different set 
for drc_set = 1:3 
    psth_summary(drc_set) = drc_summary(drc_set, neuron_num, innerIndexes, spikedata, 0); 
end 
% suptitle(['Neuron', num2str(neuron_num)]) 
% set(gcf,'position',[100,50,700,1000]) 
% savefig(['DRC', num2str(1), '_response_', '.pdf']) 
% if neuron_num < 10 
%     saveas(gcf,['DRC', '_response_0', num2str(neuron_num), '.pdf']) 
% else 
%     saveas(gcf,['DRC', '_response_', num2str(neuron_num), '.pdf']) 
% end 
% close 
%% calculate the fano factor 
clear result 
for drc_set = 1:3 
    for i = 1:length(SNR) 
        %         figure 
        %         psth_plot(psth1.(['SNR',num2str(i)]),1, 1:1250, 1) 
        fanoparam.timeOnset = 250; % DRC starts at 250 ms, and each chord last 50 ms 
        fanoparam.timeWindow = 10; 
        result(drc_set).(['SNR',num2str(i)]) = VarVsMean(psth_summary(drc_set).(['SNR',num2str(i)]).scmatrix, fanoparam, 0); 
        %         box off 
        %         set(gca,'TickDir','out') 
        %         set(gca,'fontsize',12) 
        %         set(gca,'TickLengt', [0.015 0.015]); 
        %         set(gca, 'LineWidth',1) 
        %         set(gcf,'position',[100,200,500,400]) 
    end 
end 
 
%% plot an example SNR for drc_set == 1 
% CT=cbrewer('div', 'RdYlBu', 6); % for nice color 
% colorIdx = [6,5,4]; 
% i = 1 % different SNR: 1, 3, 5, 7 
% figure 
% for drc_set = 1:3 
%     scatter(result(drc_set).(['SNR',num2str(i)]). spikecount_m, result(drc_set).(['SNR',num2str(i)]). spikecount_var, [], CT(colorIdx(drc_set),:), 'filled') 
%     hold on 
% end 
% 
% 
% spikecount_m = [result(1).(['SNR',num2str(i)]).spikecount_m, ... 
%     result(2).(['SNR',num2str(i)]).spikecount_m, ... 
%     result(3).(['SNR',num2str(i)]).spikecount_m]; 
% spikecount_var = [result(1).(['SNR',num2str(i)]).spikecount_var, ... 
%     result(2).(['SNR',num2str(i)]).spikecount_var, ... 
%     result(3).(['SNR',num2str(i)]).spikecount_var]; 
% mdl = fitlm(spikecount_m, spikecount_var, 'Intercept',false); 
% fanofactor = mdl.Coefficients.Estimate 
% 
% hold on 
% h2 = plot(spikecount_m, mdl.Fitted, 'r') 
% text(mean(spikecount_m), mean(mdl.Fitted), ['Slope = ', num2str(fanofactor)], 'fontsize', 12) 
% 
% xlabel('Mean Spike Count') 
% ylabel('Spike Count Variance') 
% box off 
% set(gca,'TickDir','out') 
% set(gca,'fontsize',12) 
% set(gca,'TickLengt', [0.015 0.015]); 
% set(gca, 'LineWidth',1) 
% set(gcf,'position',[100,200,500,400])%% calculate the cross-correlation 
% title(['SNR ',num2str(SNR(i)), ' dB']) 
 
%% calculate the fanofactor for the neuron 
clear fanofactor 
for i = 1: length(SNR) 
    spikecount_m = [result(1).(['SNR',num2str(i)]).spikecount_m, ... 
        result(2).(['SNR',num2str(i)]).spikecount_m, ... 
        result(3).(['SNR',num2str(i)]).spikecount_m]; 
    spikecount_var = [result(1).(['SNR',num2str(i)]).spikecount_var, ... 
        result(2).(['SNR',num2str(i)]).spikecount_var, ... 
        result(3).(['SNR',num2str(i)]).spikecount_var]; 
    mdl = fitlm(spikecount_m, spikecount_var, 'Intercept',false); 
    fanofactor(i).fanofactor = mdl.Coefficients.Estimate; 
    fanofactor(i).spikecount_m = spikecount_m; 
    fanofactor(i).spikecount_var = spikecount_var; 
end 
% figure; 
% plot(1:length(fanofactor), [fanofactor.fanofactor],'-k') 
% hold on 
% scatter(1:length(fanofactor), [fanofactor.fanofactor],'ok','filled') 
% xticks(1:length(fanofactor)) 
% for i = 1: 7 
%     labels{i} = num2str(SNR(i)) 
% end 
% xticklabels(labels) 
% xlabel('SNR (dB)') 
% ylabel('Fano Factor') 
% title(['Neuron', num2str(neuron_num)]) 
 
% ylim([0.5, 1]) 
% box off 
% set(gca,'TickDir','out') 
% set(gca,'fontsize',12) 
% set(gca,'TickLengt', [0.015 0.015]); 
% set(gca, 'LineWidth',1) 
% set(gcf,'position',[100,200,500,400]) 
% if neuron_num < 10 
%     saveas(gcf,['DRC', '_fanofactor_0', num2str(neuron_num), '.pdf']) 
% else 
%     saveas(gcf,['DRC', '_fanofactor_', num2str(neuron_num), '.pdf']) 
% end 
% close 
 
%% Use cross-correlation to calculate the temporal precision 
for drc_set = 1:3 
    for i = 1:length(SNR) 
         
        crossparam.timeOnset = 250; % DRC starts at 250 ms, and each chord last 50 ms 
        crossparam.timeWindow = 10; 
        result(drc_set).(['SNR',num2str(i)]).CorrCoef = cross_coeff(psth_summary(drc_set).(['SNR',num2str(i)]).scmatrix, crossparam, 0); 
         
    end 
end 
for i = 1: length(SNR) 
    corrR(i) = mean([result(1).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ... 
        result(2).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ... 
        result(3).(['SNR',num2str(i)]).CorrCoef.corrR_avg]); 
     
     
end 
 
%% save the data 
drc_result.psth_summary = psth_summary; 
drc_result.result       = result; 
drc_result.fanofactor   = fanofactor; 
drc_result.corrR        = corrR; 
drc_result.latency_p2t  = spikedata(neuron_num).latency_p2t; 
drc_result.keep         = spikedata(neuron_num).keep; 
% drc_result.example      = spikedata(neuron_num).example; 
if drc_result.corrR(1) > threshold 
    drc_result.resp = 1; 
else 
    drc_result.resp = 0; 
end 
end 
 