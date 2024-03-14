
%% signal correlation
% neuron1 = fra_summary(2);
% neuron2 = fra_summary(3);
% figure;
% scatter(mean(neuron1.tuning, 2), mean(neuron2.tuning, 2),'ko','filled')
% xlim([0,20])
% ylim([0,20])
% xlabel('Unit 6 (spike count)')
% ylabel('Unit 12 (spike count)')
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,100,400,400])
% R = corrcoef(mean(neuron1.tuning, 2),mean(neuron2.tuning, 2));
% R = R(1,2);

%% Noise correlation
% example

% figure;
% scatter(neuron1.fra(7,2,:), neuron2.fra(7,2,:),'ko','filled')
% xlim([0,10])
% ylim([0,10])
% xlabel('Unit 6 (spike count)')
% ylabel('Unit 12 (spike count)')
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,100,400,400])
% R = corrcoef(neuron1.fra(6,1,:), neuron2.fra(6,1,:));
% R = R(1,2);
%
% R_noise =[];
% for i = 1: size(neuron1.fra,1)
%     for j = 1:size(neuron1.fra,2)
%        temp = corrcoef(neuron1.fra(i,j,:), neuron2.fra(i,j,:));
%         R_noise(i,j) = temp(1,2);
%     end
% end
%
% psth = spikedata(neuron_num).clusterData.psth;
% rep  = spikedata(1).clusterData.stimData.numReps;
%% by default there are 60 trials
% figure
% psth_plot(psth,1, 1:750)

%% get the tone responsive neurons and analyze the FRAs
clear;close all;
path = pwd;
FRA_analysis_preprocess(path)
%% summary analysis
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\FRA_analysis')
genotype = {'WT', 'KO'}
for z = 1:length(genotype)
    
    switch genotype{z}
        case 'WT'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091820\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\092020\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110620\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110720\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC36\121520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC40\122620\FRA_CF'};
        case 'KO'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091620\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092320\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092720\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092820\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\103020\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC41\122920\FRA_CF'};
    end
    temp =[];
    for i = 1 : length(paths)
        load([paths{i}, '\summary_fra_responses_v2.mat'])
        for j = 1:length(fra_summary)
            fra_summary(j).paths = paths{i};
        end
        temp = [temp, fra_summary];
    end
    save(['fra_tuning_results_', genotype{z}, '.mat'],'temp', '-v7.3')
end
%% add new collected data
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\FRA_analysis')
genotype = {'WT', 'KO'}
for z = 1:length(genotype)
    
    switch genotype{z}
        case 'WT'
            paths = {
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC36\121520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC40\122620\FRA_CF'};
        case 'KO'
            paths = {
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC41\122920\FRA_CF'};
    end
    load(['E:\Ke_Chen\Processed Data\PTCHD1-Project\', 'fra_tuning_results_', genotype{z}, '_v2.mat' ])
    % temp = [];
    for i = 1 : length(paths)
        load([paths{i}, '\summary_fra_responses_v2.mat'])
        fra_summary = rmfield(fra_summary, 'keep');
        for j = 1:length(fra_summary)
            fra_summary(j).paths = paths{i};
        end
        temp = [temp, fra_summary];
    end
    save(['fra_tuning_results_', genotype{z}, '_v2.mat'],'temp', '-v7.3')
end
%% get the rlf from more than the best frequency
% edit rlf_multifreq
% clear
% load('fra_tuning_results_WT_v2.mat')
% for i = 1:length(temp)
%     [wt.rlf(i,:), wt.rlf_multifreqs(i,:), wt.spont(i), wt.fra{i}, wt.spls{i}, wt.freqs{i}] = rlf_multifreq(temp(i));
% end
% load('fra_tuning_results_KO_v2.mat')
% for i = 1:length(temp)
%     [ko.rlf(i,:), ko.rlf_multifreqs(i,:),ko.spont(i), ko.fra{i}, ko.spls{i}, ko.freqs{i}] = rlf_multifreq(temp(i));
% end
% clear temp
% % see the normalized ata
% plot_norm(wt.rlf, ko.rlf)
% plot_norm(wt.rlf_multifreqs, ko.rlf_multifreqs)
%% Let's fit the rlf with two gussian
% [data_model, r_square, spls_interp, data_interp] = rlf_model_fitting(data);
% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\FRA_analysis')
% load('fra_tuning_results_WT_v2.mat')
% for i = 1:length(temp)
%     wt.rlf(i,:) = mean(temp(i).rlf, 2); 
%     wt.spont(i) = temp(i).spont;
%     data = wt.rlf(i,:);
%     [wt.data_model{i}, wt.r_square(i), wt.best(i), wt.spls_interp{i}, wt.data_interp{i}] = rlf_model_fitting(data);
%     wt.MI(i) = (wt.data_model{i}(end) - wt.spont(i))/(max(wt.data_model{i}) - wt.spont(i));
%     %     rlf_min = min(wt.rlf(i,:));
% %     rlf_max = max(wt.rlf(i,:));
% %     wt.rlf_normal(i,:) = (wt.rlf(i,:)-rlf_min)/rlf_min;
% end

% load('fra_tuning_results_KO_v2.mat')
% for i = 1:length(temp)
%     ko.rlf(i,:) = mean(temp(i).rlf, 2);%% here rlf is divided by 4; because in the code, spike count was divided by 0.025 instead of 0.1
%     ko.spont(i) = temp(i).spont;
%     data = ko.rlf(i,:);
%     [ko.data_model{i}, ko.r_square(i), ko.best(i), ko.spls_interp{i}, ko.data_interp{i}] = rlf_model_fitting(data);
%     ko.MI(i) = (ko.data_model{i}(end) - ko.spont(i))/(max(ko.data_model{i}) - ko.spont(i));
% end
%%
%%%%%%%Let's see neurons with good fit%%%%%%%%%%%%%
% wt_ind_good = find([wt.r_square]>=0.7);
% ko_ind_good = find([ko.r_square]>=0.7);
% % normalized the firing rate
% for i = 1:length(wt_ind_good)
%     ind = wt_ind_good(i);
%     temp_n = wt.data_model{ind} - wt.spont(ind); % get the driving firing rate
%     wt.data_norm{i} = temp_n./max(temp_n);
%     wt.threshold(i) = min(find(wt.data_norm{i}>=0.2))-1;
%     [~, ind_best] = max(wt.data_model{ind}); 
%     wt.gain(i)      = (max(wt.data_model{ind}) - wt.data_model{ind}(wt.threshold(i) +1))/(ind_best-wt.threshold(i));
% end

% for i = 1:length(ko_ind_good)
%     ind = ko_ind_good(i);
%     temp_n = ko.data_model{ind} - ko.spont(ind); % get the driving firing rate
%     ko.data_norm{i} = temp_n./max(temp_n);
%     ko.threshold(i) = min(find(ko.data_norm{i}>=0.2))-1;
%     [~, ind_best] = max(ko.data_model{ind}); 
%     ko.gain(i)      = (max(ko.data_model{ind}) - ko.data_model{ind}(ko.threshold(i) +1))/(ind_best-ko.threshold(i));
% end
freqs = tuning_summary.wt(1).freqs;
[gain_first.ko, gain_sec.ko, ~, ~] = check_drifts_FRA(results.ko, freqs);
[gain_first.wt, gain_sec.wt, ~, ~] = check_drifts_FRA(results.wt, freqs);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the data


%% Let's see the population gain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Let's see the population gain
load('fra_tuning_summary.mat')

results.wt = getGain_threshold(tuning_summary.wt);
results.ko = getGain_threshold(tuning_summary.ko);
%% check the data collected at different data collection period
[gain_first.wt, gain_sec.wt, ~, ~] = check_drifts_FRA(results.wt);
[gain_first.ko, gain_sec.ko, ~, ~] = check_drifts_FRA(results.ko);
%% plot the Monotonic index
% [fig1, fig2] = ecdf_bar_plot([results.wt.MI], [results.ko.MI], 'Monotonic Index');
% set(get(fig1,'XLabel'), 'String', 'Monotonic Index');
% set(get(fig2,'YLabel'), 'String', 'Monotonic Index');
% ylim(fig2,[0,1])

%% plot the gains
% figure
% [fig1, fig2] = ecdf_bar_plot([results.wt.gains_value], [results.ko.gains_value], 'Gains');
% set(get(fig1,'XLabel'), 'String', 'Gains (\Delta sp/s per 10 dB step)');
% set(get(fig2,'YLabel'), 'String', 'Gains (\Delta sp/s per 10 dB step)');
% xlim(fig1, [0,50])
% ylim(fig2, [0,15])
%% plot the monotonic index with the gains
% figure
% h1 = scatter([results.wt.MI], [results.wt.gains_value], 'ko', 'filled');
% h1.MarkerFaceAlpha = 0.3;
% hold on
% h2 = scatter([results.ko.MI], [results.ko.gains_value], 'ro', 'filled');
% h2.MarkerFaceAlpha = 0.3
% xlabel('Monotonic Index')
% ylabel('Gains (\Delta sp/s per 10 dB step)')
% box off
% legend([h1, h2], {'WT', 'KO'})
% set(gcf, 'Color', 'w')
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,200,400,400])
% hold off

%% Let's see the population gain
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    for i = 1:length(results.(genotype{j}))
%         plots.(genotype{j}).rlf(i,:) = results.(genotype{j})(i).rlf_avg;
        plots.(genotype{j}).rlf(i,:) = results.(genotype{j})(i).rlf_avg - ...
            results.(genotype{j})(i).spont;

    end
    plots.(genotype{j}).gains       = [results.(genotype{j}).gains_value];
    plots.(genotype{j}).threshold   = [results.(genotype{j}).threshold];
    plots.(genotype{j}).latency_p2t = [results.(genotype{j}).latency_p2t];
    
    nan_indx = find(isnan(plots.(genotype{j}).threshold));
    
    plots.(genotype{j}).gains(nan_indx)=[];
    plots.(genotype{j}).threshold(nan_indx) =[];
    plots.(genotype{j}).rlf(nan_indx,:) = [];
    plots.(genotype{j}).latency_p2t(nan_indx) = [];
end

% figure;
% 
% line_errorbar_drc(plots.wt.rlf, plots.ko.rlf)
% xticks([1:8])
% xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
% xlabel('Sound Intensity(dB)')
% ylabel('Firing Rate (Hz)')
% set(gcf,'position',[100,200,400,400])
% set(gcf, 'Color', 'w')
% export_fig('Rate_Level_Function',  '-png', '-pdf')

%% Let's see the threshold
% ecdf_bar_plot(plots.wt.threshold', plots.ko.threshold')
% figure(101)
% ylim([0,50])
% xticks([1,2])
% xticklabels({'WT', 'KO'})
% ylabel('Best Frequence Response Threshold (dB)')
% set(gcf, 'Color', 'w')
% export_fig('Response_Threshold',  '-png')
%% Let's check the gains
% ecdf_bar_plot(plots.wt.gains, plots.ko.gains)
% figure(100)
% xlabel('Gains (\Delta sp/s per 10 dB step)')
% xlim([0,50])
% figure(101)
% ylabel('Gains (\Delta sp/s per 10 dB step)')
% ylim([0,20])

% export_fig('Gains',  '-png', '-pdf')
%% try cluster the rate-level function
wt_indx = 1:size(plots.wt.rlf, 1);
ko_indx = (wt_indx(end) + 1): (wt_indx(end) + size(plots.ko.rlf,1));

alldata = [plots.wt.rlf; plots.ko.rlf];
alldata_n = alldata./max(alldata, [], 2);
clusterN = 6;
[T, idx] = pca_cluster(alldata_n, [], clusterN);

wt_cluster = T(wt_indx);
ko_cluster = T(ko_indx);
% plot the cluster proportion
for i = 1:clusterN
    wt_temp = find(wt_cluster == i);
    ko_temp = find(ko_cluster == i);
    plots.wt.clusterRatio(i) = length(wt_temp);
    plots.ko.clusterRatio(i) = length(ko_temp);
    plots.wt.cluster_gains{i} = plots.wt.gains(wt_temp);
    plots.ko.cluster_gains{i} = plots.ko.gains(ko_temp);
    plots.wt.cluster_rlf{i} = plots.wt.rlf(wt_temp,:);
    plots.ko.cluster_rlf{i} = plots.ko.rlf(ko_temp,:);    
end
%% plot the cluster results
figure;
h1 = plot(plots.wt.clusterRatio/sum(plots.wt.clusterRatio), '-o', 'LineWidth', 1)
hold on
h2 = plot(plots.ko.clusterRatio/sum(plots.ko.clusterRatio), '-o', 'LineWidth', 1)

h1.Color = 'k';
h1.MarkerFaceColor = [0.5, 0.5, 0.5];
h2.Color = 'r';
h2.MarkerFaceColor = 'r';

legend([h1, h2], 'WT', 'KO')
xticks([1:1:6])
% xticklabels({'4k', '5.6k', '8k', '11.3k', '16k', '22.6k', '32k', '45.3k', '64k'})
xticklabels({'1',  '2',  '3',  '4', '5', '6'})

xlabel('Cluter #')
ylabel('Proportion')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
hold off

box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
%% plot the average for each clustter
for i = 1:clusterN
    figure
    line_errorbar_drc(plots.wt.cluster_rlf{i}, plots.ko.cluster_rlf{i})
    xticks([1:8])
    xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
    xlabel('Sound Intensity(dB)')
    ylabel('Firing Rate (Hz)')
    set(gcf,'position',[100,200,400,400])
    set(gcf, 'Color', 'w')
    title(['Cluster ', num2str(i)])
    
end
%% separate the neuron based on trough-peak interval
fs_cutoff = 0.5
pyr_cutoff = 0.6;
plots.wt.fs_indx = find(plots.wt.latency_p2t<fs_cutoff);
plots.wt.pyr_indx = find(plots.wt.latency_p2t>pyr_cutoff);
plots.wt.fs_gains = plots.wt.gains(plots.wt.fs_indx);
plots.wt.pyr_gains = plots.wt.gains(plots.wt.pyr_indx);

plots.ko.fs_indx = find(plots.ko.latency_p2t<fs_cutoff);
plots.ko.pyr_indx = find(plots.ko.latency_p2t>pyr_cutoff);
plots.ko.fs_gains = plots.ko.gains(plots.ko.fs_indx);
plots.ko.pyr_gains = plots.ko.gains(plots.ko.pyr_indx);
%% plot the gains for putative FS
ecdf_bar_plot(plots.wt.fs_gains,plots.ko.fs_gains)
figure(100)
xlabel('Gains (\Delta sp/s per 10 dB step)')
xlim([0,60])
figure(101)
ylabel('Gains (\Delta sp/s per 10 dB step)')
%% plot the gains for putative Pyramidl neurons
ecdf_bar_plot(plots.wt.pyr_gains,plots.ko.pyr_gains)
figure(100)
xlabel('Gains (\Delta sp/s per 10 dB step)')

figure(101)
ylabel('Gains (\Delta sp/s per 10 dB step)')
ylim([0,20])
%% summary the analysis
save('fra_tuning_summary.mat', 'tuning_summary', 'plots', 'results', '-v7.3')

%%
% xlswrite('Gain_wt.xlsx', wt.gains_value)
% xlswrite('Gain_ko.xlsx', ko.gains_value)
% save('fra_tuning_gain_summary.mat', 'wt', 'ko')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||%
%||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||% 
load('fra_tuning_summary.mat')
summary_plot_gain_tuning(results, plots)

%% Let's see the frequency tuning
clear
load('fra_tuning_results_WT_v2.mat')
for i = 1:length(temp)
    [temp(i).cFRA, temp(i).spontDistr, temp(i).fraDistr, temp(i).dprime, temp(i).FRAmask] = dprime_FRA(temp(i).fra);
end
tuning_summary.wt = temp;

load('fra_tuning_results_KO_v2.mat')
for i = 1:length(temp)
    [temp(i).cFRA, temp(i).spontDistr, temp(i).fraDistr, temp(i).dprime, temp(i).FRAmask] = dprime_FRA(temp(i).fra);
end
tuning_summary.ko = temp;

save('fra_tuning_summary.mat', 'tuning_summary', '-v7.3')
%%
% summary_data_wt = [tuning_summary.wt.dprime];
% summary_data_ko = [tuning_summary.ko.dprime];
% [f_wt, x_wt] = ecdf(summary_data_wt);
% [f_ko, x_ko] = ecdf(summary_data_ko);
% figure;
% h1 = plot(x_wt, f_wt, '-k', 'LineWidth', 1);
% hold on
% h2 = plot(x_ko, f_ko, '-r', 'LineWidth', 1);
% xlabel('D-prime')
% ylabel('Cumulative probablity')
% legend([h1, h2], {'WT', 'KO'})
% box off
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,200,400,400])
% 
% set(gcf, 'Color', 'w')
% export_fig('Tuning_Dprime_ecdf_MGB',  '-png')
% 
% 
% a=summary_data_wt';
% b = summary_data_ko';
% b(end+1:length(a)) =  NaN;
% figure
% [h1, h2] = bar_plot([a, b])
% h1.FaceColor = 'k';
% h2.FaceColor = 'r';
% h1.FaceAlpha = 0.3;
% h2.FaceAlpha = 0.7;
% xticks([1,2])
% xticklabels({'WT', 'KO'})
% ylabel('D-prime')
% set(gcf, 'Color', 'w')
% export_fig('D-prime_bar_MGB',  '-png')
%% update the tuning_summary and plots to plot frequency tuning refered to BF
load('fra_tuning_summary.mat')
[plots, tuning_summary] = tuning_refer_BF(tuning_summary, plots, results);

save('fra_tuning_summary.mat', 'tuning_summary', 'results', 'plots', '-v7.3')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the summary Data
load('fra_tuning_summary.mat')
summary_plot_gain_tuning(results, plots)
%%
%% organize the data for plots
genotype = {'wt', 'ko'};
for j = 1:length(genotype)


    plots.(genotype{j}).MI = [results.(genotype{j}).MI];
    nan_indx = find(isnan([results.(genotype{j}).threshold]));
    

    plots.(genotype{j}).MI(nan_indx) =[];
end

% wt_indx = find([tuning_summary.wt.dprime] > 4);
% wt_good = tuning_summary.wt(wt_indx);
% 
% ko_indx = find([tuning_summary.ko.dprime] > 4);
% ko_good = tuning_summary.ko(ko_indx);
% 
% for i = 1:length(wt_good)
%     wt_good(i).norm_tuning = mean(wt_good(i).tuning,2)./max(mean(wt_good(i).tuning,2));
% end
% for i = 1:length(ko_good)
%     ko_good(i).norm_tuning = mean(ko_good(i).tuning,2)./max(mean(ko_good(i).tuning,2));
% end
%% create the matrix to store the normalized tuning curve, where best
% frequency is centered
% tuning_best_norm.wt = NaN(length(wt_good), 17);
% tuning_best_norm.wt_n = NaN(length(wt_good), 17);
% 
% for i = 1:length(wt_good)
%     [~,I] = max(wt_good(i).norm_tuning);
%     tuning_best_norm.wt_n(i, (10-I):(18-I)) = wt_good(i).norm_tuning; % center the best frequency on 9th column
%     tuning_best_norm.wt(i, (10-I):(18-I)) = mean(wt_good(i).tuning,2); % center the best frequency on 9th column
%     
% end
% 
% tuning_best_norm.ko = NaN(length(ko_good), 17);
% tuning_best_norm.ko_n = NaN(length(ko_good), 17);
% 
% for i = 1:length(ko_good)
%     [~,I] = max(ko_good(i).norm_tuning);
%     tuning_best_norm.ko_n(i, (10-I):(18-I)) = ko_good(i).norm_tuning; % center the best frequency on 9th column
%     tuning_best_norm.ko(i, (10-I):(18-I)) = mean(ko_good(i).tuning,2); % center the best frequency on 9th column
%     
% end
%%
% figure;
% for i = 3:size(tuning_best_norm.wt, 2)-2
%     temp = tuning_best_norm.wt(:,i);
%     temp(isnan(temp))=[];
%     tuning_best_norm.wt_sum{i} = temp;
%     
%     temp = tuning_best_norm.wt_n(:,i);
%     temp(isnan(temp))=[];
%     tuning_best_norm.wt_n_sum{i} = temp;
%     
%     temp = tuning_best_norm.ko_n(:,i);
%     temp(isnan(temp))=[];
%     tuning_best_norm.ko_n_sum{i} = temp;
%     
%     temp = tuning_best_norm.ko(:,i);
%     temp(isnan(temp))=[];
%     tuning_best_norm.ko_sum{i} = temp;  
% end
%% 
% %%
% for i = 3:length(tuning_best_norm.wt_sum)
%     plot_data_m(i-2) = mean(tuning_best_norm.wt_sum{i});
%     plot_data_sem(i-2) = std(tuning_best_norm.wt_sum{i})/sqrt(length(tuning_best_norm.wt_sum{i}));
%     plot_data_m_ko(i-2) = mean(tuning_best_norm.ko_sum{i});
%     plot_data_sem_ko(i-2) = std(tuning_best_norm.ko_sum{i})/sqrt(length(tuning_best_norm.ko_sum{i}));
%     
%     plot_data_m_n(i-2) = mean(tuning_best_norm.wt_n_sum{i});
%     plot_data_sem_n(i-2) = std(tuning_best_norm.wt_n_sum{i})/sqrt(length(tuning_best_norm.wt_n_sum{i}));
%     plot_data_m_ko_n(i-2) = mean(tuning_best_norm.ko_n_sum{i});
%     plot_data_sem_ko_n(i-2) = std(tuning_best_norm.ko_n_sum{i})/sqrt(length(tuning_best_norm.ko_n_sum{i}));
%     
% end
% figure
% h1 = errorbar(plot_data_m, plot_data_sem, '-ko');
% h1 = errorbar(plot_data_m_n, plot_data_sem_n, '-ko');
% set(h1, 'LineWidth', 1)
% set(h1, 'MarkerFaceColor', 'k')
% hold on
% h2 = errorbar(plot_data_m_ko, plot_data_sem_ko, '-ro');
% % h2 = errorbar(plot_data_m_ko_n, plot_data_sem_ko_n, '-ro');
% set(h2, 'LineWidth', 1)
% set(h2, 'MarkerFaceColor', 'r')
% legend([h1, h2], {'WT', 'KO'})
% ylabel('Normalized Firing Rate (Hz)')
% xticks([1:13])
% xticklabels({'-3', '-2.5', '-2', '-1.5', '-1', '-0.5', '0', ...
%            '0.5', '1', '1.5', '2', '2.5', '3'})
% xlabel('Relative to BF (oct)')
% box off
% set(gca,'TickDir','out')
% % set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,600,400,400])
% 
% set(gcf, 'Color', 'w')
% % export_fig('Frequency_tuning_good_n',  '-png', '-pdf')
% export_fig('Frequency_tuning_good',  '-png', '-pdf')

%% save data
% save('fra_tuning_re-algin_best_MGB.mat', 'tuning_best_norm')
