% Let's analyze the Rate Level function for single units
%% Load dataset
clear
path = pwd;
RLF_analysis_preprocess(path, 'ACtx')
%% summarize WT the data
wt_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC36\121520\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC40\122620\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC43\021621\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC43\021721\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC44\031521\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC44\031621\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC46\032221\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC46\032321\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC47\032321\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC47\032421\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC48\040821\RateLevel'};
data_wt = [];
for i = 1:length(wt_file)
    cd(wt_file{i})
    load('summary_rlf_responses_v2.mat')
    badUnits = [];
    for j = 1:length(spikedata)
        if isempty(spikedata(j).resp)
            badUnits = [badUnits, j];
        end
    end
    spikedata(badUnits) = [];
    data_wt = [data_wt, spikedata];
end

%% Sumarize the KO data
ko_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC39\122320\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC41\122920\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC45\031521\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC45\031621\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC49\040721\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC49\040821\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC50\040721\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC50\040821\RateLevel'};
data_ko = [];
for i = 1:length(ko_file)
    cd(ko_file{i})
    load('summary_rlf_responses_v2.mat')
    badUnits = [];
    for j = 1:length(spikedata)
        if isempty(spikedata(j).resp)
            badUnits = [badUnits, j];
        end
    end
    spikedata(badUnits) = [];
    data_ko = [data_ko, spikedata];
end
%%
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\RLF_analysis')
save('summaryData_raw.mat', 'data_wt', 'data_ko')
%% Let's analyze excitatory neurons
load('summaryData_raw.mat')
clearvars -except data_wt data_ko
% load('summary_noise_ACtx_normalized.mat')
data_ko([data_ko.resp]==0) = [];
data_wt([data_wt.resp]==0) = [];

data_ko([data_ko.refra_violation_ratio]>0.005) =[];
data_wt([data_wt.refra_violation_ratio]>0.005) =[];

data = [data_ko, data_wt];
ko_indx = 1:length(data_ko);
wt_indx = length(data_ko)+1 : length(data);
%% summarize the data; let's first only look at excitatory responses
sign = 1; % extract exciation 
cluster1 = find([data.resp]==sign);
wt_cluster1 = intersect(cluster1, wt_indx);
ko_cluster1 = intersect(cluster1, ko_indx);
rlf =[];
results.wt.exc = gain_threshold_rlf(data(wt_cluster1));
results.ko.exc = gain_threshold_rlf(data(ko_cluster1));

%%
%% plot the Monotonic index
figure
[fig1, fig2] = ecdf_bar_plot([results.wt.exc.MI], [results.ko.exc.MI], 'Monotonic Index');
set(get(fig1,'XLabel'), 'String', 'Monotonic Index');
set(get(fig2,'YLabel'), 'String', 'Monotonic Index');
% export_fig('Monotonic Index', '-png')
%% plot the gains
figure
[fig1, fig2] = ecdf_bar_plot([results.wt.exc.gains_value], [results.ko.exc.gains_value], 'Gains');
set(get(fig1,'XLabel'), 'String', 'Gains (\Delta sp/s per 5 dB step)');
set(get(fig2,'YLabel'), 'String', 'Gains (\Delta sp/s per 5 dB step)');
% xlim(fig1, [0,20])
xlim(fig1, [-10,0])
% ylim(fig2,[0,15])
ylim(fig2,[-4, 0])

% export_fig('bar_gains',  '-png')

%% plot the monotonic index with the gains
figure
h1 = scatter([results.wt.exc.MI], [results.wt.exc.gains_value], 'ko', 'filled');
h1.MarkerFaceAlpha = 0.3;
hold on
h2 = scatter([results.ko.exc.MI], [results.ko.exc.gains_value], 'ro', 'filled');
h2.MarkerFaceAlpha = 0.3
xlabel('Monotonic Index')
ylabel('Gains (\Delta sp/s per 5 dB step)')
box off
set(gcf, 'Color', 'w')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
hold off
%% organize the data for plots
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    for i = 1:length(results.(genotype{j}).exc)
%         plots.(genotype{j}).exc.rlf(i,:) = results.(genotype{j}).exc(i).rlf_avg;
        plots.(genotype{j}).exc.rlf(i,:) = results.(genotype{j}).exc(i).rlf_avg - ...
            results.(genotype{j}).exc(i).spont;

    end
    plots.(genotype{j}).exc.gains       = [results.(genotype{j}).exc.gains_value];
    plots.(genotype{j}).exc.threshold   = [results.(genotype{j}).exc.threshold];
    plots.(genotype{j}).exc.latency_p2t = [results.(genotype{j}).exc.latency_p2t];
    plots.(genotype{j}).exc.rlf_p_value = [results.(genotype{j}).exc.rlf_p_value];
    plots.(genotype{j}).exc.MI = [results.(genotype{j}).exc.MI];
    nan_indx = find(isnan(plots.(genotype{j}).exc.threshold));
    
    plots.(genotype{j}).exc.gains(nan_indx)=[];
    plots.(genotype{j}).exc.threshold(nan_indx) =[];
    plots.(genotype{j}).exc.rlf(nan_indx,:) = [];
    plots.(genotype{j}).exc.latency_p2t(nan_indx) = [];
    plots.(genotype{j}).exc.rlf_p_value(nan_indx) = [];
    plots.(genotype{j}).exc.MI(nan_indx) =[];
end
%%
figure;

line_errorbar_drc(plots.wt.exc.rlf, plots.ko.exc.rlf)
xticks([1:2:17])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80'})
xlabel('Sound Intensity(dB)')
ylabel('Firing Rate (Hz)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
xlim([0, 18])
% export_fig('Rate_Level_Function',  '-png', '-pdf')

%% Let's see the threshold
figure
[fig1, fig2] = ecdf_bar_plot(plots.wt.exc.threshold', plots.ko.exc.threshold', 'threshold');
set(get(fig1,'XLabel'), 'String', 'Threshold (dB)');
set(get(fig2,'YLabel'), 'String', 'Threshold (dB)');
xlim(fig1, [0,80])
ylim(fig2,[0,40])

set(gcf, 'Color', 'w')
% export_fig('Response_Threshold',  '-png')
%% Let's see the MI
figure
[fig1, fig2] = ecdf_bar_plot(plots.wt.exc.MI', plots.ko.exc.MI', 'MI');
set(get(fig1,'XLabel'), 'String', 'Threshold (dB)');
set(get(fig2,'YLabel'), 'String', 'Threshold (dB)');
xlim(fig1, [0,80])
ylim(fig2,[0,40])

set(gcf, 'Color', 'w')
% export_fig('Response_Threshold',  '-png')
%% save the summary data
save('summaryData_plots_MGB.mat','plots', 'results')
save('summaryData_plots_ACtx.mat','plots', 'results')
save('summaryData_plots_ACtx_inh.mat','plots', 'results')
save('summaryData_plots_MGB_inh.mat','plots', 'results')

%% Check the gain with certain range: from 35 to 55 dB
level = 0:5:80;
range_indx = find(level>=35 & level <=55);
genotype = {'wt', 'ko'};

for i = 1:length(genotype)
    temp_rlf = plots.(genotype{i}).exc.rlf;
    plots.(genotype{i}).exc.range_rlf = temp_rlf(:, range_indx);
    plots.(genotype{i}).exc.range_gain      = plots.(genotype{i}).exc.range_rlf(:, end) - plots.(genotype{i}).exc.range_rlf(:, 1);
end

figure
[fig1, fig2] = ecdf_bar_plot(plots.wt.exc.range_gain, plots.ko.exc.range_gain, 'ranged gain');
set(get(fig1,'XLabel'), 'String', 'Gain (\Delta sp/s from 35 to 55 dB)');
set(get(fig2,'YLabel'), 'String', 'Gain (\Delta sp/s from 35 to 55 dB)');


figure
[fig1, fig2] = ecdf_bar_plot(abs(plots.wt.exc.range_gain), abs(plots.ko.exc.range_gain), 'ranged gain');
set(get(fig1,'XLabel'), 'String', 'Gain (absolute \Delta sp/s from 35 to 55 dB)');
set(get(fig2,'YLabel'), 'String', 'Gain (absolute \Delta sp/s from 35 to 55 dB)');

%% cluster the rlf
%% try cluster the rate-level function
wt_indx = 1:size(plots.wt.exc.rlf, 1);
ko_indx = (wt_indx(end) + 1): (wt_indx(end) + size(plots.ko.exc.rlf,1));

alldata = [plots.wt.exc.rlf; plots.ko.exc.rlf];
alldata_n = alldata./max(alldata, [], 2);
clusterN = 3; % A1: 4
[T, idx] = pca_cluster(alldata_n, [], clusterN);

wt_cluster = T(wt_indx);
ko_cluster = T(ko_indx);
% plot the cluster proportion
for i = 1:clusterN
    wt_temp = find(wt_cluster == i);
    ko_temp = find(ko_cluster == i);
    plots.wt.exc.clusterRatio(i)  = length(wt_temp);
    plots.ko.exc.clusterRatio(i)  = length(ko_temp);
    plots.wt.exc.cluster_gains{i} = plots.wt.exc.gains(wt_temp);
    plots.ko.exc.cluster_gains{i} = plots.ko.exc.gains(ko_temp);
    plots.wt.exc.cluster_rlf{i}   = plots.wt.exc.rlf(wt_temp,:);
    plots.ko.exc.cluster_rlf{i}   = plots.ko.exc.rlf(ko_temp,:);    
end

figure;
h1 = plot(plots.wt.exc.clusterRatio/sum(plots.wt.exc.clusterRatio), '-o', 'LineWidth', 1);
hold on
h2 = plot(plots.ko.exc.clusterRatio/sum(plots.ko.exc.clusterRatio), '-o', 'LineWidth', 1);

h1.Color = 'k';
h1.MarkerFaceColor = [0.5, 0.5, 0.5];
h2.Color = 'r';
h2.MarkerFaceColor = 'r';

legend([h1, h2], 'WT', 'KO')
xticks([1:1:6])
% xticklabels({'4k', '5.6k', '8k', '11.3k', '16k', '22.6k', '32k', '45.3k', '64k'})
xticklabels({'1',  '2',  '3',  '4', '5', '6'})

xlabel('Cluster #')
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
    line_errorbar_drc(plots.wt.exc.cluster_rlf{i}, plots.ko.exc.cluster_rlf{i})
    xticks([1:2:17])
    xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80'})
    xlabel('Sound Intensity(dB)')
    ylabel('Firing Rate (Hz)')
    set(gcf,'position',[100,200,400,400])
    set(gcf, 'Color', 'w')
    title(['Cluster ', num2str(i)])
    xlim([0, 18])
end

%% summarize the data; let's first only look at excitatory responses
sign = -1; % extract exciation 
cluster1 = find([data.resp]==sign);
wt_cluster1 = intersect(cluster1, wt_indx);
ko_cluster1 = intersect(cluster1, ko_indx);
rlf =[];
results.wt.exc = gain_threshold_rlf(data(wt_cluster1));
results.ko.exc = gain_threshold_rlf(data(ko_cluster1));

%% Fisher Information analysis
clear
% load('summaryData_raw.mat')
load('summaryData_raw_MGB.mat')

% load('summary_noise_ACtx_normalized.mat')
data_ko([data_ko.resp]==0) = [];
data_wt([data_wt.resp]==0) = [];

data_ko([data_ko.refra_violation_ratio]>0.005) =[];
data_wt([data_wt.refra_violation_ratio]>0.005) =[];

data = [data_ko, data_wt];
ko_indx = 1:length(data_ko);
wt_indx = length(data_ko)+1 : length(data);

for i = 1:length(data)
    fprintf('Process Neuron # %d\n', i)
    rlf = data(i).rlf;
    fisher_info(i,:) = fisher_information_analysis(rlf);
end

figure;
line_errorbar_drc(fisher_info(wt_indx, :), fisher_info(ko_indx, :))
xticks([1:2:13])
xticklabels({'10', '20', '30', '40', '50', '60', '70'})
xlabel('Sound Intensity(dB)')
ylabel('Fisher Information (1/dB^2)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
xlim([0, 14])
%% Check the neuro-adaptation
clearvars
load('summaryData_raw.mat')
% load('summary_noise_ACtx_normalized.mat')
data_ko([data_ko.resp]==0) = [];
data_wt([data_wt.resp]==0) = [];

data_ko([data_ko.refra_violation_ratio]>0.005) =[];
data_wt([data_wt.refra_violation_ratio]>0.005) =[];

data = [data_ko, data_wt];
ko_indx = 1:length(data_ko);
wt_indx = length(data_ko)+1 : length(data);
level = 65
for i = 1:length(data)
    fprintf('Process Neuron # %d\n', i)
    rlf = data(i).rlf;
    [adaptation{i}, boot_adaptation{i}] = rlf_adaptation(rlf);
end

%% let's find the p-value of each neuron
intensity = 0:5:80;
test = 45;
indx = find(intensity == test);
for i = 1:length(data)
    adaptation_info(i).raw = adaptation{i}(indx,:);
    adaptation_info(i).AI  = adaptation{i}(indx, end) - adaptation{i}(indx, 1);
    pseudo_adaptation = boot_adaptation{i};
    for j = 1:length(pseudo_adaptation)
        adaptation_info(i).pseudo_AI(j) = pseudo_adaptation{j}(indx, end) - pseudo_adaptation{j}(indx, 1);  
    end
    adaptation_info(i).p_value = length(find([adaptation_info(i).pseudo_AI]>= abs(adaptation_info(i).AI)))/length([adaptation_info(i).pseudo_AI]);
end

%
resp_indx = find([adaptation_info.p_value] < 0.05);
adapt_indx = find([adaptation_info.p_value] < 0.05 & [adaptation_info.AI] <0);
facil_indx = find([adaptation_info.p_value] < 0.05 & [adaptation_info.AI] >0);

wt_adapt_indx = intersect(adapt_indx, wt_indx);
wt_facil_indx = intersect(facil_indx, wt_indx);
wt_non        = intersect(wt_indx, find([adaptation_info.p_value] >= 0.05));

ko_adapt_indx = intersect(adapt_indx, ko_indx);
ko_facil_indx = intersect(facil_indx, ko_indx);
ko_non        = intersect(ko_indx, find([adaptation_info.p_value] >= 0.05));

raw_data = reshape([adaptation_info.raw], 5, [])';

figure;
line_errorbar_drc(raw_data(wt_adapt_indx, :)./ raw_data(wt_adapt_indx, 1), raw_data(ko_adapt_indx, :)./ raw_data(ko_adapt_indx, 1));

figure;
line_errorbar_drc(raw_data(wt_facil_indx, :)./ raw_data(wt_facil_indx, end), raw_data(ko_facil_indx, :)./ raw_data(ko_facil_indx, end));

% figure;
% line_errorbar_drc(raw_data(wt_non, :)./ raw_data(wt_non, 1), raw_data(ko_non, :)./ raw_data(ko_non, 1));

figure
wt_ratio = [length(wt_adapt_indx), length(wt_facil_indx), length(wt_non) ]./length(wt_indx);
ko_ratio = [length(ko_adapt_indx), length(ko_facil_indx), length(ko_non)]./length(ko_indx);

b = bar([wt_ratio; ko_ratio], 'stacked')
b(1).BarWidth = 0.5
b(1).FaceColor = [0.2, 0.5, 0.8];
b(2).FaceColor = [0.9, 0.5, 0.3];
b(3).FaceColor = 'w';
b(1).LineWidth = 1;
b(2).LineWidth = 1;
b(3).LineWidth = 1;
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,300])
xticklabels({'WT', 'KO'})
ylabel('Proportion')
% spont = [data.spont] * 0.1;
% figure;
% line_errorbar_drc(adaptation(wt_indx, :) - spont(wt_indx)', adaptation(ko_indx, :) - spont(ko_indx)');
% 
% adaptation_wt = adaptation(wt_indx, :);
% spont_wt = spont(wt_indx)';
% adaptation_ko = adaptation(ko_indx, :);
% spont_ko = spont(ko_indx)';
% 
% wt_facilitate = find(adaptation_wt(:, end) > adaptation_wt(:, 1));
% wt_habituate  = find(adaptation_wt(:, end) < adaptation_wt(:, 1));
% 
% ko_facilitate = find(adaptation_ko(:, end) > adaptation_ko(:, 1));
% ko_habituate  = find(adaptation_ko(:, end) < adaptation_ko(:, 1));
% figure;
% line_errorbar_drc(adaptation_wt(wt_facilitate, :)./adaptation_wt(wt_facilitate, end), adaptation_ko(ko_facilitate, :)./adaptation_ko(ko_facilitate, end));
% 
% figure;
% line_errorbar_drc(adaptation_wt(wt_habituate, :)./adaptation_wt(wt_habituate, 1), adaptation_ko(ko_habituate, :)./ adaptation_ko(ko_habituate, 1));


% let's put all data 



% %% group the trials based on the sound level
% innerIndexes = spikedata(neuron_num).clusterData.stimData.innerIndexes;
% inner_var = spikedata(neuron_num).clusterData.stimData.inner_variables;
%
% %%
% % sort the trials based on different level
% for i = 1:max(innerIndexes) % each inner index represents a noise level
%     indx = find(innerIndexes==i);
%     psth = spikedata(neuron_num).clusterData.psth;
%     psth1.(['Level',num2str(i)]).raster = psth.raster(indx );
%     psth1.(['Level',num2str(i)]).scmatrix = psth.scmatrix(indx ,:);
%     psth1.(['Level',num2str(i)]).stimulus = psth.stimulus;
% end
%
% %% plot the rastes and psth of the responses
% % figure
% for i = 1:2:17
%     figure
%     psth_plot(psth1.(['Level',num2str(i)]),1, 1:10:750, 1)
% end
% %% calculate the cross-correlation