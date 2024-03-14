clear
path = pwd;
region   = 'ACtx'
genotype = 'wt'
noiseburst_preprocess(path, region, genotype, 1)
%% batch process
% path = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100720\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110620\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110720\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC36\121520\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC40\122620\NoiseBursts',...
%     };

% path = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\103020\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111220\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC39\122320\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC41\122920\NoiseBursts',...
% 
%     };

% path = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\103020\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111120\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111220\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121320\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\NoiseBursts',...
%     'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\NoiseBursts',...
%     };

path = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110120\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110320\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110620\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110720\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\NoiseBursts',...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\NoiseBursts',...
    };

for i = 1:length(path)
    formatSpec = 'Processing %d out of %d files  \n';
    fprintf(formatSpec,i,length(path))
    region   = 'MGB';
    genotype = 'wt';
    noiseburst_preprocess(path{i}, region, genotype, 0)
end
%% create psth based on different bin size
clearvars -except data_ko data_wt
% load('summary_noise_wt_ACtx.mat')
% data = summarydata.wt;

load('summary_noise_ko_ACtx.mat')
data = summarydata.ko;

binsize = [2, 5, 10]; % ms
for i = 1:length(data)
    spike = data(i).clusterData.psth;
    time_ref = data(i).clusterData.stimChans{1, 1}.Delay;
    data(i).bin1.spike_scmatrix = spike.scmatrix;
    data(i).bin1.psth_avg       = mean(spike.scmatrix)/0.001;
    data(i).bin1.time           = (1:length(data(i).bin1.psth_avg)) + 0.5 -time_ref;
    data(i).bin1.spike          = spike.raster;
    data(i).bin1.stimulus       = spike.stimulus;
    for j = 1:length(binsize)
        window = 1:binsize(j):500;
        window_step = binsize(j);
        [spike_scmatrix, psth_avg, psth_avg_smooth, time] = raster_binsize(spike, window, window_step);
        data(i).(['bin', num2str(window_step)]).spike_scmatrix = spike_scmatrix;
        data(i).(['bin', num2str(window_step)]).psth_avg = psth_avg;
        data(i).(['bin', num2str(window_step)]).psth_avg_smooth = psth_avg_smooth;
        data(i).(['bin', num2str(window_step)]).time = time - time_ref;
    end
end

% normalize the firing rate of each neuron using auROC
for i = 1:length(data)
    for j = 1:length(binsize)
        window_step = binsize(j);
        spike_scmatrix = data(i).(['bin', num2str(window_step)]).spike_scmatrix;
        time = data(i).(['bin', num2str(window_step)]).time;
        baseline = find(time>=-50 & time < 0); % using the 50 ms before as baseline
        baseline_activity = spike_scmatrix(:,baseline);
        data(i).(['bin', num2str(window_step)]).PSTH_auROC = psth_auROC(spike_scmatrix, [], [], baseline_activity);
    end
    
end

% normalize the firing rate of each neuron using z-score
for i = 1:length(data)
    for j = 1:length(binsize)
        window_step = binsize(j);
        spike_scmatrix = data(i).(['bin', num2str(window_step)]).spike_scmatrix;
        psth_avg = data(i).(['bin', num2str(window_step)]).psth_avg;
        time = data(i).(['bin', num2str(window_step)]).time;
        baseline = find(time>=-100 & time < 0); % using the 50 ms before as baseline
        baseline_activity_m = mean(psth_avg(:,baseline));
        baseline_activity_std = std(psth_avg(baseline));
%         baseline_spikes = spike_scmatrix(:,baseline);
%         baseline_activity_m = mean(baseline_spikes(:));
%         baseline_activity_std = std(baseline_spikes(:));
        data(i).(['bin', num2str(window_step)]).zscore = (psth_avg - baseline_activity_m)/baseline_activity_std;
    end
    
end
%%
data_wt = data;
data_ko = data;
clearvars -except data_wt data_ko 
% save('summary_noise_ACtx_normalized.mat', 'data_ko', 'data_wt')
save('summary_noise_MGB_normalized.mat', 'data_ko', 'data_wt')

%% Try to include late onset responses; I will use consecutive bins
load('summary_noise_ACtx_normalized.mat')

for i = 1:length(data_wt)
    data_wt(i).stats = response_test_consecutiveBin(data_wt(i));
end

for i = 1:length(data_ko)
    data_ko(i).stats = response_test_consecutiveBin(data_ko(i));
end

%% replace the single bin test with consecutive bin test
for i = 1:length(data_wt)
    data_wt(i).resp = data_wt(i).stats.resp_sign * data_wt(i).stats.resp;
end

for i = 1:length(data_ko)
    data_ko(i).resp = data_ko(i).stats.resp_sign * data_ko(i).stats.resp;
end
save('summary_noise_ACtx_normalized_consecutiveBins.mat', 'data_wt', 'data_ko')
%%
% load('summary_noise_ACtx_normalized.mat')
data_ko([data_ko.resp]==0) = [];
data_wt([data_wt.resp]==0) = [];

data_ko([data_ko.refra_violation_ratio]>0.005) =[];
data_wt([data_wt.refra_violation_ratio]>0.005) =[];

data = [data_ko, data_wt];
ko_indx = 1:length(data_ko);
wt_indx = length(data_ko)+1 : length(data);
% [T, idx] = pca_cluster(data, 10, 4);
%%
figure
for sign = [1,-1] % plot exciation and inhibition
    cluster1 = find([data.resp]==sign);
    wt_cluster1 = intersect(cluster1, wt_indx);
    ko_cluster1 = intersect(cluster1, ko_indx);
    bin = 5
    psth_pop =[];
    psth_pop_t =[];
    for i = 1:length(data)
%         psth_pop(i,:) = data(i).(['bin', num2str(bin)]).PSTH_auROC;
%           psth_pop(i,:) = data(i).(['bin', num2str(bin)]).zscore;
        psth_pop(i,:) = data(i).(['bin', num2str(bin)]).psth_avg;
        psth_pop_t(i,:)= data(i).(['bin', num2str(bin)]).time;
    end
    
    switch sign
        case 1
            ko_psth_pop.exc = psth_pop(ko_cluster1, :);
            wt_psth_pop.exc = psth_pop(wt_cluster1, :);
            ko_psth_pop.time = psth_pop_t(1,:);
            wt_psth_pop.time = psth_pop_t(1,:);
            [r,c]= find(isnan(ko_psth_pop.exc));
            r = unique(r)
            ko_psth_pop.exc(r,:) = [];
            [r,c]= find(isnan(wt_psth_pop.exc));
            r = unique(r)
            wt_psth_pop.exc(r,:) = [];
            
            y1 = mean(ko_psth_pop.exc, 'omitnan');
            y2 = mean(wt_psth_pop.exc, 'omitnan');
            e1 = std(ko_psth_pop.exc, [],1, 'omitnan')/sqrt(size(ko_psth_pop.exc,1));
            e2 = std(wt_psth_pop.exc, [],1, 'omitnan')/sqrt(size(wt_psth_pop.exc,1));
        case -1
            ko_psth_pop.inh = psth_pop(ko_cluster1, :);       
            wt_psth_pop.inh = psth_pop(wt_cluster1, :);
            y1 = mean(ko_psth_pop.inh);
            y2 = mean(wt_psth_pop.inh);
            e1 = std(ko_psth_pop.inh, [],1)/sqrt(size(ko_psth_pop.inh,1));
            e2 = std(wt_psth_pop.inh, [],1)/sqrt(size(wt_psth_pop.inh,1));
    end
    x = psth_pop_t(1,:);
    color = [1, 0,0; 0, 0, 0]
    [h, p] = boundedline(x, y1, e1, x, y2, e2, 'cmap', color, 'alpha');
    h(1).LineWidth = 1;
    h(2).LineWidth = 1;
    hold on

end
legend(h, {'KO', 'WT'})
xlabel('Time (ms)')
ylabel('auROC')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])

%%
clearvars -except wt_psth_pop ko_psth_pop data bin ko_indx wt_indx
% for i = 21:30
%     figure
%     plot(ko_psth_pop(i,:))
%     
% end
t = ko_psth_pop.time;
indx = find(t>-100 & t < 300); % previous was -100 to 400
baseline = find(t>-100 & t<0);
resp_type ={'exc', 'inh'};
for i = 1:length(resp_type)
    for j = 1:size(ko_psth_pop.(resp_type{i}), 1)
        psth_temp = ko_psth_pop.(resp_type{i})(j,:);
        psth = psth_temp(indx) - mean(psth_temp(baseline));
        if j <=10 % only show the plots for the first 10 neurons
            [ACF, ACF_xcorr, lags] = xcorr_latency(psth, t(indx), bin, 1);
        else
            [ACF, ACF_xcorr, lags] = xcorr_latency(psth, t(indx), bin, 0);
        end
        ko_psth_pop.([(resp_type{i}), '_ACF'])(j,:) = ACF_xcorr';
        
    end
    ko_psth_pop.lags = lags;
end

%%
t = wt_psth_pop.time;
indx = find(t>-100 & t < 300);
baseline = find(t>-100 & t<0);
resp_type ={'exc', 'inh'};
for i = 1:length(resp_type)
    for j = 1:size(wt_psth_pop.(resp_type{i}), 1)
        psth_temp = wt_psth_pop.(resp_type{i})(j,:);
        psth = psth_temp(indx) - mean(psth_temp(baseline));
        
%         psth = wt_psth_pop.(resp_type{i})(j,indx);
        if j <=10 % only show the plots for the first 10 neurons
            [ACF, ACF_xcorr, lags] = xcorr_latency(psth, t(indx), bin, 1);
        else
            [ACF, ACF_xcorr, lags] = xcorr_latency(psth, t(indx), bin, 0);
        end
        wt_psth_pop.([(resp_type{i}), '_ACF'])(j,:) = ACF_xcorr';
        
    end
    wt_psth_pop.lags = lags;
end
lags = lags * bin;
%% summarize the lag and amplitude
[wt_psth_pop.timescale, wt_psth_pop.lags] = calculate_timescale(wt_psth_pop.exc_ACF, bin);
[ko_psth_pop.timescale, ko_psth_pop.lags] = calculate_timescale(ko_psth_pop.exc_ACF, bin);

[wt_psth_pop.inh_timescale, wt_psth_pop.inh_lags] = calculate_timescale(wt_psth_pop.inh_ACF, bin);
[ko_psth_pop.inh_timescale, ko_psth_pop.inh_lags] = calculate_timescale(ko_psth_pop.inh_ACF, bin);

wt_exc_count = length(wt_psth_pop.timescale);
wt_inh_count = length(wt_psth_pop.inh_timescale);
wt_total = wt_exc_count + wt_inh_count;
ko_exc_count = length(ko_psth_pop.timescale);
ko_inh_count = length(ko_psth_pop.inh_timescale);
ko_total = ko_exc_count + ko_inh_count;
newColor = [1, 1, 1; 0.5,0.5, 0.5]
figure
ax1 = pie([wt_exc_count, wt_inh_count])
title('WT')
ax1(1).FaceColor = [1, 1, 1];
ax1(2).FontSize = 12;
ax1(3).FaceColor = [0, 0, 0];
ax1(3).FaceAlpha = 0.4;
ax1(4).FontSize = 12;
set(gca,'fontsize',12)
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,300,300])

figure
ax2 = pie([ko_exc_count, ko_inh_count])
ax2(1).FaceColor = [1, 1, 1];
ax2(1).EdgeColor = [1, 0, 0];
ax2(2).FontSize = 12;

ax2(3).EdgeColor = [1, 0, 0];
ax2(3).FaceColor = [1, 0, 0];
ax2(3).FaceAlpha = 0.4;
ax2(4).FontSize = 12;

title('KO')
set(gca,'fontsize',12)
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,300,300])


figure
wt_ratio = [wt_exc_count, wt_inh_count]./(wt_exc_count + wt_inh_count);
ko_ratio = [ko_exc_count, ko_inh_count]./(ko_exc_count + ko_inh_count);

b = bar([wt_ratio; ko_ratio], 'stacked')
b(1).BarWidth = 0.5
b(1).FaceColor = [0.2, 0.5, 0.8];
b(2).FaceColor = [0.9, 0.5, 0.3];
b(1).LineWidth = 1;
b(2).LineWidth = 1;
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,300])
xticklabels({'WT', 'KO'})
ylabel('Proportion')

% figure
% h1 = bar(1, 1, 'w');
% hold on
% h1 = bar(1, wt_exc_count/wt_total, 'k')
% h2 = bar(2, 1, 'w')
% h2 = bar(2, ko_exc_count/ko_total, 'r')

%% convert lags into time
ko_psth_pop.lags = (ko_psth_pop.lags - 1) * bin;
ko_psth_pop.inh_lags = (ko_psth_pop.inh_lags - 1) * bin;
wt_psth_pop.lags = (wt_psth_pop.lags - 1) * bin;
wt_psth_pop.inh_lags = (wt_psth_pop.inh_lags - 1) * bin;
%% plot the results Time Scale of suppressed neurons
figure
[fig1, fig2] = ecdf_bar_plot(wt_psth_pop.inh_timescale, ko_psth_pop.inh_timescale, 'Neuron Suppressed');
set(get(fig1,'XLabel'), 'String', 'Time scale (ms)');
set(get(fig2,'YLabel'), 'String', 'Time scale (ms)');
% xlim([0,60])
% title('Neuron Suppressed')
% export_fig('Ecdf_Inh_TimeScale', '-png')

% figure(101)
% ylabel('Time scale (ms)')
% ylim([0,50])
% title('Neuron Suppressed')
% export_fig('Bar_Inh_TimeScale', '-png')

%% plot the results Lags of suppressed neurons
% ecdf_bar_plot(wt_psth_pop.inh_lags, ko_psth_pop.inh_lags)
% figure(100)
% xlabel('Lags (ms)')
% xlim([0,150])
% title('Neuron Suppressed')
% export_fig('Ecdf_Inh_lags', '-png')
% 
% figure(101)
% ylabel('Lags (ms)')
% ylim([0,100])
% title('Neuron Suppressed')
% export_fig('Bar_Inh_lags', '-png')

%% plot the results Time Scale of activated neurons
figure
[fig1, fig2] = ecdf_bar_plot(wt_psth_pop.timescale, ko_psth_pop.timescale, 'Neuron Activated');
set(get(fig1,'XLabel'), 'String', 'Time scale (ms)');
set(get(fig2,'YLabel'), 'String', 'Time scale (ms)');
% figure(100)
% xlabel('Time scale (ms)')
% xlim([0,60])
% title('Neuron Activated')
% export_fig('Ecdf_Exc_TimeScale', '-png')

% figure(101)
% ylabel('Time scale (ms)')
% ylim([0,50])
% title('Neuron Activated')
% export_fig('Bar_Exc_TimeScale', '-png')
%% plot lags for neurons activated
% ecdf_bar_plot(wt_psth_pop.lags, ko_psth_pop.lags)
% figure(100)
% xlabel('Lags (ms)')
% xlim([0,150])
% title('Neuron Activated')
% export_fig('Ecdf_Exc_lags', '-png')
% 
% figure(101)
% ylabel('Lags (ms)')
% ylim([0,100])
% title('Neuron Activated')
% export_fig('Bar_Exc_lags', '-png')
%% summarize the data for storage
summary.wt_psth_pop  = wt_psth_pop;
summary.ko_psth_pop  = ko_psth_pop;
summary.data         = data;
summary.bin          = bin;
summary.ko_indx      = ko_indx;
summary.wt_indx      = wt_indx;
summary.wt_total     = wt_total;
summary.wt_exc_count = wt_exc_count;
summary.wt_inh_count = wt_inh_count;
sumamry.ko_total     = ko_total;
summary.ko_exc_count = ko_exc_count;
summary.ko_inh_count = ko_inh_count;
summary.lags         = lags;

clearvars -except wt_psth_pop ko_psth_pop data bin ko_indx wt_indx summary lags

%% change point analysis
addpath(genpath('E:\Ke_Chen\MATLAB\ChangePointAnalysis'))
for i = 1:length(data)
    fprintf('Processing the %d/%d\n', i, length(data))
    cp_value = Timing_CP(data(i));
    data(i).cp_value = cp_value;
end
%% save the CP results, as it takes very long to run
save('summary_noise_ACtx_normalized_consecutiveBins_CP.mat', 'summary')
%% First check the response latency
for i = 1:length(data)
    if isempty(data(i).cp_value.Time_On)
        Time_On{i} = NaN;
        Time_Off{i}= NaN;
        H{i} = NaN;
    else
        Time_On{i} = data(i).cp_value.Time_On;
        Time_Off{i}= data(i).cp_value.Time_Off;
        H{i} = data(i).cp_value.H;
    end
end


for sign = [1,-1] % plot exciation and inhibition
    cluster1 = find([data.resp]==sign);
    wt_cluster1 = intersect(cluster1, wt_indx);
    ko_cluster1 = intersect(cluster1, ko_indx);
      
    switch sign
        case 1
            ko_psth_pop.exc_Time_On  = Time_On(ko_cluster1);
            ko_psth_pop.exc_Time_Off = Time_Off(ko_cluster1);
            ko_psth_pop.H            = H(ko_cluster1);
            
            wt_psth_pop.exc_Time_On  = Time_On(wt_cluster1);
            wt_psth_pop.exc_Time_Off = Time_Off(wt_cluster1);
            wt_psth_pop.H            = H(wt_cluster1);            
            
        case -1
            ko_psth_pop.inh_Time_On  = Time_On(ko_cluster1);
            ko_psth_pop.inh_Time_Off = Time_Off(ko_cluster1);
            ko_psth_pop.inh_H        = H(ko_cluster1);
            
            wt_psth_pop.inh_Time_On  = Time_On(wt_cluster1);
            wt_psth_pop.inh_Time_Off = Time_Off(wt_cluster1);
            wt_psth_pop.inh_H        = H(wt_cluster1);
    end
end

%% calculate the response onset
wt_psth_pop.exc_cp = extract_cp(wt_psth_pop.exc_Time_On, wt_psth_pop.exc_Time_Off, wt_psth_pop.H, 1);
ko_psth_pop.exc_cp = extract_cp(ko_psth_pop.exc_Time_On, ko_psth_pop.exc_Time_Off, ko_psth_pop.H, 1);

wt_psth_pop.inh_cp = extract_cp(wt_psth_pop.inh_Time_On, wt_psth_pop.inh_Time_Off, wt_psth_pop.inh_H, -1);
ko_psth_pop.inh_cp = extract_cp(ko_psth_pop.inh_Time_On, ko_psth_pop.inh_Time_Off, ko_psth_pop.inh_H, -1);

%%
clearvars -except wt_psth_pop ko_psth_pop data bin ko_indx wt_indx summary

%% plot the response onset
ecdf_bar_plot(wt_psth_pop.exc_cp.firstOnset, ko_psth_pop.exc_cp.firstOnset)
figure(100)
xlabel('Response Onset (s)')
xlim([0,0.1])
title('Neuron Activated')
export_fig('Ecdf_exc_response_onset', '-png')

figure(101)
ylabel('Response Onset (s)')
ylim([0,0.03])
title('Neuron Activated')
export_fig('Bar_Exc_response_onset', '-png')
%%
ecdf_bar_plot(wt_psth_pop.inh_cp.firstOnset, ko_psth_pop.inh_cp.firstOnset)
figure(100)
xlabel('Response Onset (s)')
xlim([0,0.1])
title('Neuron Suppressed')
export_fig('Ecdf_inh_response_onset', '-png')

figure(101)
ylabel('Response Onset (s)')
ylim([0,0.03])
title('Neuron Suppressed')
export_fig('Bar_inh_response_onset', '-png')
%% plot the first response duration
wt_exc_duration = wt_psth_pop.exc_cp.firstOffset- wt_psth_pop.exc_cp.firstOnset;
ko_exc_duration = ko_psth_pop.exc_cp.firstOffset- ko_psth_pop.exc_cp.firstOnset;
ecdf_bar_plot(wt_exc_duration, ko_exc_duration)
figure(100)
xlabel('Response Duration (s)')
xlim([0,0.25])
title('Neuron Activated')
export_fig('Ecdf_exc_response_duration', '-png')

figure(101)
ylabel('Response Duration (s)')
ylim([0,0.06])
title('Neuron Activated')
export_fig('Bar_exc_response_duration', '-png')
%%
wt_inh_duration = wt_psth_pop.inh_cp.firstOffset- wt_psth_pop.inh_cp.firstOnset;
ko_inh_duration = ko_psth_pop.inh_cp.firstOffset- ko_psth_pop.inh_cp.firstOnset;
ecdf_bar_plot(wt_inh_duration, ko_inh_duration)
figure(100)
xlabel('Response Duration (s)')
xlim([0,0.3])
title('Neuron Suppressed')
export_fig('Ecdf_inh_response_duration', '-png')

figure(101)
ylabel('Response Duration (s)')
ylim([0,0.15])
title('Neuron Suppressed')
export_fig('Bar_inh_response_duration', '-png')
%% here the duration is the continous duration including not only first resonses
ecdf_bar_plot(wt_psth_pop.exc_cp.duration, ko_psth_pop.exc_cp.duration)
%%
ecdf_bar_plot(wt_psth_pop.inh_cp.duration, ko_psth_pop.inh_cp.duration)
%% trying to separate neurons based on waveforms
latency_p2t = [data.latency_p2t];
figure; 
h1= histogram(latency_p2t, 40);
h1.FaceColor = [0.5, 0.5, 0.5];
box off
set(gcf, 'Color', 'w')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
hold off
xlabel('Trough-peak interval (ms)')
ylabel('Occurrence')
xlim([0,1.5])
hold on
plot([0.5,0.5],[0, 70], '--r')
plot([0.6,0.6],[0, 70], '--', 'color', [0.15, 0.5,0.7])
export_fig('Histogram_Trough-Peak', '-png')
%% find the neurons that are putative FS or pyramidal neurons
fs_cutoff = 0.5; 
py_cutoff = 0.6;

fs_indx = find(latency_p2t < fs_cutoff);
py_indx = find(latency_p2t > py_cutoff);

%% sumarize the data based on wt, ko
% summarize the fs neurons

% first find the indx of inhibitory and exitatory index of ko and wt
for sign = [1,-1] % plot exciation and inhibition
    cluster1 = find([data.resp]==sign);
    wt_cluster1 = intersect(cluster1, wt_indx);
    ko_cluster1 = intersect(cluster1, ko_indx);
    
    switch sign
        case 1
            ko_psth_pop.exc_indx = ko_cluster1;
            wt_psth_pop.exc_indx = wt_cluster1;

        case -1
            ko_psth_pop.inh_indx = ko_cluster1;      
            wt_psth_pop.inh_indx = wt_cluster1;
    end
end


% 2nd find the indx of fs and pyramidal neurons for ko and wt
%%
ko_psth_pop.fs_indx = intersect(fs_indx, ko_indx);
ko_psth_pop.py_indx = intersect(py_indx, ko_indx);

wt_psth_pop.fs_indx = intersect(fs_indx, wt_indx);
wt_psth_pop.py_indx = intersect(py_indx, wt_indx);

[ko_psth_pop.fs_exc_indx, ~, ko_psth_pop.ib_exc] = intersect(ko_psth_pop.fs_indx, ko_psth_pop.exc_indx);
[ko_psth_pop.fs_inh_indx, ~, ko_psth_pop.ib_inh] = intersect(ko_psth_pop.fs_indx, ko_psth_pop.inh_indx);

[wt_psth_pop.fs_exc_indx, ~, wt_psth_pop.ib_exc] = intersect(wt_psth_pop.fs_indx, wt_psth_pop.exc_indx);
[wt_psth_pop.fs_inh_indx, ~, wt_psth_pop.ib_inh] = intersect(wt_psth_pop.fs_indx, wt_psth_pop.inh_indx);
%% let's see the fs excitatory neurons
%% plot the results Lags of FS neurons
wt_exc = wt_psth_pop.ib_exc;
ko_exc = ko_psth_pop.ib_exc;
wt_inh = wt_psth_pop.ib_inh;
ko_inh = ko_psth_pop.ib_inh;
plot_timescale(wt_psth_pop, ko_psth_pop, wt_exc, ko_exc, wt_inh, ko_inh)
%% let's see the pyramidal neurons
[ko_psth_pop.pyr_exc_indx, ~, ko_psth_pop.pyr_ib_exc] = intersect(ko_psth_pop.py_indx, ko_psth_pop.exc_indx);
[ko_psth_pop.pyr_inh_indx, ~, ko_psth_pop.pyr_ib_inh] = intersect(ko_psth_pop.py_indx, ko_psth_pop.inh_indx);

[wt_psth_pop.pyr_exc_indx, ~, wt_psth_pop.pyr_ib_exc] = intersect(wt_psth_pop.py_indx, wt_psth_pop.exc_indx);
[wt_psth_pop.pyr_inh_indx, ~, wt_psth_pop.pyr_ib_inh] = intersect(wt_psth_pop.py_indx, wt_psth_pop.inh_indx);
%%
%% plot the results Lags of pyramidal neurons
wt_exc = wt_psth_pop.pyr_ib_exc;
ko_exc = ko_psth_pop.pyr_ib_exc;
wt_inh = wt_psth_pop.pyr_ib_inh;
ko_inh = ko_psth_pop.pyr_ib_inh;
plot_timescale(wt_psth_pop, ko_psth_pop, wt_exc, ko_exc, wt_inh, ko_inh)

%%
%% summarize the data for storage
summary.wt_psth_pop  = wt_psth_pop;
summary.ko_psth_pop  = ko_psth_pop;
clearvars -except wt_psth_pop ko_psth_pop data bin ko_indx wt_indx summary
save('summary_noise_ACtx_normalized_consecutiveBins_CP.mat', 'summary')
%% calculate the tau
ACF = wt_psth_pop.exc_ACF;
lags = summary.lags;
indx = wt_psth_pop.lags/bin + 1; % get the index of the ACF before the first 0 crossing
for i = 1:length(indx)
    ACF_temp = ACF(i, 1:indx(i));
    lags_temp = lags(1:indx(i));
    model = fitexp(lags_temp, ACF_temp, 0);
    wt_psth_pop.exc_model(i)= model;
    wt_psth_pop.exc_tau(i)  = model.tau_final;
    formatSpec = 'Processing %d out of %d neurons\n';
    fprintf(formatSpec,i,length(indx))
end

%%
%% calculate the tau
ACF = wt_psth_pop.inh_ACF;
lags = summary.lags;
indx = wt_psth_pop.inh_lags/bin + 1; % get the index of the ACF before the first 0 crossing
for i = 1:length(indx)
    ACF_temp = ACF(i, 1:indx(i));
    lags_temp = lags(1:indx(i));
    model = fitexp(lags_temp, ACF_temp, 0);
    wt_psth_pop.inh_model(i)= model;
    wt_psth_pop.inh_tau(i)  = model.tau_final;
    formatSpec = 'Processing %d out of %d neurons\n';
    fprintf(formatSpec,i,length(indx))
end

%% calculate the tau KO
ACF = ko_psth_pop.exc_ACF;
lags = summary.lags;
indx = ko_psth_pop.lags/bin + 1; % get the index of the ACF before the first 0 crossing
for i = 1:length(indx)
    ACF_temp = ACF(i, 1:indx(i));
    lags_temp = lags(1:indx(i));
    model = fitexp(lags_temp, ACF_temp, 0);
    ko_psth_pop.exc_model(i)= model;
    ko_psth_pop.exc_tau(i)  = model.tau_final;
    formatSpec = 'Processing %d out of %d neurons\n';
    fprintf(formatSpec,i,length(indx))
end

%%
%% calculate the tau KO
ACF = ko_psth_pop.inh_ACF;
lags = summary.lags;
indx = ko_psth_pop.inh_lags/bin + 1; % get the index of the ACF before the first 0 crossing
for i = 1:length(indx)
    ACF_temp = ACF(i, 1:indx(i));
    lags_temp = lags(1:indx(i));
    model = fitexp(lags_temp, ACF_temp, 0);
    ko_psth_pop.inh_model(i)= model;
    ko_psth_pop.inh_tau(i)  = model.tau_final;
    formatSpec = 'Processing %d out of %d neurons\n';
    fprintf(formatSpec,i,length(indx))
end
%% plot the tau

%% plot the response onset
figure
[fig1, fig2] = ecdf_bar_plot(wt_psth_pop.exc_tau, ko_psth_pop.exc_tau, 'Neuron Activated');
set(get(fig1,'XLabel'), 'String', 'Tau (ms)');
set(get(fig2,'YLabel'), 'String', 'Tau (ms)');
% figure(100)
% xlabel('Tau (ms)')
% xlim([0,100])
% title('Neuron Activated')
% export_fig('Ecdf_exc_tau', '-png')

% figure(101)
% ylabel('Tau (ms)')
% ylim([0,50])
% title('Neuron Activated')
% export_fig('Bar_Exc_tau', '-png')
%%
figure
[fig1, fig2] = ecdf_bar_plot(wt_psth_pop.inh_tau, ko_psth_pop.inh_tau, 'Neuron Suppressed');
set(get(fig1,'XLabel'), 'String', 'Tau (ms)');
set(get(fig2,'YLabel'), 'String', 'Tau (ms)');
xlim(fig1, [0, 200])

% figure(100)
% xlabel('Tau (ms)')
% xlim([0,100])
% title('Neuron Suppressed')
% export_fig('Ecdf_inh_tau', '-png')
% 
% figure(101)
% ylabel('Tau (ms)')
% ylim([0,50])
% title('Neuron Suppressed')
% export_fig('Bar_inh_tau', '-png')
%% save data for stats
filename = 'summary_noiseBurst';
xlswrite(filename,wt_psth_pop.timescale','WT_Exc_TimeScale')
xlswrite(filename,ko_psth_pop.timescale','KO_Exc_TimeScale')
xlswrite(filename,wt_psth_pop.inh_timescale','WT_Inh_TimeScale')
xlswrite(filename,ko_psth_pop.inh_timescale','KO_Inh_TimeScale')
xlswrite(filename,wt_psth_pop.exc_cp.firstOnset','WT_Exc_Onset')
xlswrite(filename,ko_psth_pop.exc_cp.firstOnset','KO_Exc_Onset')
xlswrite(filename,wt_psth_pop.exc_cp.duration','WT_Exc_duration')
xlswrite(filename,ko_psth_pop.exc_cp.duration','KO_Exc_duration')

xlswrite(filename,wt_psth_pop.inh_cp.firstOnset','WT_Inh_Onset')
xlswrite(filename,ko_psth_pop.inh_cp.firstOnset','KO_Inh_Onset')
xlswrite(filename,wt_psth_pop.inh_cp.duration','WT_Inh_duration')
xlswrite(filename,ko_psth_pop.inh_cp.duration','KO_Inh_duration')
%%
for sign = [1,-1] % plot exciation and inhibition

switch sign
    case 1
        y1 = mean(ko_psth_pop.exc, 'omitnan');
        y2 = mean(wt_psth_pop.exc, 'omitnan');
        e1 = std(ko_psth_pop.exc, [],1, 'omitnan')/sqrt(size(ko_psth_pop.exc,1));
        e2 = std(wt_psth_pop.exc, [],1, 'omitnan')/sqrt(size(wt_psth_pop.exc,1));
    case -1
        y1 = mean(ko_psth_pop.inh);
        y2 = mean(wt_psth_pop.inh);
        e1 = std(ko_psth_pop.inh, [],1)/sqrt(size(ko_psth_pop.inh,1));
        e2 = std(wt_psth_pop.inh, [],1)/sqrt(size(wt_psth_pop.inh,1));
end
x = ko_psth_pop.time(1,:);
color = [1, 0,0; 0, 0, 0]
[h, p] = boundedline(x, y1, e1, x, y2, e2, 'cmap', color, 'alpha');
h(1).LineWidth = 1;
h(2).LineWidth = 1;
hold on
end