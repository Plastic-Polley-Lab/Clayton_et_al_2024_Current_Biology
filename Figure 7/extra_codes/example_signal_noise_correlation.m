
%%
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\FRA_CF_50')
load('summary_fra_responses_v2.mat')

fra_neuron1 = fra_summary(23).fra;
fra_neuron2 = fra_summary(26).fra;
fra_neuron1_tuning = mean(fra_neuron1,3); % the 1st column is 50 dB, 2nd column is 70 dB
fra_neuron2_tuning = mean(fra_neuron2,3); % the 1st column is 50 dB, 2nd column is 70 dB
% get the best freq for each neurons
% [~, sig_cor.cor70Best(iter,1)] = max(fra_neuron1_tuning(:,2));
% [~, sig_cor.cor70Best(iter,2)] = max(fra_neuron2_tuning(:,2));


temp = corrcoef(fra_neuron1_tuning(:,2), fra_neuron2_tuning(:,2)); % 70 dB
sig_cor_cor70 = temp(1,2);
%%
figure;
hold on
plot(1:9, fra_neuron2_tuning(:,2), '-k')
scatter(1:9, fra_neuron2_tuning(:,2), 36, 'ko', 'MarkerFaceColor', 'w')
hold on
scatter(ones(1,50)*3, squeeze(fra_neuron2(3, 2, :)), 12, 'ko', 'MarkerFaceColor', 'k')

freqs = fra_summary(1).freqs;

x_ticks_indx = 1:1:length(freqs);
xticks(x_ticks_indx)
labels = {};
for i = 1:length(x_ticks_indx)
    labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
end
xticklabels(labels)
%     yticks([1:1:length(spls)])
%     for i = 1:length(spls)
%         labels_spls{i} = num2str(spls(i));
%     end
%     yticklabels(flip(labels_spls(1:1:end)))
xlabel('Frequency (kHz)')
ylabel('Spike Counts')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])
%%


figure; scatter(fra_neuron1_tuning(:,2), fra_neuron2_tuning(:,2), 36, 'ko', 'filled')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlabel('Spike Count Unit28')
ylabel('Spike Count Unit31')
temp = corrcoef(fra_neuron1_tuning(:,2), fra_neuron2_tuning(:,2)); % 70 dB
sig_cor_cor70 = temp(1,2)

figure; scatter(fra_neuron1(4,2,:), fra_neuron2(4,2,:), 36, 'ko', 'filled')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlabel('Spike Count Unit28')
ylabel('Spike Count Unit31')

for i = 1:9
    temp = corrcoef(fra_neuron1(i,2,:), fra_neuron2(i,2,:)); % 70 dB
    noise_cor70(i) = temp(1,2)
end


fra_neuron1_tuning70 = squeeze(fra_neuron1(:,2,:));
fra_neuron1_tuning70_zscore = zscore(fra_neuron1_tuning70, 0, 2);

fra_neuron2_tuning70 = squeeze(fra_neuron2(:,2,:));
fra_neuron2_tuning70_zscore = zscore(fra_neuron2_tuning70, 0, 2);

figure; scatter(fra_neuron1_tuning70(:), fra_neuron2_tuning70(:), 'ko', 'filled' )
figure; scatter(fra_neuron1_tuning70_zscore(:), fra_neuron2_tuning70_zscore(:), 36, 'ko', 'filled')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlabel('Z-Scored Spike Counts Unit28')
ylabel('Z-Scored Spike Counts Unit31')

temp = corrcoef(fra_neuron1_tuning70(:), fra_neuron2_tuning70(:)); % 70 dB
noise_cor70all = temp(1,2)

temp = corrcoef(fra_neuron1_tuning70_zscore(:), fra_neuron2_tuning70_zscore(:)); % 70 dB
noise_cor70all = temp(1,2)
%% re-calulate the signal noise correlation, data are z-scored before summary
%% calculate the singal and noise correlation
% for i = 1:length(ko_file)
%     cd(ko_file{i})
%     load('summary_fra_responses_v2.mat')
%     %signal correlation
%     sig_noise_process(fra_summary, 1)
% end
%%
figure; scatter(summaryData.ko(3).sig_cor70, summaryData.ko(3).noise_cor70all, 24, 'ro', 'filled')
ylim([-0.4, 0.8])
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlabel('Signal Correlation Coefficient')
ylabel('Noise Correlation Coefficient')
%%
figure; scatter([summaryData.wt.sig_cor70], [summaryData.wt.noise_cor70all], 12, 'ko', 'filled')
ylim([-0.4, 0.8])
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlabel('Signal Correlation Coefficient')
ylabel('Noise Correlation Coefficient')
% for i = 1:9
%     temp = corrcoef(fra_neuron1_tuning70_zscore(i,:), fra_neuron2_tuning70_zscore(i,:)); % 70 dB
%     noise_cor70all(i) = temp(1,2)
% end
% 
% temp = corrcoef(mean(fra_neuron1_tuning,2), mean(fra_neuron2_tuning,2)); % average tuning based on 50 and 70 dB
% sig_cor.corall(iter) = temp(1,2);
% iter = iter + 1;

%%
%% re-do the analysis
% if in at lease one single bin, the firing rate is significant bigger than
% the base, then that neuron is responsive
load('Summary_raw_data_files.mat')
%%

n_bins = 5; % there are 5 bins in total (e.g., 100 ms bin size, 500 ms duration)
files  = wt_file;


% stats parameter
alpha_value = 0.01;
adj_alpha = 1 - (1-alpha_value)^(1/n_bins);

% batch process
for i = 1:length(files)
    cd(files{i})
    load('summary_fra_responses_v2.mat')
    
    for j = 1:length(fra_summary)
        stats = fra_summary(j).stats;
        resp_indx = find(stats.sign == 1);
        if isempty(resp_indx)
        else
            exc_p  = stats.p(resp_indx);
            if isempty(find(exc_p < adj_alpha))
            else
                fra_summary(j).resp = 1;
            end
        end
    end
    sig_noise_process(fra_summary, 1) 
end

%%
%% summarize MGB data
clear summary
path_tosave = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Signal_Noise_Correlation';
for i = 1:length(ko_file_MGB)
    cd(ko_file_MGB{i})
    load('sumamry_sig_noise_exc_data_v2.mat')
    file_parts = strsplit(wt_file{i}, '\');
    filename = [file_parts{end-2}, '_', file_parts{end-1}];
    summary(i).filename = filename;
    summary(i).genotype = 'KO';
    summary(i).sig_cor  = sig_cor;
    summary(i).noise_cor = noise_cor;
end
cd(path_tosave)
save('summary_ko_MGB.mat', 'summary', '-v7.3')