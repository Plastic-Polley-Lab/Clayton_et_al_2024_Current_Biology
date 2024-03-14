%% This script is used to analyze the singal-noise correlation
% load the data
wt_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\092020\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092420\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110720\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC36\121520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC40\122620\FRA_CF_50';...
    };
%%
ko_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091620\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092320\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC41\122920\FRA_CF_50';...
    };
%% load MGB data
wt_file_MGB = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC29\101520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110720\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\FRA_CF_50';...    
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\FRA_CF_50'};
%% 
ko_file_MGB = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC30\101920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111120\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111220\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121320\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\FRA_CF_50';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\FRA_CF_50';...
    };
%% preprocess the data; This part is run session by session; 
% you can also run in a batch such as: path = ko_file;
clear
path = pwd;
FRA_sigNoise_preprocess(path)
%% calculate the singal and noise correlation; this code is also run session by session
% you can also run in a batch as previous session
load('summary_fra_responses_v2.mat')
%signal correlation
sig_noise_process(fra_summary, 1)
%% summarize wt data
clear
load('Summary_raw_data_files.mat')

path_tosave = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Signal_Noise_Correlation';
for i = 1:length(wt_file)
    cd(wt_file{i})
    load('sumamry_sig_noise_exc_data_v3.mat')
    file_parts = strsplit(wt_file{i}, '\');
    filename = [file_parts{end-2}, '_', file_parts{end-1}];
    summary(i).filename = filename;
    summary(i).genotype = 'WT';
    summary(i).sig_cor  = sig_cor;
    summary(i).noise_cor = noise_cor;
end
cd(path_tosave)
save('summary_wt.mat', 'summary', '-v7.3')
%% summarize ko data
clear summary
clearvars -except ko_file
path_tosave = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Signal_Noise_Correlation';
for i = 1:length(ko_file)
    cd(ko_file{i})
    load('sumamry_sig_noise_exc_data_v3.mat')
    file_parts = strsplit(ko_file{i}, '\');
    filename = [file_parts{end-2}, '_', file_parts{end-1}];
    summary(i).filename = filename;
    summary(i).genotype = 'KO';
    summary(i).sig_cor  = sig_cor;
    summary(i).noise_cor = noise_cor;
end
cd(path_tosave)
save('summary_ko.mat', 'summary', 'ko_file', '-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% summarize wt and ko data in a easier plot format
clear
path = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Signal_Noise_Correlation';
% summaryData = summarize_signal_noise(path, 'ACtx')
summaryData = summarize_signal_noise(path, 'MGB');

%% plots the data
%% plot the signal-noise correlation
% load('summary_data_all_MGB.mat')
load('summary_data_all.mat')
signal_noise_correlation_plots(groupData, summaryData, freqs);

%% see the correlation between noise correlation and signal correlation
% ecdf_bar_plot([summaryData.wt.slope], [summaryData.ko.slope])
% ecdf_bar_plot([summaryData.wt.intercept], [summaryData.ko.intercept])
% ecdf_bar_plot([summaryData.wt.Rsquare], [summaryData.ko.Rsquare])
% ecdf_bar_plot([summaryData.wt.adjustRsquare], [summaryData.ko.adjustRsquare])
% ecdf_bar_plot([summaryData.wt.corr], [summaryData.ko.corr])

% ecdf_bar_plot([summaryData.wt.sig_cor_avg], [summaryData.ko.sig_cor_avg])
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')
% 
% ecdf_bar_plot([summaryData.wt.noise_cor_all_avg], [summaryData.ko.noise_cor_all_avg])
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')

% ecdf_bar_plot([summaryData.wt.noise_cor_avg], [summaryData.ko.noise_cor_avg])
% 
% 
% ecdf_bar_plot([summaryData.wt.positive_ratio], [summaryData.ko.positive_ratio])
% 
% ecdf_bar_plot([summaryData.wt.sig_cor_positive], [summaryData.ko.sig_cor_positive])
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')

% ecdf_bar_plot([summaryData.wt.noise_cor_all_positive], [summaryData.ko.noise_cor_all_positive])
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')

% ecdf_bar_plot([summaryData.wt.noise_cor_positive], [summaryData.ko.noise_cor_positive])

% ecdf_bar_plot([summaryData.wt.sig_cor_negative], [summaryData.ko.sig_cor_negative])
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')
% ecdf_bar_plot([summaryData.wt.noise_cor_all_negative], [summaryData.ko.noise_cor_all_negative])
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')
% 
% ecdf_bar_plot([summaryData.wt.noise_cor_negative], [summaryData.ko.noise_cor_negative])

% ecdf_bar_plot([summaryData.wt.sameBest_sig_cor_avg], [summaryData.ko.sameBest_sig_cor_avg])
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')

% ecdf_bar_plot([summaryData.wt.sameBest_noise_cor_avg], [summaryData.ko.sameBest_noise_cor_avg])
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')

% ecdf_bar_plot([summaryData.wt.sameShank_sig_cor_avg], [summaryData.ko.sameShank_sig_cor_avg])
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')
% 
% ecdf_bar_plot([summaryData.wt.sameShank_noise_cor_avg], [summaryData.ko.sameShank_noise_cor_avg])
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')



% figure
% line_sem_plot(1:9, reshape([summaryData.wt.noise_cor_positive], [], length(summaryData.wt))', '-k')
% line_sem_plot(1:9, reshape([summaryData.ko.noise_cor_positive], [], length(summaryData.ko))', '-r')
% 
% 
% figure
% line_sem_plot(1:9, reshape([summaryData.wt.noise_cor_negative], [], length(summaryData.wt))', '--k')
% line_sem_plot(1:9, reshape([summaryData.ko.noise_cor_negative], [], length(summaryData.ko))', '--r')



%% assemble all data together

%% see the population signal correlation
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    label = genotype{j};
    groupData.(label).all_sig_cor =[];
    groupData.(label).all_noise_cor =[];
    groupData.(label).all_SameBest_sig = [];
    groupData.(label).all_SameBest_noise = [];
    groupData.(label).noise_cor =[];
    
    for i = 1:length(summaryData.(label))
        temp1 = summaryData.(label)(i).sig_cor70';
        temp2 = summaryData.(label)(i).noise_cor70all';
        temp3 = summaryData.(label)(i).sameBest_sig_cor';
        temp4 = summaryData.(label)(i).sameBest_noise_cor';
        temp5 = summaryData.(label)(i).noise_cor70;
        
        groupData.(label).all_sig_cor = [groupData.(label).all_sig_cor; temp1];
        groupData.(label).all_sig_cor_avg(i) = mean(temp1);
        groupData.(label).all_noise_cor =[groupData.(label).all_noise_cor; temp2];
        groupData.(label).all_noise_cor_avg(i) = mean(temp2);
        
        groupData.(label).all_SameBest_sig =[groupData.(label).all_SameBest_sig; temp3];
        groupData.(label).all_SameBest_sig_avg(i) = mean(temp3);
        
        groupData.(label).all_SameBest_noise =[groupData.(label).all_SameBest_noise; temp4];
        groupData.(label).all_SameBest_noise_avg(i) = mean(temp4);
        
        groupData.(label).noise_cor =[groupData.(label).noise_cor; temp5];
        groupData.(label).noise_cor_avg(i,:) = mean(temp5, 1, 'omitnan');        
            
    end
end


%% lineaer fit for the group data
% predictor      = groupData.ko.all_sig_cor;
% resp = groupData.ko.all_noise_cor;
% mdl       = fitlm(predictor, resp, 'linear');
% figure
% scatter(predictor, resp, 6, 'ro', 'filled', 'MarkerFaceAlpha', 0.5)
% hold on
% plot(predictor, mdl.Fitted)
% 
% predictor      = groupData.wt.all_sig_cor;
% resp = groupData.wt.all_noise_cor;
% mdl       = fitlm(predictor, resp, 'linear');
% scatter(predictor, resp, 6, 'ko', 'filled', 'MarkerFaceAlpha', 0.5)
% hold on
% plot(predictor, mdl.Fitted, 'k')
% 
% xlabel('Signal Correlation Coefficient at 70 dB')
% ylabel('Collapse Noise Correlation Coefficient at 70 dB')
% ylim([-1, 1])
% xlim([-1, 1])
% title(label)
%% separate positive and negative signal correlation
groupData.wt.postive_indx = find(groupData.wt.all_sig_cor>0);
groupData.wt.negative_indx = find(groupData.wt.all_sig_cor<0);

groupData.ko.postive_indx = find(groupData.ko.all_sig_cor>0);
groupData.ko.negative_indx = find(groupData.ko.all_sig_cor<0);

%%
% wt_positive = groupData.wt.postive_indx;
% ko_positive = groupData.ko.postive_indx;
% 
% ecdf_bar_plot(groupData.wt.all_sig_cor(wt_positive), groupData.ko.all_sig_cor(ko_positive))
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')
% 
% ecdf_bar_plot(groupData.wt.all_noise_cor(wt_positive), groupData.ko.all_noise_cor(ko_positive))
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')

% wt_negative = groupData.wt.negative_indx;
% ko_negative = groupData.ko.negative_indx;

% ecdf_bar_plot(groupData.wt.all_sig_cor(wt_negative), groupData.ko.all_sig_cor(ko_negative))
% figure(100)
% xlabel('Signal Correlation Coefficient')
% figure(101)
% ylabel('Signal Correlation Coefficient')
% 
% ecdf_bar_plot(groupData.wt.all_noise_cor(wt_negative), groupData.ko.all_noise_cor(ko_negative))
% figure(100)
% xlabel('Noise Correlation Coefficient')
% figure(101)
% ylabel('Noise Correlation Coefficient')


% figure;
% scatter(groupData.wt.all_sig_cor(wt_positive), groupData.wt.all_noise_cor(wt_positive), 6, 'ko', 'filled', 'MarkerFaceAlpha', 0.5 )
% hold on
% scatter(groupData.ko.all_sig_cor(ko_positive), groupData.ko.all_noise_cor(ko_positive), 6, 'ro', 'filled', 'MarkerFaceAlpha', 0.5 )
%% check the noise correlaiton at different frequencies
% figure
% h_wt = line_sem_plot(1:9, groupData.wt.noise_cor(wt_positive,:),'k' )
% hold on
% h_wt_n = line_sem_plot(1:9, groupData.wt.noise_cor(wt_negative,:),'--k' )
% 
% h_ko = line_sem_plot(1:9, groupData.ko.noise_cor(ko_positive,:),'r' )
% hold on
% h_ko_n = line_sem_plot(1:9, groupData.ko.noise_cor(ko_negative,:),'--r' )
% xticks(x_ticks_indx)
% xticklabels(labels)
% xlabel('Frequency (kHz)')
% ylabel('Noise Correlation Coefficient')
% box off
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,600,400,350])
% legend([h_wt, h_wt_n, h_ko, h_ko_n], {'WT: positive', 'WT: negative', 'KO: positive', 'KO: negative'})
%% check the frequency distribution
load('summary_ko.mat')
freqs = summary(1).sig_cor.exc_data(1).freqs;
wt_bestFreq = [summaryData.wt.bestFreq];
ko_bestFreq = [summaryData.ko.bestFreq];
ko_bestFreq_count = [];
wt_bestFreq_count =[];

for i = 1:length(freqs)
    wt_bestFreq_count(i) = length(find(wt_bestFreq == freqs(i)));
    ko_bestFreq_count(i) = length(find(ko_bestFreq == freqs(i)));
end

groupData.wt.bestFreq = wt_bestFreq;
groupData.ko.bestFreq = ko_bestFreq;
groupData.wt.bestFreq_count = wt_bestFreq_count;
groupData.ko.bestFreq_count = ko_bestFreq_count;

figure;
h1 = plot(groupData.wt.bestFreq_count/length(groupData.wt.bestFreq), '-ok', 'LineWidth',1, 'MarkerFaceColor', 'w')
hold on
h2 = plot(groupData.ko.bestFreq_count/length(groupData.ko.bestFreq), '-or', 'LineWidth',1, 'MarkerFaceColor', 'w')
hold on

% re-set the x-tick labels
x_ticks_indx = 1:1:length(freqs);
xticks(x_ticks_indx)
labels = {};
for i = 1:length(x_ticks_indx)
    labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
end
xticklabels(labels)
xlabel('Best Frequency (kHz)')
ylabel('Proportion')
legend([h1, h2], {'WT', 'KO'})
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,350])
%% check the frequency distribution for each session
for i = 1:length(genotype)
    label = genotype{i};
    for j = 1:length(summaryData.(label))
        bestFreq = summaryData.(label)(j).bestFreq;
        bestFreq_ratio=[];
        for k = 1:length(freqs)
            bestFreq_ratio(k) = length(find(bestFreq == freqs(k)))/length(bestFreq);
        end
        
        groupData.(label).bestFreq_ratio_session(j,:) = bestFreq_ratio;
    end
    
    
end

% figure
% line_errorbar_drc(groupData.wt.bestFreq_ratio_session, groupData.ko.bestFreq_ratio_session)
% % re-set the x-tick labels
% x_ticks_indx = 1:1:length(freqs);
% xticks(x_ticks_indx)
% labels = {};
% for i = 1:length(x_ticks_indx)
%     labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
% end
% xticklabels(labels)
% xlabel('Best Frequency (kHz)')
% ylabel('Proportion')
% legend([h1, h2], {'WT', 'KO'})
% box off
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,600,400,350])

%% check whether the increased noise correlation was due to increased signal correlation
% figure;
% n_space = 10;
% step = 1/n_space;
% start = 0;
% temp_data_wt = groupData.wt.all_sig_cor(wt_positive);
% temp_data_noise_wt = groupData.wt.all_noise_cor(wt_positive);
% 
% for i = 1:n_space
%     indx_temp = find(temp_data_wt > start & temp_data_wt <= start+ step);
%     groupData.wt.positive_sig{i} = temp_data_wt(indx_temp);
%     groupData.wt.positive_noise{i} = temp_data_noise_wt(indx_temp);
%     start = start+ step;
% end
% 
% start = 0;
% temp_data_ko = groupData.ko.all_sig_cor(ko_positive);
% temp_data_noise_ko = groupData.ko.all_noise_cor(ko_positive);
% 
% for i = 1:n_space
%     indx_temp = find(temp_data_ko > start & temp_data_ko <= start+ step);
%     groupData.ko.positive_sig{i} = temp_data_ko(indx_temp);
%     groupData.ko.positive_noise{i} = temp_data_noise_ko(indx_temp);
%     start = start+ step;
% end

%%
% figure
% for i = 1:n_space
%     group_signal_wt(i) =  mean(groupData.wt.positive_sig{i});
%     group_signal_ko(i) =  mean(groupData.ko.positive_sig{i});
%     group_noise_wt(i) =  mean(groupData.wt.positive_noise{i});
%     group_noise_ko(i) =  mean(groupData.ko.positive_noise{i});
%     group_wt_proportion(i) = length(groupData.wt.positive_sig{i})/length(groupData.wt.postive_indx);
%     
%     group_signal_wt_sem(i) =  std(groupData.wt.positive_sig{i})/sqrt(length(groupData.wt.positive_sig{i}));
%     group_signal_ko_sem(i) =  std(groupData.ko.positive_sig{i})/sqrt(length(groupData.ko.positive_sig{i}));
%     group_noise_wt_sem(i)  =  std(groupData.wt.positive_noise{i})/sqrt(length(groupData.wt.positive_noise{i}));
%     group_noise_ko_sem(i)  =  std(groupData.ko.positive_noise{i})/sqrt(length(groupData.ko.positive_noise{i}));
%     group_ko_proportion(i) = length(groupData.ko.positive_sig{i})/length(groupData.ko.postive_indx);
% 
%     
% end
% figure
% h1 = errorbar(1:n_space,group_signal_wt, group_signal_wt_sem, '--k')
% hold on
% h2 = errorbar(1:n_space,group_signal_ko, group_signal_ko_sem, '--r')
% h3 = errorbar(1:n_space,group_noise_wt,  group_noise_wt_sem, '-k')
% h4 = errorbar(1:n_space,group_noise_ko,  group_noise_ko_sem, '-r')
% hold on
% h5 = plot(group_wt_proportion, '-k')
% h6 = plot(group_ko_proportion, '-r')
% legend([h1, h2, h3, h4, h5, h6], ....
%     {'wt-signal', 'ko-signal', 'wt-noise', 'ko-noise', 'wt-proportion', 'ko-proportion'})

%% here trying to separate the positive/negative signal correlation for each session


%%



% figure;
% h = hist3([groupData.ko.all_sig_cor, groupData.ko.all_noise_cor], [20,20] );
% rotate3d on
% xlabel('signal correlation')
% ylabel('noise correlation')
%%
% neuron1 = fra_summary(1);
% neuron2 = fra_summary(2);
% figure;
% scatter(mean(neuron1.tuning, 2), mean(neuron2.tuning, 2),'ko','filled')
% xlim([0,5])
% ylim([0,5])
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
%         figure;
%         scatter(neuron1.fra(i,j,:), neuron2.fra(i,j,:),'ko','filled')
%         temp = corrcoef(neuron1.fra(i,j,:), neuron2.fra(i,j,:));
%         R_noise(i,j) = temp(1,2);
%     end
% end


