function summaryData = summary_paired_pulse_analysis(path, region)
% load('drc_results_KO.mat')
cd(path)
genotype = {'wt', 'ko'}
for k = 1:length(genotype)
    switch region
        case 'ACtx'
            load(['summary_',  genotype{k}, '.mat'])
        case 'MGB'
            load(['summary_',  genotype{k}, '_MGB.mat'])
    end
    resp_indx = find([summary_data.resp] ==1);
    temp = summary_data(resp_indx);
    for i = 1: length(temp)
        temp(i).performance_avg = mean(temp(i).performance', 1);
    end
    summaryData.(genotype{k}).data = temp; % save all neurons with excitatory responses

    iter = [];
    for i = 1:length(temp)
        if isempty(temp(i).resp_decoding)
        else
            iter = [iter, i] ;
        end
    end
    summaryData.(genotype{k}).data_stats = temp(iter); % save all neurons with significant decoding at 256 ms
end
summaryData.region   = region;

summaryData.IPI = [0, 2.^([0:0.5:9])]; % inter pulse intervals

% load('summary_data_all.mat')
% load('paired_pulse_results_wt_spont.mat')
% load('paired_pulse_results_ko_spont.mat')
% 
% for i = 1:length(summary_data)
%     summaryData.wt.data(i).spont = summary_data(i);
% end
    


% save('drc_results_summary005', 'summary', 'threshold')