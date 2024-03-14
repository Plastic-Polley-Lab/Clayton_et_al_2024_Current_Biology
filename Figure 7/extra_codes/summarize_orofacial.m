function summaryData = summarize_orofacial(path)
cd(path)
genotype = {'wt', 'ko'};
region = 'posterior';
% region = 'anterior';

switch region
    case 'posterior'
        label_save = 'Oro_post';
        label_exist = 'Mtrace_posterior';
    case 'anterior'
        label_save = 'Oro_ant';
        label_exist = 'Mtrace_anterior';
        
end

for j = 1:length(genotype)
    label = genotype{j};
    load(['summary_', label, '.mat'])  %A1 analysis
%     load(['summary_', label, '_long.mat'])  %timecourse analysis

    for i = 1:length(summary)
        orofacial = summary(i).orofacial.(label_exist);
        video_summary = summary(i).video_summary.(label_exist);
        
        summaryData.(label)(i).filename                    = summary(i).filename;
        summaryData.(label)(i).genotype                    = summary(i).genotype;
        summaryData.(label)(i).(label_save)                = orofacial;
        summaryData.(label)(i).video_summary               = video_summary;
        summaryData.(label)(i).time                        = summary(i).time;
        summaryData.(label)(i).event_time                  = summary(i).event_time;
        summaryData.(label)(i).level                       = summary(i).level;
        % start unwrap the data
        % step 1: get the spontaneous pupil, in other words, the pupil
        % before stimuli
        outliers = find(orofacial>30 | orofacial ==0); % this usually occurs at the first frame;
        if ~isempty(outliers)
            for kkk = 1:length(outliers)
                fprintf('There are outfliers at frame %d with Intensity %f \n', outliers(kkk), orofacial(outliers(kkk)))
            end
            orofacial(outliers) = NaN;
        end
        
        orofacial_new=inpaint_nans(orofacial, 5);

        figure
        orofacial = smooth(orofacial_new, 5);
        plot(orofacial)
        min_avg = min(orofacial);
        max_avg = max(orofacial);
        
        
        edges = 0:0.05:30;
        h = histogram(orofacial, edges);

        [~, bin_indx] = max(h.BinCounts);
        min_base = h.BinEdges(bin_indx);
        sponta = summary(i).video_summary.sponta_post;
        sponta(sponta>30 | sponta == 0) = NaN;
        summaryData.(label)(i).sponta = sponta;

%         summaryData.(label)(i).sponta_n = (sponta - min_avg)/(max_avg - min_avg);
        summaryData.(label)(i).sponta_n = (sponta - min_base)/ min_base;

        summaryData.(label)(i).sponta_n_avg = mean(summaryData.(label)(i).sponta_n(:), 'omitnan');
%         summaryData.(label)(i).orofacial_n = (orofacial-min_avg)/(max_avg - min_avg);
        summaryData.(label)(i).orofacial_n = (orofacial-min_base)/min_base;

        % get the evoked pupil area:
        level = {'level35', 'level45', 'level55', 'level65', ...
            'level75', 'level85', 'level95'};
        for k = 1:length(level)
            summaryData.(label)(i).evoke_baseline{k} = (video_summary.(level{k})-min_base)/min_base;

%             summaryData.(label)(i).evoke_baseline{k} = (video_summary.(level{k})-min_avg)/(max_avg - min_avg);
            summaryData.(label)(i).evoke_n{k}        = video_summary.([level{k},'_norm']);
            summaryData.(label)(i).evoke_avg(k,:)    = video_summary.([level{k}, '_avg']);
        end
        
        % re-construct the timestamps for the sound events
        step = 0.033;
        t= -0.033 * 60: step: 0;
        t_post = 0.033:step:300*step;
        t = [t, t_post];
        summaryData.(label)(i).t = t;
        
    end
end

%% keep updating the summaryData
% for j = 1:length(genotype)
%     label = genotype{j};
%     for i = 1:length(summaryData.(label))
%         
%         
%         summaryData.(label)(i).sponta_pxx_avg = mean(summaryData.(label)(i).sponta_pxx,2)';
%         summaryData.(label)(i).sponta_pxx_log_avg = mean(10*log10(summaryData.(label)(i).sponta_pxx),2)';
%         
%         
%         
%         
% %         summaryData.(label)(i).noise_cor_all_avg = mean(noise_cor_all, 'omitnan');
% %         summaryData.(label)(i).noise_cor_avg = mean(noise_cor, 'omitnan');
% %         
% %         summaryData.(label)(i).positive_ratio  = positive_ratio;
% %         summaryData.(label)(i).sig_cor_positive = mean(sig_cor_positive, 'omitnan');
% %         summaryData.(label)(i).sig_cor_negative = mean(sig_cor_negative, 'omitnan');
% %         summaryData.(label)(i).noise_cor_all_positive = mean(noise_cor_all_positive, 'omitnan');
% %         summaryData.(label)(i).noise_cor_all_negative = mean(noise_cor_all_negative, 'omitnan');
% %         summaryData.(label)(i).noise_cor_positive = mean(noise_cor_positive, 'omitnan');
% %         summaryData.(label)(i).noise_cor_negative = mean(noise_cor_negative, 'omitnan');
% 
% 
%         
%     end
% end