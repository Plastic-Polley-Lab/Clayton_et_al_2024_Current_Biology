function summaryData = summarize_pupil(path)

cd(path)
genotype = {'wt', 'ko'};

for j = 1:length(genotype)
    label = genotype{j};
    load(['summary_', label, '.mat'])  %A1 analysis
    
    for i = 1:length(summary)
        pupil_area = summary(i).pupil_area;
        video_summary = summary(i).video_summary;
        
        summaryData.(label)(i).filename                    = summary(i).filename;
        summaryData.(label)(i).genotype                    = summary(i).genotype;
        summaryData.(label)(i).pupil_area                  = pupil_area;
        summaryData.(label)(i).video_summary               = video_summary;
        summaryData.(label)(i).time                        = summary(i).time;
        summaryData.(label)(i).event_time                  = summary(i).event_time;
        summaryData.(label)(i).level                       = summary(i).level;
        
        
        
        % start unwrap the data
        % step 1: get the spontaneous pupil, in other words, the pupil
        % before stimuli
        %filter the pupil_area data to smooth the cure
        pupil_area = smooth(pupil_area, 5);
        sponta = video_summary.sponta;
        summaryData.(label)(i).sponta = sponta;
        min_avg = min(pupil_area);
        max_avg = max(pupil_area);
        summaryData.(label)(i).sponta_n = (sponta - min_avg)/(max_avg - min_avg);
        summaryData.(label)(i).sponta_n_avg = mean(summaryData.(label)(i).sponta_n(:));
        
        summaryData.(label)(i).pupil_area_n = (pupil_area-min_avg)/(max_avg - min_avg);
        % get the power spectrum of the pupil fluctuation during the
        % baseline
        spont_n = summaryData.(label)(i).sponta_n;
        spont_n_center = spont_n - mean(spont_n,2);
        n_fft = size(spont_n_center, 2);
        Fs = 30;
        [pxx,f] = periodogram(spont_n_center',[], n_fft,Fs);
        summaryData.(label)(i).sponta_pxx = pxx;
        summaryData.(label)(i).sponta_f = f;
        
        
        % get the evoked pupil area:
        level = {'level35', 'level45', 'level55', 'level65', ...
            'level75', 'level85', 'level95'};
        for k = 1:length(level)
            summaryData.(label)(i).evoke_baseline{k} = (video_summary.(level{k})-min_avg)/(max_avg - min_avg);
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
for j = 1:length(genotype)
    label = genotype{j};
    for i = 1:length(summaryData.(label))
        
        
        summaryData.(label)(i).sponta_pxx_avg = mean(summaryData.(label)(i).sponta_pxx,2)';
        summaryData.(label)(i).sponta_pxx_log_avg = mean(10*log10(summaryData.(label)(i).sponta_pxx),2)';
        
        
        
        
%         summaryData.(label)(i).noise_cor_all_avg = mean(noise_cor_all, 'omitnan');
%         summaryData.(label)(i).noise_cor_avg = mean(noise_cor, 'omitnan');
%         
%         summaryData.(label)(i).positive_ratio  = positive_ratio;
%         summaryData.(label)(i).sig_cor_positive = mean(sig_cor_positive, 'omitnan');
%         summaryData.(label)(i).sig_cor_negative = mean(sig_cor_negative, 'omitnan');
%         summaryData.(label)(i).noise_cor_all_positive = mean(noise_cor_all_positive, 'omitnan');
%         summaryData.(label)(i).noise_cor_all_negative = mean(noise_cor_all_negative, 'omitnan');
%         summaryData.(label)(i).noise_cor_positive = mean(noise_cor_positive, 'omitnan');
%         summaryData.(label)(i).noise_cor_negative = mean(noise_cor_negative, 'omitnan');


        
    end
end