function summaryData = summarize_signal_noise(path, area)

cd(path)
genotype = {'wt', 'ko'};

for j = 1:length(genotype)
    label = genotype{j};
    switch area
        case 'ACtx'
            load(['summary_', label, '.mat'])  %A1 analysis
        case 'MGB'
            load(['summary_', label, '_MGB.mat'])  %MGB analysis
    end
    for i = 1:length(summary)
        resp      = summary(i).noise_cor.cor70all;
        predictor = summary(i).sig_cor.cor70;
        
        
        summaryData.(label)(i).filename                    = summary(i).filename;
        summaryData.(label)(i).genotype                    = summary(i).genotype;
        summaryData.(label)(i).sig_cor70                   = predictor';
        summaryData.(label)(i).noise_cor70all              = resp';
        summaryData.(label)(i).noise_cor70                 = summary(i).noise_cor.cor70;
        
        
        if length(resp)<=2
            summaryData.(label)(i).slope                       = [];
            summaryData.(label)(i).intercept                   = [];
            summaryData.(label)(i).Rsquare                     = [];
            summaryData.(label)(i).adjustRsquare               = [];
            summaryData.(label)(i).corr                        = [];
        else
            mdl       = fitlm(predictor, resp, 'linear');
            %         figure
            %         plot(mdl)
            %         xlabel('Signal Correlation Coefficient at 70 dB')
            %         ylabel('Collapse Noise Correlation Coefficient at 70 dB')
            %         ylim([-1, 1])
            %         xlim([-1, 1])
            %         title(label)
            slope = table2array(mdl.Coefficients(2,1));
            intercept = table2array(mdl.Coefficients(1,1));
            Rsquare   = mdl.Rsquared.Ordinary;
            adjustRsqure = mdl.Rsquared.Adjusted;
            corr = corrcoef(predictor, resp);
            corr = corr(1,2);
            summaryData.(label)(i).slope                       = slope;
            summaryData.(label)(i).intercept                   = intercept;
            summaryData.(label)(i).Rsquare                     = Rsquare;
            summaryData.(label)(i).adjustRsquare               = adjustRsqure;
            summaryData.(label)(i).corr                        = corr;
        end

        switch area
            case 'ACtx'
                shank = summary(i).sig_cor.shank;
                sameShank_indx = find(shank(:, 1) == shank(:,2));
                diffShank_indx = find(shank(:, 1) ~= shank(:,2));
                
                summaryData.(label)(i).sameShank_sig_cor           = summary(i).sig_cor.cor70(sameShank_indx)';
                summaryData.(label)(i).sameShank_noise_cor         = summary(i).noise_cor.cor70all(sameShank_indx)';
                
                summaryData.(label)(i).diffShank_sig_cor           = summary(i).sig_cor.cor70(diffShank_indx)';
                summaryData.(label)(i).diffShank_noise_cor         = summary(i).noise_cor.cor70all(diffShank_indx)';
            case 'MGB'
                
        end


        % get the best frequency for each neurons for each session
        n_neuron = length(summary(i).sig_cor.exc_data);
        freq     = summary(i).sig_cor.exc_data(1).freqs;
        bestFreq_indx = [summary(i).sig_cor.cor70Best(1, 1); summary(i).sig_cor.cor70Best(1:(n_neuron-1), 2)]; 
        bestFreq = freq(bestFreq_indx);
        summaryData.(label)(i).bestFreq = bestFreq;
        
        % get the pairs sharing the same best frequency
        summaryData.(label)(i).sameBest_sig_cor            = summary(i).sig_cor.sameBest_sig_corr';
        summaryData.(label)(i).sameBest_noise_cor          = summary(i).noise_cor.sameBest_noise_corr;

    end
end

%% keep updating the summaryData
for j = 1:length(genotype)
    label = genotype{j};
    for i = 1:length(summaryData.(label))
        sig_cor = summaryData.(label)(i).sig_cor70;
        noise_cor_all =  summaryData.(label)(i).noise_cor70all;
        noise_cor     =  summaryData.(label)(i).noise_cor70;
        
        sig_cor_positive = sig_cor(sig_cor > 0);
        sig_cor_negative = sig_cor(sig_cor < 0);
        
        noise_cor_all_positive = noise_cor_all(sig_cor > 0);
        noise_cor_all_negative = noise_cor_all(sig_cor < 0);
        positive_ratio         = length(sig_cor_positive)/length(sig_cor);
        
        noise_cor_positive     = noise_cor(sig_cor > 0, :);
        noise_cor_negative     = noise_cor(sig_cor < 0, :);
        
        summaryData.(label)(i).sig_cor_avg = mean(sig_cor, 'omitnan');
        summaryData.(label)(i).noise_cor_all_avg = mean(noise_cor_all, 'omitnan');
        summaryData.(label)(i).noise_cor_avg = mean(noise_cor, 'omitnan');
        
        summaryData.(label)(i).positive_ratio  = positive_ratio;
        summaryData.(label)(i).sig_cor_positive = mean(sig_cor_positive, 'omitnan');
        summaryData.(label)(i).sig_cor_negative = mean(sig_cor_negative, 'omitnan');
        summaryData.(label)(i).noise_cor_all_positive = mean(noise_cor_all_positive, 'omitnan');
        summaryData.(label)(i).noise_cor_all_negative = mean(noise_cor_all_negative, 'omitnan');
        summaryData.(label)(i).noise_cor_positive = mean(noise_cor_positive, 'omitnan');
        summaryData.(label)(i).noise_cor_negative = mean(noise_cor_negative, 'omitnan');
  
        switch area
            case 'ACtx'
                summaryData.(label)(i).sameShank_sig_cor_avg    = mean(summaryData.(label)(i).sameShank_sig_cor);
                summaryData.(label)(i).sameShank_noise_cor_avg  = mean(summaryData.(label)(i).sameShank_noise_cor);
                
                summaryData.(label)(i).diffShank_sig_cor_avg    = mean(summaryData.(label)(i).diffShank_sig_cor);
                summaryData.(label)(i).diffShank_noise_cor_avg  = mean(summaryData.(label)(i).diffShank_noise_cor);
            case 'MGB'
        end
        
        
        summaryData.(label)(i).sameBest_sig_cor_avg    = mean(summaryData.(label)(i).sameBest_sig_cor);
        summaryData.(label)(i).sameBest_noise_cor_avg  = mean(summaryData.(label)(i).sameBest_noise_cor);

        
    end
end