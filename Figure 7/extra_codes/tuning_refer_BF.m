function [plots, tuning_summary] = tuning_refer_BF(tuning_summary, plots, results)
% updating the tuning_summary (re-algin the tuning curve to BF)
genotype = {'wt', 'ko'};
for i = 1:length(genotype)
    label = genotype{i};
    for j = 1:length(tuning_summary.(label))
        tuning_avg = mean(tuning_summary.(label)(j).tuning, 2);
        
        tuning_avg_center = NaN(1, 17);  % BF is centered at the 9th column
        [~, I] = max(tuning_avg);
        tuning_avg_center((10-I) : (18-I)) = tuning_avg;
        tuning_avg_center_n = tuning_avg_center./max(tuning_avg_center);
        tuning_summary.(label)(j).tuning_avg_center   = tuning_avg_center;
        tuning_summary.(label)(j).tuning_avg_center_n = tuning_avg_center_n;
        
    end
end

% updating the plots
for i = 1:length(genotype)
    label = genotype{i};
    plots.(label).tuning_avg_center = reshape([tuning_summary.(label).tuning_avg_center], 17, [])';
    plots.(label).tuning_avg_center_n = reshape([tuning_summary.(label).tuning_avg_center_n], 17, [])';
    plots.(label).dprime = [tuning_summary.(label).dprime]';
    BestF = [results.(label).BestF];
    gains = [results.(label).gains_value];
    nan_indx = find(isnan(gains));
    BestF(nan_indx) = [];
    freqs = tuning_summary.wt(1).freqs;
    for j = 1:length(freqs)
        plots.(label).BestFs_counts(j) = length(find(BestF == freqs(j)));
    end
    
end

