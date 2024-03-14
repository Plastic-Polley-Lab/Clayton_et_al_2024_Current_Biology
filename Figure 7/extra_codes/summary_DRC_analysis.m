function summaryData = summary_DRC_analysis(path, threshold, region)
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
    
    iter = 1;
    CrossCoef =[];
    
    for i = 1:length(temp)
        if ~isempty(temp(i).corrR)
            %         if temp(i).resp ==1
            if temp(i).corrR(1) >= threshold
                CrossCoef(iter, :) = temp(i).corrR;
                iter = iter + 1;
            end
        end
    end
    keep = [temp.keep];
    resp_ratio = size(CrossCoef, 1)/length(keep);
    summaryData.(genotype{k}).CrossCoef = CrossCoef;
    summaryData.(genotype{k}).CrossCoef_n = CrossCoef./CrossCoef(:,1);
    summaryData.(genotype{k}).keep = keep;
    summaryData.(genotype{k}).resp_ratio = resp_ratio;
    summaryData.(genotype{k}).threshold = threshold;
    summaryData.(genotype{k}).region   = region;
end
% save('drc_results_summary005', 'summary', 'threshold')