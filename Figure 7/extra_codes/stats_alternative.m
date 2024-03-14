function fra_summary = stats_alternative(fra_summary)
% here, I define the neurons with at lease one bin significant excited as
% responsive neuorn

n_bins = 5; % there are 5 bins in total (e.g., 100 ms bin size, 500 ms duration)
% files  = wt_file;


% stats parameter
alpha_value = 0.01;
adj_alpha = 1 - (1-alpha_value)^(1/n_bins);

% batch process


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