%%
function [gain_first, gain_sec, BestFs_counts_first, BestFs_counts_sec]= check_drifts_FRA(data, unique_BestFs)

BestFs = [data.BestF];
gains = [data.gains_value];
% unique_BestFs = unique(BestFs);
% let's divide data into two halves
half_points =floor(length(gains)/2);
gain_first = gains(1:half_points);
gain_sec   = gains(half_points+1 : end);


BestFs_counts_total =[];
BestFs_counts_first =[];
BestFs_counts_sec   =[];
for i = 1:length(unique_BestFs)
    BestFs_counts_total(i) = length(find(BestFs == unique_BestFs(i)));
    BestFs_counts_first(i) = length(find(BestFs(1:half_points) == unique_BestFs(i)));
    BestFs_counts_sec(i) = length(find(BestFs(half_points+1 :end) == unique_BestFs(i)));

end
%% let's plot the distribution of the best Frequence sampled
figure
hold on
h1 = plot(BestFs_counts_first/sum(BestFs_counts_first), '-o', 'LineWidth', 1);
h2 = plot(BestFs_counts_sec/sum(BestFs_counts_sec), '-*', 'LineWidth', 1);
legend([h1, h2], '1st half', '2nd half')
xticks([1:2:9])
% xticklabels({'4k', '5.6k', '8k', '11.3k', '16k', '22.6k', '32k', '45.3k', '64k'})
xticklabels({'4k',  '8k',  '16k',  '32k', '64k'})

xlabel('Best Frequency')
ylabel('Proportion')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
hold off

box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,400])

% try remove 64k
indx = find(BestFs == 64000);

