function rlf_plot(rlf, spls, chan)
% INPUT:
%       rlf: from rlf_implae
%       spls:  from rlf_impale
%       chan: which channel to plot the fra
data = squeeze(rlf(chan,:,:));
err = std(data,0,2)/sqrt(size(data,2)); % standard error mean
CT=cbrewer('div', 'RdYlBu', 6); % for nice color
e = errorbar(mean(data,2),err, '-s','MarkerSize',6,'MarkerFaceColor', CT(6,:), 'Color',CT(6,:), 'LineWidth',1);
hold on
xticks([1:2:length(spls)])
for i = 1:length(spls)
    labels_spls{i} = num2str(spls(i));
end
xticklabels(labels_spls(1:2:end))
ylabel('Spike Counts')
xlabel('Level (dB SPL)')
title(['Channel ', num2str(chan)])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,300])
hold off