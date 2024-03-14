function fra_plot(fra, freqs, spls, chan)
% INPUT:
%       fra: from fra_impale
%       freqs: from fra_impale
%       spls:  from fra_impale
%       chan: which channel to plot the fra
max_fra = max(max(squeeze(fra(chan,:,:))));
imagesc(rot90(squeeze(fra(chan,:,:))/max_fra))
colormap(hot)
x_ticks_indx = 1:5:length(freqs);
xticks(x_ticks_indx)
labels = {};
for i = 1:length(x_ticks_indx)
    labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
end
xticklabels(labels)

yticks([1:2:length(spls)])
for i = 1:length(spls)
    labels_spls{i} = num2str(spls(i));
end
yticklabels(flip(labels_spls(1:2:end)))
xlabel('Frequency (kHz)')
ylabel('Level (dB SPL)')
title(['Channel ', num2str(chan)])
box off
set(gca,'TickDir','out')
% set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,600,400,300])