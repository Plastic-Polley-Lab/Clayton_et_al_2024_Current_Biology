%
function drc_parallel_response_fanofactor(neuron_num, drc_result, SNR)
for i = 1:length(SNR)
%     SNR = 20:-5:-10;
    psth = drc_result(neuron_num).psth_summary(1).(['SNR', num2str(i)]);
    fanofactor = [drc_result(neuron_num).fanofactor.fanofactor];
    corrR      = drc_result(neuron_num).corrR;
    drc_plot_withfanofactor(psth, fanofactor, 1, 1:1250, i)
end

subaxis(14,2, 2:2:14 , 'sh', 0.2, 'sv', 0.1)
plot(1:length(fanofactor), [fanofactor],'-k')
hold on
scatter(1:length(fanofactor), [fanofactor],'ok','filled')
xticks(1:length(fanofactor))
for i = 1: length(SNR)
    labels{i} = num2str(SNR(i));
end
xticklabels(labels)
xlabel('SNR (dB)')
ylabel('Fano Factor')

% ylim([0.5, 1])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)

subaxis(14,2, 16:2:28 , 'sh', 0.2, 'sv', 0.1)
plot(1:length(corrR), corrR,'-k')
hold on
scatter(1:length(corrR), corrR,'ok','filled')
xticks(1:length(fanofactor))
for i = 1: length(SNR)
    labels{i} = num2str(SNR(i));
end
xticklabels(labels)
xlabel('SNR (dB)')
ylabel('CorrCoef')

% ylim([0.5, 1])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,1000,800])




end