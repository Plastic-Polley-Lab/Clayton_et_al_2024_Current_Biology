%% Load ABR data for Figure 7 
load('\\apollo\Polley_Lab\KeChen\Data_BackUP\E_drive\Processed Data\PTCHD1-Project\ABR analysis\ABR_summary.mat')

%% plot the data
figure
data1 = summary.wt_sum.summary(:, 1:end-1);
h1 = plot(data1', 'k');
for i = 1:length(h1)
h1(i).Color = [0, 0, 0, 0.3];
end
hold on
h1 = errorbar(mean(data1,1, 'omitnan'), std(data1, 0, 1, 'omitnan')./sqrt(size(data1, 1)), '-ko');
set(h1, 'LineWidth', 1)
set(h1, 'MarkerFaceColor', 'k')
hold on

data2 = summary.ko_sum.summary(:, 1:end-1);
h2 = plot(data2', 'r');
for i = 1:length(h2)
h2(i).Color(4) = 0.3;
end

h2 = errorbar(mean(data2,1, 'omitnan'), std(data2, 0, 1, 'omitnan')./sqrt(size(data2, 1)), '-ro');
set(h2, 'LineWidth', 1)
set(h2, 'MarkerFaceColor', 'r')
legend([h1, h2], {'WT', 'KO'})
box off
xticks([0, 1:5, 6])
xlim([0,6])
xticklabels({'','8','11.3', '16', '22.6', '32',''})
ylabel('ABR Threshold (dB)')
xlabel('Sound Frequency (kHz)')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
title('ABR threshold')
set(gcf, 'Color', 'w')

