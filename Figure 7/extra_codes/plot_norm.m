% you need two inputs; one is data_wt, wt data, and the other is data_ko
function plot_norm(data_wt, data_ko)

data_wt_n = data_wt./max(data_wt,[],2); % normalize the wt data by its max value
data_ko_n = data_ko./max(data_ko,[],2); % normalize the ko data by its max value
data_wt_n_h = data_wt./data_wt(:,end);  % normalize the wt data by the highest level
data_ko_n_h = data_ko./data_ko(:,end);  % normalize the ko data by the highest level
% plot the population average data
figure;
line_errorbar_drc(data_wt, data_ko) % plot the raw data
ylabel('Firing Rate (Hz)')
xticks([1:8])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlabel('Sound Intensity(dB)')

figure;
line_errorbar_drc(data_wt_n, data_ko_n) % plot the data normalized by max value
ylabel('Normalized Firing Rate (to max FR)')
xticks([1:8])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlabel('Sound Intensity(dB)')

figure;
line_errorbar_drc(data_wt_n_h, data_ko_n_h) % plot the data normalized by highest level
ylabel('Normalized Firing Rate (to max Level)')
xticks([1:8])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlabel('Sound Intensity(dB)')

figure;
line_errorbar_drc(diff(data_wt,1,2), diff(data_ko,1,2))
ylabel('Gain per step (Hz)')
xticks([1:7])
xticklabels({'0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70'})
xlabel('Sound Intensity(dB)')

figure;
line_errorbar_drc(diff(data_wt_n, 1, 2), diff(data_ko_n, 1, 2))
ylabel('Normalized Gain per step (to max FR)')
xticks([1:7])
xticklabels({'0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70'})
xlabel('Sound Intensity(dB)')

figure;
line_errorbar_drc(diff(data_wt_n_h,1,2), diff(data_ko_n_h,1,2))
ylabel('Normalized Gain per step (to max Level)')
xticks([1:7])
xticklabels({'0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70'})
xlabel('Sound Intensity(dB)')