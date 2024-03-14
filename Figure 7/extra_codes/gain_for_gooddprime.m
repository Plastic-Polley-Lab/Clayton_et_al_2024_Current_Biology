function gain_for_gooddprime(wt_good, ko_good)
temp = wt_good;
for i = 1:length(temp)
    rlf = temp(i).rlf;
    for j = 2:size(rlf,1)
        [h, p] = ttest2(rlf(1,:), rlf(j,:));
        if p< 0.05 && mean(rlf(j,:))> mean(rlf(1,:))
            wt.threshold(i) = temp(1).spls(j);
            wt.thrshold_ind(i) = j;
            wt.forgain(i,:) = mean(rlf,2)';
            [~, wt.threshold_peak(i)] = max(mean(rlf,2));
            break
        end
    end
end


% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
temp = ko_good;
for i = 1:length(temp)
    rlf = temp(i).rlf; 
    for j = 2:size(rlf,1)
        [h, p] = ttest2(rlf(1,:), rlf(j,:));
        if p< 0.05 && mean(rlf(j,:))> mean(rlf(1,:))
            ko.threshold(i) = temp(1).spls(j);
            ko.thrshold_ind(i) = j;
            ko.forgain(i,:) = mean(rlf,2)';
            [~, ko.threshold_peak(i)] = max(mean(rlf,2));
            
            break
        end
    end
end
% get rid of neurons not modulated by level
idx = find(wt.threshold ==0);
wt.threshold(idx)=[];
wt.thrshold_ind(idx)=[];
wt.forgain(idx,:)     =[];
wt.threshold_peak(idx)=[];

idx = find(ko.threshold ==0);
ko.threshold(idx)=[];
ko.thrshold_ind(idx)=[];
ko.forgain(idx,:)     =[];
ko.threshold_peak(idx)=[];


% figure;
% h1 = histogram(wt.threshold,[5:10:75],'Normalization', 'probability')
% h1.FaceColor = [0.5, 0.5, 0.5];
% hold on
% h2 =  histogram(ko.threshold,[5:10:75],'Normalization', 'probability')
% h2.FaceColor = [1 0 0];
% xlabel('Threshold Level (dB SPL)')
% ylabel('Probability')
% legend([h1, h2], 'WT', 'KO')
% box off
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,200,400,400])
a= wt.threshold';
b = ko.threshold';
b(end+1:length(a)) =  NaN;
figure
[h1, h2] = bar_plot([a, b])
h1.FaceColor = 'k';
h2.FaceColor = 'r';
h1.FaceAlpha = 0.3;
h2.FaceAlpha = 0.7;
ylim([0,50])
xticks([1,2])
xticklabels({'WT', 'KO'})
ylabel('Best Frequence Response Threshold (dB)')
% set(gcf, 'Color', 'w')
% export_fig('Response_Threshold',  '-png')
%% Let's check the gains
for i = 1: length(wt.thrshold_ind)
    start = wt.thrshold_ind(i)-1;
    stop  = wt.threshold_peak(i);
    gains = wt.forgain(i, start:stop);
    if length(gains)< 2
        warning('There is something wrong')
        wt.gains_value(i) = NaN;
    else
        wt.gains_value(i) = mean(diff(gains));
    end
end

for i = 1: length(ko.thrshold_ind)
    start = ko.thrshold_ind(i)-1;
    stop  = ko.threshold_peak(i);
    gains = ko.forgain(i, start:stop);
    if length(gains)< 2
        warning('There is something wrong')
        ko.gains_value(i) = NaN;
    else
        ko.gains_value(i) = mean(diff(gains));
    end
end
%% remove nan and save the data
wt.gains_value(isnan(wt.gains_value)) =[];
ko.gains_value(isnan(ko.gains_value)) =[];

figure;
[wt.f, wt.x] = ecdf(wt.gains_value);
hold on
[ko.f, ko.x] = ecdf(ko.gains_value);
h_wt = plot(wt.x, wt.f, '-k', 'LineWidth',1);
h_ko = plot(ko.x, ko.f, '-r', 'LineWidth',1);
legend([h_wt, h_ko], {'WT', 'KO'})
xlabel('Gain (\Delta sp/s per 10 dB step)')
ylabel('Cumulative probability')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlim([0,50])

a= wt.gains_value';
b = ko.gains_value';
b(end+1:length(a)) =  NaN;
figure
[h1, h2] = bar_plot([a, b])
h1.FaceColor = 'k';
h2.FaceColor = 'r';
h1.FaceAlpha = 0.3;
h2.FaceAlpha = 0.7;
ylim([0,20])
xticks([1,2])
xticklabels({'WT', 'KO'})
ylabel('Gain (\Delta sp/s per 10 dB step)')

[h, p] = ttest2(wt.gains_value, ko.gains_value)