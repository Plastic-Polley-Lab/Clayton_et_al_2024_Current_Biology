function summaryData = spont_movement(summaryData, cutOff)

% cutOff = 0.8;

wt_indx = 10;
data = summaryData.wt(wt_indx).orofacial_n;
time = summaryData.wt(wt_indx).time;
event_time = summaryData.wt(wt_indx).event_time;

% data = summaryData.ko(wt_indx).orofacial_n;
% time = summaryData.ko(wt_indx).time;
% event_time = summaryData.ko(wt_indx).event_time;
% data_zscore = zscore(data);
deviation = std(data);
% figure; plot(data_zscore)
% hold on
% plot(1:length(data_zscore), ones(size(data_zscore)))
figure;
plot(time, data)
hold on
plot(time, cutOff * deviation*ones(size(data)))
xlim([120, 180])

for i = 1:length(summaryData.wt)
    wt_indx = i;
    for j = 1:size(summaryData.wt(wt_indx).sponta_n, 1)
        spont_smooth = smooth(summaryData.wt(wt_indx).sponta_n(j,:),5);
        deviation = std(spont_smooth);
        resp_indx = find(spont_smooth > deviation* cutOff);
        if isempty(resp_indx)
            movements(j,1) = 0;
            movements(j,2) = 0;
        else
            movements(j,1) = sum(length(resp_indx));
            rise = diff(resp_indx);
            rise_indx = find(rise>10);    % here at least 10 frame apart to set different events
            if isempty(rise_indx)
                movements(j,2) = 1;
            else
                movements(j,2) = length(rise_indx) + 1;
            end
        end
    end
    summaryData.wt(wt_indx).spont_move = movements;
    summaryData.wt(wt_indx).spont_move_avg = (mean(movements))';
    clear movements
end

clear movements
for i = 1:length(summaryData.ko)
    ko_indx = i;
    for j = 1:size(summaryData.ko(ko_indx).sponta_n, 1)
        spont_smooth = smooth(summaryData.ko(ko_indx).sponta_n(j,:),5);
        deviation = std(spont_smooth);
        resp_indx = find(spont_smooth > deviation * cutOff);
        if isempty(resp_indx)
            movements(j,1) = 0;
            movements(j,2) = 0;
        else
            movements(j,1) = sum(length(resp_indx));
            rise = diff(resp_indx);
            rise_indx = find(rise>10);    % here at least 10 frame apart to set different events
            if isempty(rise_indx)
                movements(j,2) = 1;
            else
                movements(j,2) = length(rise_indx) + 1;
            end
        end
    end
    summaryData.ko(ko_indx).spont_move = movements;
    summaryData.ko(ko_indx).spont_move_avg = (mean(movements))';

    clear movements

end