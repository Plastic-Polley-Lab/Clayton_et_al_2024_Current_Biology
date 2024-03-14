function results = gain_threshold_rlf(data)

n_level = size(data(1).rlf, 1);
n_rep   = size(data(1).rlf, 2);

forgain = zeros(length(data), n_level);
for i = 1:length(data)
    rlf_raw = data(i).rlf;
    rlf_avg(i,:) =  mean(rlf_raw,2)';
    for j = 2:size(rlf_raw,1)
        rlf_p_value(i) = kruskalwallis(rlf_raw', [], 'off');
        [h, p] = ttest2(rlf_raw(1,:), rlf_raw(j,:));
%         if p< 0.05 && mean(rlf_raw(j,:))> mean(rlf_raw(1,:))
        if p< 0.05 && mean(rlf_raw(j,:)) < mean(rlf_raw(1,:)) % also get the threshold for inhibitory responses

            threshold(i) = data(1).spls(j);
            threshold_ind(i) = j;
            forgain(i,:) = mean(rlf_raw,2)';
%             [~, threshold_peak(i)] = max(mean(rlf_raw,2));
            [~, threshold_peak(i)] = min(mean(rlf_raw,2));
            

            break
        else
            threshold(i) = NaN;
            threshold_ind(i) = NaN;
            forgain(i,:) = NaN;
            threshold_peak(i) = NaN;
        end
    end
end



% Let's check the gains
for i = 1: length(threshold_ind)
    if isnan(threshold_ind(i))
        gains_value(i) = NaN;
    else
        start = threshold_ind(i)-1;
        stop  = threshold_peak(i);
        gains = forgain(i, start:stop);
        if length(gains)< 2
            warning('There is something wrong')
            gains_value(i) = NaN;
        else
            gains_value(i) = mean(diff(gains));
        end
    end
end

% Let's calculate the Monotonic Index
% MI = (Rate_Max_Level - Rate_spont)/(Rate_MaxFiring - Rate_spont)
latency_p2t = [data.latency_p2t];
spont = [data.spont];
% MI = (rlf_avg(:,end) -spont')./(max(rlf_avg, [], 2)- spont'); %excitation
MI = (rlf_avg(:,end) -spont')./(min(rlf_avg, [], 2)- spont'); %inhibition

% store the data
for i = 1:length(gains_value)

    results(i).threshold   = threshold(i);
    results(i).gains_value = gains_value(i);
    results(i).rlf_avg     = rlf_avg(i,:);
    results(i).spont       = spont(i);
    results(i).MI          = MI(i);
    results(i).rlf_p_value = rlf_p_value(i);
    results(i).gains       = gains;
    results(i).latency_p2t = latency_p2t(i);
end