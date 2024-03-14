
function results = getGain_threshold(data)

%set parameters
spk_violation_cutoff = 0.005;
freqs = data(1).freqs;

% get the best frequency for each neuron
for i = 1:length(data)
    tuning(i,:) = mean(data(i).tuning,2);
    spk_violation(i,:) = data(i).refra_violation_ratio;
end
violation_indx = find(spk_violation>=spk_violation_cutoff);
tuning(violation_indx,:) =[];

[~, bestF_indx] = max(tuning,[],2);
for i = 1:length(bestF_indx)
    BestF(i) = freqs(bestF_indx(i));
end

% get the d-prime for each neuron
dprime = [data.dprime];
dprime(violation_indx) =[];

% get the threshold
for i = 1:length(data)
    rlf_avg(i,:) = mean(data(i).rlf, 2); 
    rlf_min = min(rlf_avg(i,:));
    rlf_max = max(rlf_avg(i,:));
    rlf_avg_normal(i,:) = rlf_avg(i,:)/rlf_max;
end

rlf_avg(violation_indx,:) =[];
rlf_avg_normal(violation_indx,:) =[];


forgain = zeros(length(data), size(rlf_avg,2));
for i = 1:length(data)
    rlf_raw = data(i).rlf;
    for j = 2:size(rlf_raw,1)
        rlf_p_value(i) = kruskalwallis(rlf_raw', [], 'off');
        [h, p] = ttest2(rlf_raw(1,:), rlf_raw(j,:));
        if p< 0.05 && mean(rlf_raw(j,:))> mean(rlf_raw(1,:))
            threshold(i) = data(1).spls(j);
            threshold_ind(i) = j;
            forgain(i,:) = mean(rlf_raw,2)';
            [~, threshold_peak(i)] = max(mean(rlf_raw,2));
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

gains_value(violation_indx) =[];
threshold(violation_indx) =[];

% get the spontaneous rate
spont = [data.spont];
spont(violation_indx) =[];
% Let's calculate the Monotonic Index
% MI = (Rate_Max_Level - Rate_spont)/(Rate_MaxFiring - Rate_spont)
MI = (rlf_avg(:,end) -spont')./(max(rlf_avg, [], 2)- spont');


rlf_p_value(violation_indx) = [];

latency_p2t =[data.latency_p2t]; 
latency_p2t(violation_indx) = [];

data(violation_indx) =[];

for i = 1:length(gains_value)
    results(i).BestF = BestF(i);
    results(i).dprime = dprime(i);
    results(i).threshold = threshold(i);
    results(i).gains_value = gains_value(i);
    results(i).tuning      = tuning(i,:);
    results(i).rlf_avg     = rlf_avg(i,:);
    results(i).spont       = spont(i);
    results(i).MI          = MI(i);
    results(i).rlf_p_value = rlf_p_value(i);
    results(i).latency_p2t = latency_p2t(i);
    results(i).fra         = data(i).fra;
end






