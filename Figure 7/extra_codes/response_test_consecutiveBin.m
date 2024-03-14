function stats = response_test_consecutiveBin(data)

spike = data(1).bin1.spike_scmatrix;
delay = data(1).bin1.stimulus.delay;
baseWindow = 50; 
consecutive_binsize = 50; % use two consecutive bins, each bin is 50 ms
n_consecutive = 2; % two consecutive bins
baseline_window = (delay-baseWindow + 1) :delay;
for i = 1:n_consecutive
    evoke_window{i} = (delay+1 + consecutive_binsize*(i-1)):(delay + consecutive_binsize *i);
end

%stats parameter
alpha_value = 0.01;
adj_alpha = 1 - (1-alpha_value)^(1/n_consecutive);

for i = 1:n_consecutive
    baseline_activity = sum(spike(:,baseline_window),2)/length(baseline_window)* 1000;
    evoke_activity = sum(spike(:,evoke_window{i}),2)/length(evoke_window{i})* 1000;
    [p(i), h(i)] = ranksum(baseline_activity, evoke_activity);
   
    if mean(evoke_activity)> mean(baseline_activity)
        sign(i) = 1;
    elseif mean(evoke_activity)< mean(baseline_activity)
        sign(i)   = -1;
    else
        sign(i) = 0;
        warning('Evoked responses are the same the spontaneous')
    end
end
significant = find(p < adj_alpha);
if isempty(significant)
    resp = 0;
    resp_sign = 0;
else
    resp = 1;
    resp_sign = sign(significant(1)); % choose the first significant response to define the sign
end

stats.p = p;
stats.sign = sign;
stats.resp = resp;
stats.resp_sign = resp_sign;