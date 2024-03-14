function t_jitter = interval_jitter_swap(t_jitter,reg_interval, Variable)
CycleLength = Variable.CycleLength
for i = 1:CycleLength-2
    stay_idx = randperm(CycleLength,i); % index for keep fixed.
    swap_idx = setdiff(1:CycleLength, stay_idx);
    for k = 1:Variable.NumCycles
        swap_interval = zeros(size(reg_interval));
        swap_interval(swap_idx) = reg_interval(swap_idx(randperm(length(swap_idx))));
        swap_interval(stay_idx) = reg_interval(stay_idx);
        t_jitter(i,(k-1)*CycleLength + 1: k*CycleLength) = swap_interval;
    end
end

    