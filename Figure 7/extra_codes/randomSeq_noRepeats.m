function reg_interval = randomSeq_noRepeats(intervals, Variable)
if Variable.Seed < 0,
    rng('shuffle');
else
    rng(Variable.Seed + Variable.CycleLength*Variable.NumCycles);
end

cycle_duration = mean(intervals)*Variable.CycleLength;
for i = 1:Variable.CycleLength
    if (i == 1)
        reg_interval(1,i) = intervals(randi(length(intervals)));
    else
        sum_rem = cycle_duration - sum(reg_interval(1,1:i-1));
        if (i == Variable.CycleLength)
            reg_interval(1,i) = sum_rem;
            if ~isempty(find(reg_interval(1,1:i-1)== reg_interval(1,i)))
                error('There are two same intervals')  % add by ke to remove possible same intervals
            end
        else
            if(sum_rem-intervals(end)*(Variable.CycleLength-i)<intervals(1))
                lower = intervals(1);
            else
                lower = sum_rem-intervals(end)*(Variable.CycleLength-i);
            end
            if (sum_rem-intervals(1)*(Variable.CycleLength-i)>intervals(end))
                upper = intervals(end);
            else
                upper = sum_rem-intervals(1)*(Variable.CycleLength-i);
            end
            new_intervals = setdiff(lower:20:upper,reg_interval); % add by ke to remove possible same intervals
            reg_interval(1,i) = new_intervals(randi(length(new_intervals)));
        end
    end
end
     