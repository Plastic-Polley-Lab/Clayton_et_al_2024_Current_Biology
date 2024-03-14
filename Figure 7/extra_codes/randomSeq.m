function reg_interval = randomSeq(intervals, CycleLength)
%INPUT:
%     intervals: the available values to pick from,
%     CycleLength: the total number of value to pick
%OUTPUT:
%    reg_interval:
%                the random sequences with the mean as intervals.
rng('shuffle');
% intervals = [-10:10];
% CycleLength = 4;
cycle_duration = mean(intervals)*CycleLength;
reg_interval = [];
for i = 1:CycleLength
    if (i == 1)
        reg_interval(1,i) = intervals(randi(length(intervals)));
    else
        sum_rem = cycle_duration - sum(reg_interval(1,1:i-1));
        if (i == CycleLength)
            reg_interval(1,i) = sum_rem;
%             if ~isempty(find(reg_interval(1,1:i-1)== reg_interval(1,i)))
%                 error('There are two same intervals')  % add by ke to remove possible same intervals
%             end
        else
            if(sum_rem-intervals(end)*(CycleLength-i)<intervals(1))
                lower = intervals(1);
            else
                lower = sum_rem-intervals(end)*(CycleLength-i);
            end
            if (sum_rem-intervals(1)*(CycleLength-i)>intervals(end))
                upper = intervals(end);
            else
                upper = sum_rem-intervals(1)*(CycleLength-i);
            end
            new_intervals = lower:1:upper; % add by ke to remove possible same intervals
            reg_interval(1,i) = new_intervals(randi(length(new_intervals)));
        end
    end
end