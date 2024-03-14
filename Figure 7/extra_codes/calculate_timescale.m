function [timescale, lags] = calculate_timescale(wt_exc_ACF, bin)

timescale = [];
lags =[];
iter = 1;
for i = 1:size(wt_exc_ACF,1)
    iter = 1;
    while wt_exc_ACF(i,iter)>0
        iter = iter + 1;
        if iter == size(wt_exc_ACF, 2)
            break
        end
    end
    iter = iter-1; % iter is the first point fall below zero
    if iter <=0
        warnning('something is wrong')
        timescale(i) =NaN;
        lags(i)      = NaN;
    else
        timescale(i) = sum(wt_exc_ACF(i,2:iter));
        lags(i)      = iter;

    end
end
timescale = timescale * bin;