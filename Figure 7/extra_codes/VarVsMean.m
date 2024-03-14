function result = VarVsMean(scmatrix, fanoparam, fig)
% by default there are 20 chords, and each chord last 50 ms
% let's get the index of each chord
idx = 1: 50 : 20 * 50 ;
% let's store the mean spike count and variane for each chord
spikecount_m = [];
spikecount_var = [];
timeOnset = fanoparam.timeOnset;
timeWindow = fanoparam.timeWindow;
for i = 1: length(idx)
    temp_idx = (timeOnset + idx(i)) : (timeOnset + idx(i) + timeWindow);
    temp = sum(scmatrix(:,temp_idx),2);
    result.spikecount_m(i) = mean(temp);
    result.spikecount_var(i) = var(temp);
end

if fig ==1
    figure
    h1 = scatter(result.spikecount_m, result.spikecount_var, 'k', 'filled')
    xlabel('Mean Spike Count')
    ylabel('Spike Count Variance')
end

mdl = fitlm(result.spikecount_m,result.spikecount_var, 'Intercept',false);
result.fanofactor = mdl.Coefficients.Estimate;

if fig ==1
    hold on
    h2 = plot(result.spikecount_m, mdl.Fitted, 'r')
    text(mean(result.spikecount_m), mean(mdl.Fitted), ['Slope = ', num2str(result.fanofactor)])
end
