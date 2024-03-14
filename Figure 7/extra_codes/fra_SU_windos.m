function fra = fra_SU_windos(neuron_num, spikedata, window )

psth = spikedata(neuron_num).clusterData.psth;
rep  = spikedata(1).clusterData.stimData.numReps;


spls = spikedata(neuron_num).clusterData.stimData.inner_variables; % inner seq
n_spls = length(spls);
freqs = floor(spikedata(neuron_num).clusterData.stimData.outer_variables/1e2)*1e2; %outer seq %just rounding off to octaves recognizable
n_freqs = length(freqs);
innerIndexs = spikedata(neuron_num).clusterData.stimData.innerIndexes;
outerIndexs = spikedata(neuron_num).clusterData.stimData.outerIndexes;
n_rep = spikedata(1).clusterData.stimData.numTrials/ (n_spls * n_freqs);
%%
fra = zeros(n_freqs, n_spls, n_rep);
k = 1
for i = 1: n_freqs
    for j = 1:n_spls
        idx = find(outerIndexs == i & innerIndexs == j);
        scmatrix = psth.scmatrix(idx,:);
        %         subplot(n_freqs,n_spls,k)
        %         fr_avg = mean(scmatrix);
        %         psth_avg_smooth = smoothts(fr_avg,'b',10);
        %         plot(psth_avg_smooth/0.001)
        %         ylim([0,500])
        %         xlim([0, 100])
        %         k = k+1
        %         title([num2str(freqs(i)), '-', num2str(spls(j))])
        scmatrix_sum = sum(scmatrix);
        delay = psth.stimulus.delay;
        duration = psth.stimulus.width;
        analysis_window = window + delay;
        fra(i,j,:) = sum(scmatrix(:,analysis_window),2);
    end
end