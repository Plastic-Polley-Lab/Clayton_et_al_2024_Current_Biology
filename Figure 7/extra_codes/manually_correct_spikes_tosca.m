function manually_correct_spikes_tosca(path)

load(path)
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
n_spike = length(spikedata);
for i = 1:n_spike
    figure(100)
%     subplot(2,1,1)
%     imagesc(spikedata(i).clusterData.psth.scmatrix)
%     colormap(BW)
%     subplot(2,1,2)
%     plot(smoothts(mean(spikedata(i).clusterData.psth.scmatrix,1), 'b',10))
    try
    psth = spikedata(i).psth;
    
    
    psth_plot_tosca(psth, i, 0)
    fr_avg= (sum(psth.scmatrix(:))/size(psth.scmatrix, 1))/(size(psth.scmatrix, 2)/1000);
    fprintf('Mean firing rate:%4.2f\n', fr_avg)
    set(gcf,'position',[100,200,600,800])
    x = input(['Keep or Discard', ' ', num2str(i), '/', num2str(n_spike), ': 1 or 0\n'])
    switch x
        case 1
            spikedata(i).keep = 1;
        case 0
            spikedata(i).keep = 0;
    end
    catch
        spikedata(i).keep = 0;
        warning('Something is wrong')
    end
end
save([path(1:end-4), '-adj.mat'], 'spikedata')