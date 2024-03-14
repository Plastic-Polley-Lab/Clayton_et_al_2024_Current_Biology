function psth1 = drc_summary(drc_set, neuron_num, innerIndexes, spikedata, fig)
%INPUT:
%      fig: 1 is to plot the figure, 0 is not to plot the figure
%% for only KeC08; in which the DRCs are not interleaved
% drc1_indx = [];
% switch drc_set
%     case 1
%         temp = 1:7;
%     case 2
%         temp = 8:14;
%     case 3
%         temp = 15:21;
% end
% for i = 1:20
%     drc1_indx = [drc1_indx,temp + 21 * (i-1)]; 
% end


% drc1_indx = drc_set: 3: length(innerIndexes);
drc1_indx = drc_set:3:60;
drc1_indx_temp = [];
for i = 1: length(drc1_indx)
     drc1_indx_temp = [drc1_indx_temp, find(spikedata(1).clusterData.stimData.repIndexes == drc1_indx(i))];
end
drc1_indx = drc1_indx_temp;
%%
drc1_innerIndex = innerIndexes(drc1_indx);
% SNR = 20:-5:-10;
if fig ==1
    figure
end
for i = 1:max(drc1_innerIndex) % each inner index represents a noise level
    indx = find(drc1_innerIndex==i);
    psth_idx = drc1_indx(indx);
    psth = spikedata(neuron_num).clusterData.psth;
    psth1.(['SNR',num2str(i)]).raster = psth.raster(psth_idx );
    psth1.(['SNR',num2str(i)]).scmatrix = psth.scmatrix(psth_idx ,:);
    psth1.(['SNR',num2str(i)]).stimulus = psth.stimulus;
    if fig ==1
        DRC_plot(psth1.(['SNR',num2str(i)]),1, 1:1250, i)
    end
    % psth_plot(psth1.SNR0,1, 1:1250)
    % saveas(gcf, [['DRC1', '_response_', num2str(inner_var(i)), '.pdf']])
    % savefig(['DRC3', '_response_', num2str(inner_var(i)), '.fig'])
    % close all
end