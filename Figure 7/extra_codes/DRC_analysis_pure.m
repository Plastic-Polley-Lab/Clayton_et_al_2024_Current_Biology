clear
mydir  = pwd;
file = dir('*-adj.mat');
load(file.name)
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)
keep = find([spikedata.keep] ==1);

for i = 1:length(keep)
    neuron_num = keep(i);
    drc_result(neuron_num) = drc_analysis_pure_DRC(neuron_num, spikedata, 0);
end
save('summary_drc.mat','drc_result')

% plot the drc responses
threshold = 0.05
for neuron_num = 1:length(drc_result)
    if drc_result(neuron_num).corrR > threshold
        raster_pure_drcs(drc_result, neuron_num)
        suptitle([animalID, '-Unit', num2str(neuron_num)])
        if neuron_num < 10
            print(['Resp','_',animalID, '-DRC', '_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
            %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
        else
            print(['Resp','_',animalID, '-DRC', '_response_', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
            %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
        end
        close
        
    end
    
    
end
%% Get the correlation-coefficient of the pure DRCs
function drc_result = drc_analysis_pure_DRC(neuron_num, spikedata, fig)


innerIndexes = spikedata(neuron_num).clusterData.stimData.innerIndexes;
inner_var = spikedata(neuron_num).clusterData.stimData.inner_variables;
drc1_indx = [];
for drc_set = 1:3
    psth_summary(drc_set) = drc_summary(drc_set, neuron_num, innerIndexes, spikedata, 0);
end

% Use cross-correlation to calculate the temporal precision
for drc_set = 1:3
    crossparam.timeOnset = 250; % DRC starts at 250 ms, and each chord last 50 ms
    crossparam.timeWindow = 10;
    result(drc_set).SNR1.CorrCoef = cross_coeff(psth_summary(drc_set).SNR1.scmatrix, crossparam, 1);
    
end
i = 1;
corrR(i) = mean([result(1).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ...
    result(2).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ...
    result(3).(['SNR',num2str(i)]).CorrCoef.corrR_avg]);
drc_result.psth_summary = psth_summary;
drc_result.result       = result;
drc_result.corrR        = corrR;
drc_result.latency_p2t  = spikedata(neuron_num).latency_p2t;
drc_result.keep         = spikedata(neuron_num).keep; 
end
%% Plot the rasters for pure drcs
function raster_pure_drcs(drc_result, neuron_num)
figure;
drc = [1, 2, 5]; % plot locations
for i = 1:length(drc)
    subaxis(4,2, drc(i) , 'sh', 0.2, 'sv', 0.0)
    raster_plot_simple(drc_result(neuron_num).psth_summary(i).SNR1, 1:1250)
    subaxis(4,2, drc(i)+2 , 'sh', 0.2, 'sv', 0.05)
    psth_plot_simple(drc_result(neuron_num).psth_summary(i).SNR1, 1:1250)
end

subaxis(4,2, [6,8] , 'sh', 0.2, 'sv', 0.05)
plot([drc_result(neuron_num).result(1).SNR1.CorrCoef.corrR_avg;... 
    drc_result(neuron_num).result(2).SNR1.CorrCoef.corrR_avg;
    drc_result(neuron_num).result(3).SNR1.CorrCoef.corrR_avg], '-o')
ylabel('CorrCoef')
xticks([0,1,2,3,4])
xlim([0,4])
xticklabels({'', 'DRC1', 'DRC2', 'DRC3', ''})
set(gcf)
set(gcf,'position',[100,200,800,800])
end
