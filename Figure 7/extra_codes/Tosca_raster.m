%%settings
clear
pre = -100;
post = 500;
binsize = 1;


%% load files
path = uigetdir('E:\Ke_Chen\Processed Data\PTCHD1-Project\*.mat','Select the folder that stores sorted spikes')
cd(path)
file = dir('*.mat');
fs = 24414.0625;

%%
for j = 1:length(file)
    load(file(j).name)
    for i = 1:length(spikedata)
        try
            fprintf('processing session # %d spike # %d\n', j, i)
            % calculate psth
            spikes = spikedata(i).spiketimes;
            event  = spikedata(i).stimData.event.onset;
            unit = spike2eventRasteandPSTH_NP(spikes/1000, event, binsize, pre, post);
            spikedata(i).psth = unit;
            template_indx = spikedata(i).templateMin;
            waveform = spikedata(i).spatiotemporalTemplate(:,template_indx);
            spikedata(i).waveform = waveform;
            [~,trough] = min(waveform);
            [~,peak] = max(waveform(trough:end));
            spikedata(i).latency_p2t = (peak-1)/fs * 1000;
            
            
        catch
            warning(['problems ocurs for spike # ', num2str(i)])
        end
    end
    fileName = split(file(j).name, '.');
    save([path, '\',char(fileName(1)),'_raster.mat'], 'spikedata')
end
%%


