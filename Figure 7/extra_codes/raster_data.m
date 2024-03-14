%%
function raster_data(path)
% path = 'E:\Ke_Chen\Processed Data\Rach Recording\KeC08\MGB\061420';
cd(path)
file = dir('*.mat');
fs = 24414.0625;
%% start to organize all spike data

spikedata = [];
load(file(1).name)
n_run = length(clusterData.run)
for j = 1:n_run
    
    for i = 1:length(file)
        try
            fprintf('processing spike # %d\n', i)
            spikedata(i).name = [file(i).folder,'\',file(i).name];
            load(file(i).name)
            f = clusterData.run(j).filename;
            f_s1 = split(f,'\'); % 1st time split
            f_s2 = split(char(f_s1(end)),'-');
            f_runname = char(f_s2(end));
            spikedata(i).runName = f_runname(1:end-4);
            if exist(spikedata(i).runName) == 7
            else
                mkdir(spikedata(i).runName)
            end
            clusterData = psth_spike(clusterData,j);
            spikedata(i).clusterData = clusterData.run(j);
            template_indx = clusterData.run(1).templateMin;
            waveform = clusterData.run(1).spatiotemporalTemplate(:,template_indx);
            spikedata(i).waveform = waveform;
            [~,trough] = min(waveform);
            [~,peak] = max(waveform(trough:end));
            spikedata(i).latency_p2t = (peak-1)/fs * 1000;
        
        
        catch
%             spikedata(i).clusterData = clusterData;
%             template_indx = clusterData.run(1).templateMin;
%             waveform = clusterData.run(1).spatiotemporalTemplate(:,template_indx);
%             spikedata(i).waveform = waveform;
%             [~,trough] = min(waveform);
%             [~,peak] = max(waveform(trough:end));
%             spikedata(i).latency_p2t = (peak-1)/fs * 1000;
            warning(['problems ocurs for run # ', num2str(j)])
        end
    end
    save([path, '\',spikedata(i).runName,'\', char(f_s2(1)),'-', char(f_s2(end))], 'spikedata')
    
end