%%
function raster_data_pattern(path)
% path = 'E:\Ke_Chen\Processed Data\Rach Recording\KeC08\MGB\061420';
cd(path)
file = dir('*.mat');
fs = 24414.0625
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
            if contains(spikedata(i).runName, 'jitter') || contains(spikedata(i).runName, 'scale')
                f_split = split(spikedata(i).runName, '_');
                id_1 = char(f_split(1)); % jitter or scale
                id_2 = char(f_split(2)); % ascending or descending
                id_3 = char(f_split(end)); % which set
                load(['E:\Ke_Chen\Behavior\Patterns\test\', id_3, '\final_cyc4_reg_rand_',id_3, '.mat'])
                switch id_2
                    case 'ascending'
                       clusterData.run(j).event_interval = intervals.([id_1, '_', 'ascending']); 
                    case 'de'
                       clusterData.run(j).event_interval = intervals.([id_1, '_', 'descending']); 
                    case 'rand'
                       clusterData.run(j).event_interval = intervals.([id_1, '_', 'random']);
                end
                clusterData.run(j).event_temp          = clusterData.run(j).ticks(1) + cumsum([0, clusterData.run(j).event_interval]);
                %align the event to the tdt; the time in PXI is different
                %from the time in tdt; correct for each Impale ticks
                m = 1;
                for k = 1:length(clusterData.run(j).event_temp)
                    if clusterData.run(j).event_temp(k) >= clusterData.run(j).ticks(m) && clusterData.run(j).event_temp(k) < clusterData.run(j).ticks(m+1)
                    clusterData.run(j).event(k) = clusterData.run(j).event_temp(k) - (m-1)*1000 + (clusterData.run(j).ticks(m) - clusterData.run(j).ticks(1));
                    else
                        m = m+ 1;
                        clusterData.run(j).event(k) = clusterData.run(j).event_temp(k) - (m-1)*1000 + (clusterData.run(j).ticks(m) - clusterData.run(j).ticks(1));

                    end
                end
                clusterData = psth_spike_pattern(clusterData,j,100);
            elseif contains(spikedata(i).runName, 'freq')
                load('E:\Ke_Chen\Behavior\Patterns\test\frequence_patterns\frequency_reg_rand_set1_cyc6_rep_10_trial100.mat')
                clusterData.run(j).event_freq = reg_rand_frequency;
                clusterData.run(j).event = clusterData.run(j).ticks(1) + [0:3000:297000];
                clusterData = psth_spike_pattern(clusterData,j,3000);
            elseif contains(spikedata(i).runName, 'regrand') || contains(spikedata(i).runName, 'randrand')
                f_split = split(spikedata(i).runName, '_');
                id_1 = char(f_split(1)); % regrand or randrand
                id_2 = char(f_split(2)); % ascending or descending
                id_3 = char(f_split(end)); % which set
                load(['E:\Ke_Chen\Behavior\Patterns\test\', id_3, '\final_cyc4_regrand_',id_3, '_jitter0.mat'])
                switch id_1
                    case 'regrand'
                        clusterData.run(j).event_interval = t;
                    case 'randrand'
                        clusterData.run(j).event_interval = t_rand_rand;
                end
                clusterData.run(j).event_temp          = clusterData.run(j).ticks(1) + cumsum([0, clusterData.run(j).event_interval]);
                %align the event to the tdt; the time in PXI is different
                %from the time in tdt; correct for each Impale ticks
                m = 1;
                for k = 1:length(clusterData.run(j).event_temp)
                    if clusterData.run(j).event_temp(k) >= clusterData.run(j).ticks(m) && clusterData.run(j).event_temp(k) < clusterData.run(j).ticks(m+1)
                        clusterData.run(j).event(k) = clusterData.run(j).event_temp(k) - (m-1)*1000 + (clusterData.run(j).ticks(m) - clusterData.run(j).ticks(1));
                    else
                        m = m+ 1;
                        clusterData.run(j).event(k) = clusterData.run(j).event_temp(k) - (m-1)*1000 + (clusterData.run(j).ticks(m) - clusterData.run(j).ticks(1));
                        
                    end
                end
                clusterData = psth_spike_pattern(clusterData,j,100);
                
            else
                clusterData = psth_spike(clusterData,j);
            end
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