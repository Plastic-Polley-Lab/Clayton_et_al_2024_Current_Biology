function Generate_Raster(clusterData,run_num)
%load('Rach_Ephys-Sorted_Data-ACtx-cluster_2.mat')
clear raster
% run_num = 1; % the block number
num_trials = clusterData.run(run_num).stimData.numTrials; % number of trials for each block
trial_times = clusterData.run(run_num).ticks; % start time of each trials (may bigger than actual trial numbers)
trial_start_times = trial_times(1:num_trials);% start time of each trials
trial_duration = round(trial_start_times(2) - trial_start_times(1)); % duration of trials
spike_times = clusterData.run(run_num).spiketimes; % timestamps for each spike
raster=[]; % for storing raster
for j=1:length(trial_start_times)
   if j<length(trial_start_times)
       current_tick = trial_start_times(j); % get each trial start time
       next_tick = trial_start_times(j+1);  % get next trial start time
       idxs = find(spike_times>=current_tick & spike_times<next_tick); % get index of spikes within a trial
   else
       current_tick = trial_start_times(j); % get last trial start time
       idxs = find(spike_times>=current_tick); % get spikes after last trial starts
   end
   chan_spike_times = spike_times(idxs) - trial_start_times(j); % align spike time to start time of each trial
   spikes = histcounts(chan_spike_times,0:1:trial_duration);  % generate the spike raster structure for each trial
   if isempty(spikes) % no spikes in that trial
       raster(j,:) = zeros(trial_duration,1); % pad zeros
   else
       raster(j,:) = spikes;
   end
end
%% calculate the fra
load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
figure;
subplot(2,1,1)
imagesc(raster),colormap(BW)
xlim([0,500])
box off; axis off
title(['ClusterID ', num2str(clusterData.run(1).clusterID)])
subplot(2,1,2)
spike_counts = sum(raster,1)/size(raster,1);
stem(spike_counts,'k')
onset = clusterData.run(run_num).stimChans{1, 1}.Delay;
duration = clusterData.run(run_num).stimChans{1, 1}.Width
hold on
rectangle('Position',[onset,0,duration,max(spike_counts)],'EdgeColor','r')
hold off
xlim([0,500])

% %%
% num_trials = clusterData.run(run_num).stimData.numTrials; % number of trials for each block
% trial_times = clusterData.run(run_num).ticks; % start time of each trials (may bigger than actual trial numbers)
% trial_start_times = trial_times(1:num_trials);
% trial_duration = round(trial_start_times(2) - trial_start_times(1));
% spike_times = clusterData.run(run_num).spiketimes;
% load('F:\KeChen\RawData\Rach_Ephys\KeC022520\Impale\KeC022520-1-32-FRA.mat');
% clusterData.run(run_num).stimData.outerIndexes = [SCL.outerIndex];
% clusterData.run(run_num).stimData.innerIndexes = [SCL.innerIndex];
% %     end
% spls = innerSeq.master.values; % get the SPL
% num_spls = length(spls);
% freqs = outerSeq.master.values; % get the frequency
% num_freqs = length(freqs);
% fra = zeros(num_freqs, num_spls);  % generate an empty matrix for store spike count for ploting fra
% for j=1:length(trial_start_times) % loop through each trial
%     %     if (length(clusterData.run(run_num).stimData.outerIndexes)<num_trials) % trial number should match
%     %         warning('something unexpected')
%     freq = clusterData.run(run_num).stimData.outerIndexes(j);
%     spl = clusterData.run(run_num).stimData.innerIndexes(j);
%     if j<length(trial_start_times)
%         current_tick = trial_start_times(j);
%         next_tick = trial_start_times(j+1);
%         idxs = find(spike_times>=current_tick & spike_times<next_tick); % get the index of spike time that is inside the trial window
%     else
%         current_tick = trial_start_times(j);
%         idxs = find(spike_times>=current_tick);
%     end
%     chan_spike_times = spike_times(idxs) - trial_start_times(j); % re-align to trial start
%     fra(freq,spl) = fra(freq,spl)+ histcounts(chan_spike_times,105:50:155); % count spikes within [105, 155]
% end
% 
% fra = fra./max(fra(:)); % normalized to max spike count
% figure; imagesc(fra); colormap(hot)
% figure; imagesc(flipud(fra')); colormap(hot)
% % xticks([1:length(freqs)])
% % for i = 1:length(freqs)
% %     labels{i} = num2str(freqs(i));
% % end
% % xticklabels(labels)
