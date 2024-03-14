function clusterData = psth_spike(clusterData, run_num)
%INPUT:
%       clusterData: data after align spike from phy to impale
%       run_num:   number of which run

% run_num = 5; % the block number
num_trials = clusterData.run(run_num).stimData.numTrials; % number of trials for each block
trial_times = clusterData.run(run_num).ticks; % start time of each trials (may bigger than actual trial numbers or may smaller when crashed)
if num_trials <= length(trial_times) % this is what usually happends
    trial_start_times = trial_times(1:num_trials);% start time of each trials
else
    % this is when system crahsed during recording
    warning('Recording system seems crashed during recording')
    trial_start_times = trial_times(1:(length(trial_times)-1)); % here the last trial was discared
end
trial_duration = round(trial_start_times(2) - trial_start_times(1)); % duration of trials
spike_times = clusterData.run(run_num).spiketimes; % timestamps for each spike
raster=[]; % for storing raster
stimChans = clusterData.run(run_num).stimChans;
psth = [];
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
   psth.raster(j).ts = chan_spike_times;
   spikes = histcounts(chan_spike_times,0:1:trial_duration);  % generate the spike raster structure for each trial
   if isempty(spikes) % no spikes in that trial
       raster(j,:) = zeros(trial_duration,1); % pad zeros
   else
       raster(j,:) = spikes;
   end
   psth.scmatrix = raster;
   for k = 1:length(stimChans)
       psth.stimulus.delay(k) = stimChans{k}.Delay;
       psth.stimulus.width(k) = stimChans{k}.Width;
   end
end
clusterData.run(run_num).psth = psth;