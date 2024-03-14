function cp_value = Timing_CP(data)
i = 1;
num_trials = data(i).clusterData.stimData.numTrials; % number of trials for each block
trial_times = data(i).clusterData.ticks; % start time of each trials (may bigger than actual trial numbers or may smaller when crashed)
if num_trials <= length(trial_times) % this is what usually happends
    events = trial_times(1:num_trials);% start time of each trials
else
    % this is when system crahsed during recording
    warning('Recording system seems crashed during recording')
    events = trial_times(1:(length(trial_times)-1)); % here the last trial was discared
end

spike_times = data(i).clusterData.spiketimes; % timestamps for each spike
delay = data(i).clusterData.stimChans{1, 1}.Delay;
% convert the unit of time from ms to s
events = (events + delay)/1000;
spike_times = spike_times/1000;

p1 = spike2eventRasteandPSTH_NP(spike_times,events,5,-100,400);
[Time_On, Time_Off, H, p, tol, badCP, CP]= CP_detection_lz (p1, -0.1,1);
cp_value.Time_On = Time_On;
cp_value.Time_Off = Time_Off;
cp_value.H = H;
cp_value.p = p;
cp_value.Time_On = Time_On;
cp_value.tol = tol;
cp_value.badCP = badCP;
cp_value.CP = CP;

% figure;
% plot(p1.timepoint, p1.FR_avg)