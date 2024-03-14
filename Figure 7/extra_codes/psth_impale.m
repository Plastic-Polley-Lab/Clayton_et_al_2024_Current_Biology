%% FRA from Impale
function psth = psth_impale(path, window)
% INPUT:
%       path: data directory
load(path)
addpath(genpath('F:\KeChen\MATLAB\MATLAB-Noise'));
%%
spls = innerSeq.master.values; % inner seq
n_spls = length(spls);
n_rep = measParam.stopVal;
n_tot_trials = dacInt.numTotalTrials; %total number of trials - n_spls*n_freqs*n_rep
% n_tot_trials = 520; %error in impale
%channels and spike times
maxCodes = 6; % i think max no. of clusters you could have if you online sorted
for i = 1:n_tot_trials
    ch = SCL(i).ch;
    tdtChan = floor((ch-1) / maxCodes) + 1; % this takes care if you do online sorting
    tdtCode = mod((ch-1), maxCodes) + 1; % if you didn't do online sorting, it would all be 1
    SCL(i).tdtChan = tdtChan;
    SCL(i).tdtCode = tdtCode;
    SCL(i).spike_times = SCL(i).t -SCL(i).t(1);
end
%% for calculating fra
%%for each spl,count spikes between 110 and 210 ms (for now)
n_chans = 64;
psth = [];
for i = 1:n_tot_trials
    if ~isempty(SCL(i).spike_times)
        for j = 1:64
            psth(j).raster(i).ts = SCL(i).spike_times(find(SCL(i).tdtChan ==j));
            psth(j).scmatrix(i,:) = histcounts(psth(j).raster(i).ts, window);  % Here the trial length is set to 1 s. 
            for k = 1:length(stimChans)
                psth(j).stimulus.delay(k) = stimChans{k}.Gate.Delay;
                psth(j).stimulus.width(k) = stimChans{k}.Gate.Width;
            end
        end
    end
end

% rlf = rlf./0.1; % convert to spike rate
%% which channel to plot fra
% chan = 33;
% figure; imagesc(squeeze((fra(chan,:,:))))

%% for raster plot
% n_chans = 16;
% 
% for i=1:length(n_tot_trials)
% 
%     for j = 1:length(SCL(i).spike_times)  % loop through all spikes for the trial
%         if (SCL(i).spike_times(j)>=110 && SCL(i).spike_times(j)<=145)
%             fra(SCL(i).tdtChan(j),freq,spl) = fra(SCL(i).tdtChan(j),freq,spl)+1;
%         end
%     end
%    chan_spike_times = spike_times(idxs) - trial_start_times(j); % align spike time to start time of each trial
%    spikes = histcounts(chan_spike_times,0:1:trial_duration);  % generate the spike raster structure for each trial
%    if isempty(spikes) % no spikes in that trial
%        raster(j,:) = zeros(trial_duration,1); % pad zeros
%    else
%        raster(j,:) = spikes;
%    end
% end