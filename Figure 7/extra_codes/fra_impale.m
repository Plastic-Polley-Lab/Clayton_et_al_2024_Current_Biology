%% FRA from Impale
function [fra, freqs, spls] = fra_impale(path)
% INPUT:
%       path: data directory
%       chan: which channel to plot
load(path)
addpath(genpath('F:\KeChen\MATLAB\MATLAB-Noise'));
%%
spls = innerSeq.master.values; % inner seq
n_spls = length(spls);
freqs = floor(outerSeq.master.values/1e2)*1e2; %outer seq %just rounding off to octaves recognizable
n_freqs = length(freqs);
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
%%for each spl and freq, channel wise accumulate spikes between 110 and 175 ms (for now)
n_chans = 64;
fra = zeros(n_chans,n_freqs, n_spls);
for i = 1:n_tot_trials
    freq = SCL(i).outerIndex;
    spl = SCL(i).innerIndex;
    for j = 1:length(SCL(i).spike_times)
        if (SCL(i).spike_times(j)>=105 && SCL(i).spike_times(j)<=145)  % MGB: using 105-145; ACtx: using 110-175
            fra(SCL(i).tdtChan(j),freq,spl) = fra(SCL(i).tdtChan(j),freq,spl)+1;
        end
    end
end
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