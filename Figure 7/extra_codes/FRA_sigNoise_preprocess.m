function FRA_sigNoise_preprocess(path)
cd(path)
mydir  = pwd;
file = dir('*-FRA_CF_50-adj.mat');
if isempty(file)
    sprintf('No manually adjustment')
    file = dir('*-FRA_CF_50.mat');
    load(file.name)
    keep = 1:length(spikedata);
    
else
    load(file.name)
    keep = find([spikedata.keep] ==1);
    
end
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)
% keep = 1:length(spikedata);
resp_indx = [];
iter = 1;
for i = 1:length(keep)
    neuron_num = keep(i);
    isi = diff(spikedata(neuron_num).clusterData.spiketimes);
    spikedata(neuron_num).refra_violation = length(find(isi<1));
    spikedata(neuron_num).refra_violation_ratio = spikedata(neuron_num).refra_violation/length(isi);
    spikedata(neuron_num).stats = response_test(spikedata(neuron_num));
    spikedata(neuron_num).resp  = spikedata(neuron_num).stats.resp_sign;
    if spikedata(neuron_num).stats.resp ==1
        resp_indx(iter) = neuron_num;
        iter = iter + 1;
    end
end
%%
% analysis_window = [6:500];
analysis_window = [1:50]; % for MGB
analysis_window = [1:500]; % for cortex

iter = 1;
fig = 1; % plot the data or not
% fra_summary =[];
for i = 1:length(resp_indx)
    try
        neuron_num = resp_indx(i);
        [fra, spls, freqs, tuning, rlf, spont] = fra_SU_analysis(neuron_num, spikedata, analysis_window,fig);
        temp_file = spikedata(neuron_num);
        temp_file.fra = fra;
        temp_file.spls = spls;
        temp_file.freqs = freqs;
        temp_file.tuning = tuning;
        temp_file.rlf = rlf;
        temp_file.spont = spont;
        fra_summary(iter) = temp_file;
        %         fra_summary(iter).fra = fra;
        %         fra_summary(iter).spls = spls;
        %         fra_summary(iter).freqs = freqs;
        %         fra_summary(iter).tuning = tuning;
        %         fra_summary(iter).rlf = rlf;
        iter = iter + 1;
    catch
        warning('some errors happend, discard that units')
        continue
    end
    if fig ==1
        suptitle([animalID, '-Unit', num2str(neuron_num)])
        if neuron_num < 10
            %         print([animalID, '-Freq', '_response_stats_0', num2str(neuron_num)],'-dpdf','-bestfit')
            set(gcf, 'Color', 'w')
            export_fig([animalID, '-Freq', '_response_0', num2str(neuron_num)],  '-pdf')
        else
            %         print([animalID, '-Freq', '_response_stats_', num2str(neuron_num)],'-dpdf','-bestfit')
            set(gcf, 'Color', 'w')
            export_fig([animalID, '-Freq', '_response_', num2str(neuron_num)],  '-pdf')
        end
    end
    close
end
save('summary_fra_responses_v2.mat','fra_summary')

%% calculate the responses based on multiple bins
function stats = response_test(spikedata)
base_window = 100; % for cortex; let's choose 100 ms; for MGB: let's choose 50 ms
test_window = 500; % for cortex; let's choose 500 ms; for MGB: let's choose 50 ms

spike = spikedata(1).clusterData.psth.scmatrix;


delay = spikedata(1).clusterData.stimChans{1, 1}.Delay;
% baseWindow = 50;
consecutive_binsize = 100; % use two consecutive bins, each bin is 50 ms
n_consecutive = 5; % 5 consecutive bins to cover the duration of the responses
baseline_window = (delay-base_window + 1) :delay;
for i = 1:n_consecutive
    evoke_window{i} = (delay+1 + consecutive_binsize*(i-1)):(delay + consecutive_binsize *i);
end

%stats parameter
alpha_value = 0.01;
adj_alpha = 1 - (1-alpha_value)^(1/n_consecutive);

for i = 1:n_consecutive
    baseline_activity = sum(spike(:,baseline_window),2)/length(baseline_window)* 1000;
    evoke_activity = sum(spike(:,evoke_window{i}),2)/length(evoke_window{i})* 1000;
    [p(i), h(i)] = ranksum(baseline_activity, evoke_activity);
    
    if mean(evoke_activity)> mean(baseline_activity)
        sign(i) = 1;
    elseif mean(evoke_activity)< mean(baseline_activity)
        sign(i)   = -1;
    else
        sign(i) = 0;
        warning('Evoked responses are the same the spontaneous')
    end
end
significant = find(p < adj_alpha);
if isempty(significant)
    resp = 0;
    resp_sign = 0;
else
    resp = 1;
    resp_sign = sign(significant(1)); % choose the first significant response to define the sign
end

stats.p = p;
stats.sign = sign;
stats.resp = resp;
stats.resp_sign = resp_sign;
