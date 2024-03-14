function RLF_analysis_preprocess(path, region)
% cd(path)
mydir  = pwd;
file = dir('*-RateLevel-adj.mat');
if isempty(file)
    sprintf('No manually adjustment')
    file = dir('*-RateLevel.mat');
    load(file.name)
    keep = 1:length(spikedata);

else
    load(file.name)
    keep = find([spikedata.keep] ==1);

end
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)
% keep = 1:length(spikedata);
switch region
    case 'ACtx'
        base_window = 100; % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms
        test_window = 100; % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms
    case 'MGB'
        base_window = 50; % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms
        test_window = 50; % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms
end

cutoff = 0.001;

for i = 1:length(keep)
    indx = keep(i);
    isi = diff(spikedata(indx).clusterData.spiketimes);
    spikedata(indx).refra_violation = length(find(isi<1));
    spikedata(indx).refra_violation_ratio = spikedata(indx).refra_violation/length(isi);
    delay = spikedata(1).clusterData.stimChans{1, 1}.Delay;
    baseline_window = (delay-base_window+1):delay;
    evoke_window    = (delay+1):(delay + test_window);
    % get the firing rate during baseline and evoke activity
    if isnan(spikedata(indx).refra_violation)
        spikedata(indx).pvalue = 1;
        spikedata(indx).resp   = 0;
    else
        scmatrix = spikedata(indx).clusterData.psth.scmatrix;
        baseline        = 1000 * sum(scmatrix(:,baseline_window),2)/length(baseline_window);
        evoke_activity  = 1000 * sum(scmatrix(:,evoke_window),2)/length(evoke_window);
        [p, h] = ranksum(baseline, evoke_activity);
        if p < cutoff
            if mean(evoke_activity)> mean(baseline)
                spikedata(indx).resp   = 1;
                spikedata(indx).pvalue = p;
            elseif mean(evoke_activity)< mean(baseline)
                spikedata(indx).resp   = -1;
                spikedata(indx).pvalue = p;
            else
                error('Evoked responses are the same the spontaneous')
            end
        else
            spikedata(indx).pvalue = p;
            spikedata(indx).resp   = 0;
        end
    end
end
%%
iter = 1;
keep = [];
for i = 1:length(spikedata)
%     if spikedata(i).resp ==1 
    if isempty(spikedata(i).resp) % try to include both exc and inh
    else
        if (spikedata(i).resp ==1) || (spikedata(i).resp ==-1)
            keep(iter) = i;
            iter = iter + 1;
        end
    end
end

switch region
    case 'ACtx'
        analysis_window = [1:100];  % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms 
    case 'MGB'
        analysis_window = [1:50];  % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms 
end
% analysis_window = [6:500];
% analysis_window = [1:50];  % for cortex; let's choose 100 ms; for MGB; let's choose 50 ms 
iter = 1;
% fra_summary =[];
for i = 1:length(keep)
    try
        neuron_num = keep(i);
        [spikedata(neuron_num).spls,  spikedata(neuron_num).rlf, spikedata(neuron_num).spont, spikedata(neuron_num).psth_summary] = rlf_SU_analysis(neuron_num, spikedata, analysis_window);
%         temp_file = spikedata(neuron_num);
%         temp_file.fra = fra;
%         temp_file.spls = spls;
%         temp_file.freqs = freqs;
%         temp_file.tuning = tuning;
%         temp_file.rlf = rlf;
%         temp_file.spont = spont;
%         fra_summary(iter) = temp_file;
%         %         fra_summary(iter).fra = fra;
%         %         fra_summary(iter).spls = spls;
%         %         fra_summary(iter).freqs = freqs;
%         %         fra_summary(iter).tuning = tuning;
%         %         fra_summary(iter).rlf = rlf;
%         iter = iter + 1;
    catch
        warning('some errors happend, discard that units')
        continue
    end
    suptitle([animalID, '-Unit', num2str(neuron_num)])
    if neuron_num < 10
        print([animalID, '-RLF', '_response_stats_0', num2str(neuron_num)],'-dpdf','-bestfit')
        %             set(gcf, 'Color', 'w')
        %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
    else
        print([animalID, '-RLF', '_response_stats_', num2str(neuron_num)],'-dpdf','-bestfit')
        %             set(gcf, 'Color', 'w')
        %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
    end
    close
end
save('summary_rlf_responses_v2.mat','spikedata')
