function pair_pulse_PreProcess_Auto(base_window, test_window, cutoff)

file = dir('*-adj.mat');
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)

load(file.name)
keep = find([spikedata.keep] ==1);
clear summary
for i = 1:length(keep)
    neuron_num = keep(i);
    summary(neuron_num) = paired_pulse_preprocess(spikedata, neuron_num, 0);
    
%     suptitle([animalID, '-Unit', num2str(neuron_num)])
    if neuron_num < 10
%         print([animalID, '-paired_pulse', '_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
        %             set(gcf, 'Color', 'w')
    else
%         print([animalID, '-paired_pulse', '_response_', num2str(neuron_num)],'-dpdf','-bestfit')
        %             set(gcf, 'Color', 'w')
        %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
    end
    close
end
save('summary.mat', 'summary')

%% Test whether neuron is responsive to noise burst, I am using the noise burst ITI = 512, 362, 256 ms; 75 trials
clearvars -except summary spikedata animalID base_window test_window cutoff
% here I use 50 ms before the stimuli as baseline, 50 ms after the stimuli
% as evoke activity; stats were done with ranksum with p < 0.01
if nargin < 3
    base_window = 50;
    test_window = 50; 
    cutoff = 0.01;
end
for i = 1:length(spikedata)
    isi = diff(spikedata(i).clusterData.spiketimes);
    spikedata(i).refra_violation = length(find(isi<1));
    spikedata(i).refra_violation_ratio = spikedata(i).refra_violation/length(isi);
    delay = spikedata(1).clusterData.stimChans{1, 1}.Delay;
    baseline_window = (delay-base_window+1):delay;
    evoke_window    = (delay+1):(delay + test_window);
    % get the firing rate during baseline and evoke activity
    if isnan(spikedata(i).refra_violation) || spikedata(i).keep ==0
        spikedata(i).pvalue = 1;
        spikedata(i).resp   = 0;
    else
        % here using the noise burst ITI = 512, 362, 256 ms
        f = fieldnames(summary);
        scmatrix = [summary(i).(f{end-2}).scmatrix; summary(i).(f{end-1}).scmatrix; summary(i).(f{end}).scmatrix];
        baseline        = 1000 * sum(scmatrix(:,baseline_window),2)/length(baseline_window);
        evoke_activity  = 1000 * sum(scmatrix(:,evoke_window),2)/length(evoke_window);
        [p, h] = ranksum(baseline, evoke_activity);
        if p < cutoff
            if mean(evoke_activity)> mean(baseline)
                spikedata(i).resp   = 1;
                spikedata(i).pvalue = p;
            elseif mean(evoke_activity)< mean(baseline)
                spikedata(i).resp   = -1;
                spikedata(i).pvalue = p;
            else
                error('Evoked responses are the same the spontaneous')
            end
        else
            spikedata(i).pvalue = p;
            spikedata(i).resp   = 0;
        end
    end
end
%% only perform decoding for neurons with noiseburst responses
clearvars -except summary spikedata animalID
innerVar = spikedata(1).clusterData.stimData.inner_variables; % for ploting
decoding_results = [];
for i = 1 : length(spikedata)
    if spikedata(i).resp ==1
        figure
        subplot(1,3,1)
        paired_pulse_preprocess(spikedata, i, 1);        
        subplot(1,3,2)
        options.binsize = 1;
        options = single_neuron_decoding_PSTH(summary(i), options);
        decoding_results = [decoding_results, options];
        yticks([1:2:20])
        labels = cellstr(num2str(flipud(innerVar')));
        yticklabels(labels(1:2:end))
        xticks([2:2:20])
        labels = cellstr(num2str(innerVar'));
        xticklabels(labels(2:2:end))
        box off
        set(gca,'TickDir','out')
        set(gca,'fontsize',12)
        set(gca,'TickLengt', [0.015 0.015]);
        set(gca, 'LineWidth',1)
%         set(gcf,'position',[100,200,1000,800])
        subplot(1,3,3)
        plot(1:length(options.performance), options.performance,'-k')
        hold on
        scatter(1:length(options.performance), options.performance,'ok','filled')
        xticks([2:2:20])
        xticklabels(labels(2:2:end))
        xlabel('Interval (ms)')
        ylabel('Decoding Performance')
        set(gca,'TickDir','out')
        set(gca,'fontsize',12)
        set(gca,'TickLengt', [0.015 0.015]);
        set(gca, 'LineWidth',1)
        set(gcf,'position',[100,200,1200,600])

        neuron_num = i;
        suptitle([animalID, '-Unit', num2str(neuron_num)])
        if neuron_num < 10
            print([animalID, '-paired_pulse_decoding', '_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
        else
            print([animalID, '-paired_pulse_decoding', '_response_', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
            %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
        end
        close
        
    end
end
%
save('summary.mat', 'spikedata', 'summary', 'decoding_results')