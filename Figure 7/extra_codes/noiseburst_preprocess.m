function noiseburst_preprocess(path,region,genotype, fig)
cd(path)
mydir  = pwd;
file = dir('*-NoiseBursts-adj.mat');
if isempty(file)
    sprintf('No manually adjustment')
    file = dir('*-NoiseBursts_CF.mat');
    load(file.name)
    keep = 1:length(spikedata);

else
    load(file.name)
    keep = find([spikedata.keep] ==1);

end
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)
% noiseburst or noiseCSD analysis
% check the refractory period violation, and whether the neuron is
% responsive to noise burst or noise CSD
% here I use 50 ms before the stimuli as baseline, 50 ms after the stimuli
% as evoke activity; stats were done with ranksum with p < 0.01
base_window = 50;
test_window = 50; % I use single bin or consecutive bins
cutoff = 0.01; % criteria for responsiveness
% only use the data that is good

spikedata = spikedata(keep);

for i = 1:length(spikedata)
    spikedata(i).animalID = animalID;
    isi = diff(spikedata(i).clusterData.spiketimes);
    spikedata(i).refra_violation = length(find(isi<1));
    spikedata(i).refra_violation_ratio = spikedata(i).refra_violation/length(isi);
    delay = spikedata(1).clusterData.stimChans{1, 1}.Delay;
    baseline_window = (delay-base_window+1):delay;
    evoke_window    = (delay+1):(delay + test_window);
    % get the firing rate during baseline and evoke activity
    if isnan(spikedata(i).refra_violation)
        spikedata(i).pvalue = 1;
        spikedata(i).resp   = 0;
    else
        baseline        = 1000 * sum(spikedata(i).clusterData.psth.scmatrix(:,baseline_window),2)/length(baseline_window);
        evoke_activity  = 1000 * sum(spikedata(i).clusterData.psth.scmatrix(:,evoke_window),2)/length(evoke_window);
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
save('summary_noise', 'spikedata')
%% let's plot all the responsive and non-responsive units
cd(mydir)
if fig == 1
    resp = [spikedata.resp];
    responsive_idx = find(resp ~=0);
    nonresponsive_idx = find(resp ==0);
    plot_noiseburst(spikedata,responsive_idx, 'Responsive')
    plot_noiseburst(spikedata,nonresponsive_idx, 'Non-Responsive')
end

%% Let's summarize the data in one files

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\NoiseBurst')

to_save = ['summary_noise_', genotype, '_' region, '.mat'];
if exist(to_save) ==0
    summarydata.(genotype) = spikedata;
    save(to_save,'summarydata')
elseif exist(to_save) ==2
    load(to_save)
    summarydata.(genotype) = [summarydata.(genotype),spikedata];
    lastwarn('') % This will clear any last warnings so there is no  confusion
    save(to_save,'summarydata') % attempt to save the file in -v7 format
    [msg,id]=lastwarn('');     % This will save the last warning
    if strcmp(id,'MATLAB:save:sizeTooBigForMATFile')    % compare the id of the warning with the one we are looking for
        save(to_save,'-v7.3','summarydata');    % save the variable in a -v7.3 MAT-file
    end
    
    
end
cd(mydir)

% psth = spikedata(neuron_num).clusterData.psth;
%%function for plots
function plot_noiseburst(spikedata, indx, type)
for i = 1: length(indx)
    neuron_num = indx(i);
    psth = spikedata(neuron_num).clusterData.psth;
    psth_plot(psth,neuron_num, 1:size(psth.scmatrix,2), 1)
    peak_firing = max(mean(psth.scmatrix, 1))/0.001;
    text(270,peak_firing/2, sprintf('p value:%1.3f', spikedata(neuron_num).pvalue),'fontsize',12)
    if spikedata(neuron_num).refra_violation_ratio> 0.005
        text(270,peak_firing/2-5, sprintf('ISI violtation:%1.3f', spikedata(neuron_num).refra_violation_ratio),'fontsize',12)
    end
    set(gcf,'position',[100,200,600,600])
    set(gcf, 'Color', 'w')
    if neuron_num < 10
        export_fig([type,'_Noise_response_0', num2str(neuron_num)],  '-pdf')
    else
        export_fig([type, '_Noise_response_', num2str(neuron_num)],  '-pdf')
    end
    close
end
   
