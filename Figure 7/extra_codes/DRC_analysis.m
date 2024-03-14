%% calculate the sound level for the DRC audio file
% load('Test1-1-4-DRC_set2_no_noise.mat')
% mask= stimChans{1}.Source.Mask;
% for i = 1:size(mask, 1)
%     level = mask(i,(mask(i,:) > 0 ));
%     level_avg(i) = 20*log10(sqrt(sum((10.^(level/20)).^2))); % calculate the overall sound intensity
% end
% % there is 250 ms gating, the first 5 time slices are skipped
% for i = 1:3 % there are three DRCs
%     drc(i).idx = 6 + (i-1) * 25 :  i * 25;
%     drc(i).level = level_avg(drc(i).idx);
%     drc(i).level_max = max(drc(i).level);
% end
%% DRC pure
% edit DRC_analysis_pure


%% fano factor analysis
% clear
% animalID = 'KeC17'
clear
mydir  = pwd;
file = dir('*-adj.mat');
load(file.name)
idcs2 = strfind(file.name,'_');
animalID = file.name(1:idcs2(1)-1)
keep = find([spikedata.keep] ==1);
%
% resp.Resp = intersect(keep, resp_idx);
% resp.NonResp = intersect(keep, setdiff(1:length(spikedata), resp_idx));
% type = {'Resp', 'NonResp'};
%     drc_result = [];
threshold = 0.05;
keep_idx = keep;
SNR = 20:-5:-10;
for i = 1:length(keep_idx)
    neuron_num = keep_idx(i);
    drc_result(neuron_num) = DRC_process(neuron_num, spikedata, threshold, SNR);
    if drc_result(neuron_num).resp == 1
        drc_parallel_response_fanofactor(neuron_num, drc_result, SNR)
        %     if drc_result(neuron_num).resp == 1
        type = 'Resp';
        suptitle([animalID, '-Unit', num2str(neuron_num)])
        if neuron_num < 10
            print([type,'_',animalID, '-DRC', '_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
            %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
        else
            print([type,'_',animalID, '-DRC', '_response_', num2str(neuron_num)],'-dpdf','-bestfit')
            %             set(gcf, 'Color', 'w')
            %             export_fig([animalID, '-nDRC', '_response_', num2str(neuron_num)],  '-png')
        end
        close
    else
        %         type = 'NonResp';
    end
    
    save('summary_drc.mat','drc_result')
end
clear drc_result


%% summary analysis
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
genotype = {'wt', 'ko'}
for k = 1:length(genotype)
    switch genotype{k}
        case 'wt'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091820\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\092020\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092420\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092520\DRC_set2_noise', ...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100720\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110620\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110720\DRC_set1_noise'};
            
            
        case 'ko'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091620\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091720\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092120\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092220\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092320\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092620\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092720\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092820\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\103020\DRC_set1_noise'};
    
    
    end
    temp =[];
    clear drc_result
    for i = 1 : length(paths)
        load([paths{i}, '\summary_drc.mat'])
        for j = 1:length(drc_result)
            drc_result(j).paths = paths{i};
        end
        temp = [temp, drc_result];
    end
    save(['drc_results_', genotype{k},'.mat'],'temp', '-v7.3')
end
%% Let's keep load data
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
genotype = {'wt', 'ko'}
for i = 1:length(genotype)
    switch genotype{i}
        case 'wt'
            newData = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\DRC_set1_noise\', ...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\DRC_set1_noise'};
        case 'ko'
            newData = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\DRC_set1_noise\', ...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111220\DRC_set1_noise'};
    end
    load(['drc_results_', genotype{i},'.mat'])
    
    for j = 1:length(newData)
        
        load([newData{j}, '\summary_drc.mat'])
        for k = 1:length(drc_result)
            drc_result(k).paths = newData{j};
        end
        temp = [temp, drc_result];
        clear drc_result
    end
    save(['drc_results_', genotype{i},'.mat'],'temp', '-v7.3')
    clear temp
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||%
%||%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%||%
%% summarize the data for plots
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\DRC_analysis')
path = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\DRC_analysis';

% summaryData = summary_DRC_analysis(path, 0.05, 'ACtx');
summaryData = summary_DRC_analysis(path, 0.1, 'MGB');

%% Plot the raw data
figure
line_errorbar_drc(summaryData.wt.CrossCoef, summaryData.ko.CrossCoef)
set(gcf, 'Color', 'w')


% export_fig('DRC_summary_MGB01',  '-png', '-pdf')

% i = 1
% [h, p] = ttest2(summary.wt.CrossCoef(:,1), summary.ko.CrossCoef(:,1))
% % plot the normalized data, normalize to the highest SNR
% figure;
% line_errorbar_drc(summary.wt.CrossCoef_n, summary.ko.CrossCoef_n)
% set(gcf, 'Color', 'w')
% export_fig('DRC_normalized_summary_MGB01',  '-png', '-pdf')
% 
% % plot the slope of the changes
% figure;
% line_errorbar_drc(diff(summary.wt.CrossCoef, 1,2), diff(summary.ko.CrossCoef,1,2))
% ylabel('\Delta CorrCoef')
% set(gcf, 'Color', 'w')
% export_fig('DRC_slop_summary_MGB01',  '-png', '-pdf')

%%
% xlswrite(filename, diff(CrossCoef,1, 2))
% xlswrite(filename, CrossCoef_norm)

%%
% for i = 2
%     switch i
%         case 1
%             animalID = 'KeC08'
%             load('E:\Ke_Chen\Processed Data\Rach Recording\KeC08\MGB\061420\DRC_set1_noise\summary_drc.mat')
%             neuron_num = [8,12,21];
%             neuron_num = 8
%
%         case 2
%             animalID = 'KeC09'
%             neuron_num = [20,24]
% %             neuron_num = [4,7,8,12,15,16,20,21,22,24,26,30,31,33,34,35,37,39];
%             load('E:\Ke_Chen\Processed Data\Rach Recording\KeC09\061920\DRC_set2_noise\summary_drc.mat')
%
%         case 3
%             animalID = 'KeC10'
%             neuron_num = [16, 17, 18, 19,25];
%             load('E:\Ke_Chen\Processed Data\Rach Recording\KeC10\062720\DRC_set1_noise\summary_drc.mat')
%
%         case 4
%             animalID = 'KeC12'
%             neuron_num = [14, 15, 37, 40, 41, 43, 45, 47, 49, 50, 51]
%             load('E:\Ke_Chen\Processed Data\Rach Recording\KeC12\070920\DRC_set2_noise\summary_drc.mat')
%
%     end
%     for i = 1:length(neuron_num)
%         figure
%         drc_parallel_response_fanofactor(neuron_num(i), drc_result)
%         suptitle([animalID, '-Unit', num2str(neuron_num(i))])
%         if neuron_num(i) < 10
% %             print([animalID, '-DRC', '_response_0', num2str(neuron_num(i))],'-dpdf','-bestfit')
%             set(gcf, 'Color', 'w')
%             export_fig([animalID, '-DRC', '_response_', num2str(neuron_num(i))],  '-png')
%         else
% %             print([animalID, '-DRC', '_response_', num2str(neuron_num(i))],'-dpdf','-bestfit')
%             set(gcf, 'Color', 'w')
%             export_fig([animalID, '-DRC', '_response_', num2str(neuron_num(i))],  '-png')
%         end
% %         close
%     end
% end
%%





%% Get the correlation-coefficient of the pure DRCs
function drc_result = drc_analysis_pure_DRC(neuron_num, spikedata, fig)


innerIndexes = spikedata(neuron_num).clusterData.stimData.innerIndexes;
inner_var = spikedata(neuron_num).clusterData.stimData.inner_variables;
drc1_indx = [];
for drc_set = 1:3
    psth_summary(drc_set) = drc_summary(drc_set, neuron_num, innerIndexes, spikedata, 1);
end

% Use cross-correlation to calculate the temporal precision
for drc_set = 1:3
    crossparam.timeOnset = 250; % DRC starts at 250 ms, and each chord last 50 ms
    crossparam.timeWindow = 10;
    result(drc_set).SNR1.CorrCoef = cross_coeff(psth_summary(drc_set).SNR1.scmatrix, crossparam, 1);
    
end
i = 1
corrR(i) = mean([result(1).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ...
    result(2).(['SNR',num2str(i)]).CorrCoef.corrR_avg, ...
    result(3).(['SNR',num2str(i)]).CorrCoef.corrR_avg]);
drc_result.psth_summary = psth_summary;
drc_result.result       = result;
drc_result.corrR        = corrR;
drc_result.latency_p2t  = spikedata(neuron_num).latency_p2t;
drc_result.keep         = spikedata(neuron_num).keep;
end