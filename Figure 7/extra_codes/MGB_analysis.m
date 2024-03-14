%% Frequency tuning analysis

path = pwd;
FRA_analysis_preprocess(path)
%% summary analysis
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\FRA_analysis')
genotype = {'WT', 'KO'}
for z = 1:length(genotype)
    
    switch genotype{z}
        case 'WT'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100820\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC29\101520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110320\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110620\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110720\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\FRA_CF'};
        case 'KO'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100620\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100720\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC28\101420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC30\101920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\103020\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111120\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111220\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121320\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\FRA_CF',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\FRA_CF'};
    end
    
    temp =[];
    for i = 1 : length(paths)
        load([paths{i}, '\summary_fra_responses_v2.mat'])
        if isfield(fra_summary, 'keep')
        fra_summary = rmfield(fra_summary, 'keep');
        end
        for j = 1:length(fra_summary)
            fra_summary(j).paths = paths{i};
        end
        temp = [temp, fra_summary];
    end
    save(['fra_tuning_results_MGB_', genotype{z}, '.mat'],'temp', '-v7.3')
end
%% Let's see the frequency tuning
clear
load('fra_tuning_results_MGB_WT.mat')
for i = 1:length(temp)
    [temp(i).cFRA, temp(i).spontDistr, temp(i).fraDistr, temp(i).dprime, temp(i).FRAmask] = dprime_FRA(temp(i).fra);
end
tuning_summary.wt = temp;

load('fra_tuning_results_MGB_KO.mat')
for i = 1:length(temp)
    [temp(i).cFRA, temp(i).spontDistr, temp(i).fraDistr, temp(i).dprime, temp(i).FRAmask] = dprime_FRA(temp(i).fra);
end
tuning_summary.ko = temp;

save('fra_tuning_summary_MGB.mat', 'tuning_summary', '-v7.3')
clearvars -except tuning_summary
%% Let's see the population gain
results.wt = getGain_threshold(tuning_summary.wt);
results.ko = getGain_threshold(tuning_summary.ko);
%% Let's see the population gain

figure;
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    for i = 1:length(results.(genotype{j}))
        plots.(genotype{j}).rlf(i,:) = results.(genotype{j})(i).rlf_avg;
    end
    plots.(genotype{j}).gains = [results.(genotype{j}).gains_value];
    plots.(genotype{j}).threshold = [results.(genotype{j}).threshold];
    nan_indx = find(isnan(plots.(genotype{j}).threshold));
    plots.(genotype{j}).gains(nan_indx)=[];
    plots.(genotype{j}).threshold(nan_indx) =[];
    plots.(genotype{j}).rlf(nan_indx,:) = [];
end



line_errorbar_drc(plots.wt.rlf, plots.ko.rlf)
xticks([1:8])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70'})
xlabel('Sound Intensity(dB)')
ylabel('Firing Rate (Hz)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
export_fig('Rate_Level_Function_MGB',  '-png', '-pdf')

%% Let's see the threshold
figure;

a= plots.wt.threshold';
b = plots.ko.threshold';
b(end+1:length(a)) =  NaN;
figure
[h1, h2] = bar_plot([a, b])
h1.FaceColor = 'k';
h2.FaceColor = 'r';
h1.FaceAlpha = 0.3;
h2.FaceAlpha = 0.7;
ylim([0,50])
xticks([1,2])
xticklabels({'WT', 'KO'})
ylabel('Best Frequence Response Threshold (dB)')

set(gcf, 'Color', 'w')
export_fig('Response_Threshold_MGB',  '-png')
%% Let's check the gains
figure;
[wt.f, wt.x] = ecdf(plots.wt.gains);
hold on
[ko.f, ko.x] = ecdf(plots.ko.gains);
h_wt = plot(wt.x, wt.f, '-k', 'LineWidth', 1);
h_ko = plot(ko.x, ko.f, '-r', 'LineWidth', 1);
legend([h_wt, h_ko], {'WT', 'KO'})
xlabel('Gain (\Delta sp/s per 10 dB step)')
ylabel('Cumulative probability')
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
xlim([0,100])
set(gcf, 'Color', 'w')
export_fig('Gains_MGB',  '-png', '-pdf')
%% Bar plots

a= plots.wt.gains';
b = plots.ko.gains';
[h, p] = ttest2(a, b)
b(end+1:length(a)) =  NaN;
figure
[h1, h2] = bar_plot([a, b])
h1.FaceColor = 'k';
h2.FaceColor = 'r';
h1.FaceAlpha = 0.3;
h2.FaceAlpha = 0.7;
ylim([0,20])
xticks([1,2])
xticklabels({'WT', 'KO'})
ylabel('Gain (\Delta sp/s per 10 dB step)')
set(gcf, 'Color', 'w')
export_fig('Gain_bar_MGB',  '-png', '-pdf')
%%
xlswrite('Gain_MGB_wt.xlsx', wt.gains_value)
xlswrite('Gain_MGB_ko.xlsx', ko.gains_value)
save('fra_tuning_gain_summary_MGB.mat', 'wt', 'ko')

%% MGB DRCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% summary analysis
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\DRC_analysis')
genotype = {'wt', 'ko'}
for k = 1:length(genotype)
    switch genotype{k}
        case 'wt'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100120\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100220\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100820\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100920\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC29\101520\DRC_set1_noise', ...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110120\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110320\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110620\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110720\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\DRC_set2_noise'};
                
            
        case 'ko'
            paths = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100620\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100720\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC28\101420\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC30\101920\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\103020\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111120\DRC_set1_noise', ...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC34\111220\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\DRC_set2_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\DRC_set1_noise',...
                'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\DRC_set2_noise'};
    
    
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
    save(['drc_results_MGB_', genotype{k},'.mat'],'temp', '-v7.3')
end
%%
%%
clear
% load('drc_results_KO.mat')

genotype = {'wt', 'ko'}
for k = 1:length(genotype)
    switch genotype{k}
        case 'wt'
            load('drc_results_MGB_wt.mat')
%             filename = 'drc_wt005_n.xlsx';
        case 'ko'
            load('drc_results_MGB_ko.mat')
%             filename = 'drc_ko005_n.xlsx';
    end
    
    iter = 1;
    CrossCoef =[];
    
    threshold = 0.05;
    
    for i = 1:length(temp)
        if ~isempty(temp(i).corrR)
            %         if temp(i).resp ==1
            if temp(i).corrR(1) >= threshold
                CrossCoef(iter, :) = temp(i).corrR;
                iter = iter + 1;
            end
        end
    end
    keep = [temp.keep];
    resp_ratio = size(CrossCoef, 1)/length(keep);
    summary.(genotype{k}).CrossCoef = CrossCoef;
    summary.(genotype{k}).CrossCoef_n = CrossCoef./CrossCoef(:,1);
    summary.(genotype{k}).keep = keep;
    summary.(genotype{k}).resp_ratio = resp_ratio;
end
save('drc_results_summary_MGB005', 'summary', 'threshold')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%|
%% MGB paired pulse data
clear
genotype = {'KO', 'WT'}
g = 2;
switch genotype{g}
    case 'KO'
        % These are the PTCHD1 KO data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100720\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC28\101420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC30\101920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\paired_pulse_noisebursts\'};
        
    case 'WT'
        % These are the PTCHD1 WT data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100820\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC29\101520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\paired_pulse_noisebursts\';
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\paired_pulse_noisebursts'};
end
summary_data = [];
decoding_results =[];
for i = 1:length(files)
    cd(files{i})
    load('summary.mat')
    
    for j = 1:length(spikedata)
        if spikedata(j).resp ==1
            fprintf('Processing neuron # %d \n', j)
            for z = 1:10 % 10 repetition of decoding
                options.binsize = 1;
                temp = single_neuron_decoding_PSTH_v2(summary(j), options);
                decoding_run(z) = temp;
                clear temp
            end
            decoding_sum.confMatrix = reshape([decoding_run.confMatrix], size([decoding_run.confMatrix],1), size([decoding_run.confMatrix],1),[]);
            decoding_sum.performance = reshape([decoding_run.performance],size([decoding_run.confMatrix],1),[]);
            decoding_sum.run = decoding_run;
            summary_data = [summary_data, decoding_sum];
            clear decoding_sum 
        end
    end
end

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse')
save(['paired_pulse_results_MGB',genotype{g},'.mat'], 'summary_data', 'files')
clearvars -except files g genotype

summary_data = [];
for i = 1:length(files)
    cd(files{i})
    load('summary_decoding_stats.mat')
    %         for j = 1:length(decoding_results)
    %             decoding_results(j).files = files{i};
    %
    %         end
    iter = 1;
    
    for j = 1:length(spikedata)
        if spikedata(j).resp == 1 && spikedata(j).decoding.results.all.pvalue_stimulus(end) < 0.05
            spikedata(j).resp_decoding = 1;
            decoding_results(iter).performance = diag(spikedata(j).decoding.results.all.data.confusion);
            decoding_results(iter).files = files{i};
            iter = iter + 1;
        end
    end
    save('summary_decoding_stats.mat', 'spikedata')
    load('summary.mat')
    load('summary_decoding_stats.mat')
    clear decoding_results
    decoding_results = [];
    for j = 1:length(spikedata)
        
        if spikedata(j).resp_decoding == 1
            fprintf('Processing neuron # %d \n', j)
            options.binsize = 1;
            for z = 1:10 % 10 repetition of decoding
                temp = single_neuron_decoding_PSTH_v2(summary(j), options);
                decoding_run(z) = temp;
                clear temp
            end
            decoding_sum.confMatrix = reshape([decoding_run.confMatrix], size([decoding_run.confMatrix],1), size([decoding_run.confMatrix],1),[]);
            decoding_sum.performance = reshape([decoding_run.performance],size([decoding_run.confMatrix],1),[]);
            decoding_sum.run = decoding_run;
            decoding_results = [decoding_results, decoding_sum];
            clear decoding_sum
        end
    end
    summary_data = [summary_data, decoding_results];
end

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse')
save(['paired_pulse_results_MGB', genotype{g}, '_stats.mat'], 'summary_data', 'files')

%% make some plots
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse')
genotype = {'ko', 'wt'}
for j = 1:length(genotype)
    switch genotype{j}
        case 'ko'
%                     load('paired_pulse_results_MGBKO.mat')
            load('paired_pulse_results_MGBKO_stats.mat')
        case 'wt'
            load('paired_pulse_results_MGBWT_stats.mat')
%                      load('paired_pulse_results_MGBWT.mat')
            
    end
    
    performance =[];
    for i = 1: length(summary_data)
        summary.(genotype{j}).performance(i,:) = mean(summary_data(i).performance');
    end
end
% save('paired_pulse_summary_MGB.mat', 'summary')
save('paired_pulse_summary_MGB_stats.mat', 'summary')
%%
% figure;
load('paired_pulse_summary_MGB.mat')
% load('paired_pulse_summary_MGB_stats.mat')
figure;
line_errorbar_drc(summary.wt.performance, summary.ko.performance)
xticks([2:2:20])
xticklabels({'1', '2', '4', '8', '16', '32', '64', '128', '256', '512'})
xlabel('Paired Pulse Interval (ms)')
ylabel('Performance')
set(gcf, 'Color', 'w')
% export_fig('Paired_pulse_MGB_stats',  '-png')
export_fig('Paired_pulse_MGB',  '-png')

%% calculate spontaneous activity
clear
genotype = {'KO', 'WT'}
g = 2;
switch genotype{g}
    case 'KO'
        % These are the PTCHD1 KO data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC24\100720\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC28\101420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC30\101920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC31\102920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC37\121520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\paired_pulse_noisebursts\'};
        
    case 'WT'
        % These are the PTCHD1 WT data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC23\100220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC27\100820\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC29\101520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC32\110320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC33\110620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC35\111420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\paired_pulse_noisebursts\';
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\paired_pulse_noisebursts'};
end

summary_data = [];
for i = 1:length(files)
    cd(files{i})
    load('summary_decoding_stats.mat')
    %         for j = 1:length(decoding_results)
    %             decoding_results(j).files = files{i};
    %
    %         end
    iter = 1;
    spont_rate =[];
    for j = 1:length(spikedata)
        if spikedata(j).resp == 1 
            spont = spikedata(j).clusterData.psth.scmatrix(:, 1:250);
            spont_rate(iter) = mean(sum(spont,2))/0.25;
            iter = iter + 1;
        end
    end
    
    summary_data = [summary_data, spont_rate];
end

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse')
save(['paired_pulse_results_MGB_', genotype{g}, '_spont.mat'], 'summary_data', 'files')
%%
load('paired_pulse_results_MGB_WT_spont.mat')
summary_data_wt = summary_data;
load('paired_pulse_results_MGB_KO_spont.mat')
summary_data_ko = summary_data;
[f_wt, x_wt] = ecdf(summary_data_wt);
[f_ko, x_ko] = ecdf(summary_data_ko);
figure;
h1 = plot(x_wt, f_wt, '-k', 'LineWidth', 1);
hold on
h2 = plot(x_ko, f_ko, '-r', 'LineWidth', 1);
xlabel('Spontaneous firing rate (Hz)')
ylabel('Cumulative probablity')
legend([h1, h2], {'WT', 'KO'})
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,500,500])

set(gcf, 'Color', 'w')
export_fig('Spontanous_firing_ecdf',  '-png')

a=summary_data_wt';
b = summary_data_ko';
b(end+1:length(a)) =  NaN;
figure
[h1,h2]=bar_plot([a,b])
h1.FaceColor = 'k';
h2.FaceColor = 'r';

xticks([1, 2])
xticklabels({'WT', 'KO'})
ylabel('Spontaneous Firing Rate (Hz)')
set(gcf, 'Color', 'w')
export_fig('Spontanous_firing_bar_MGB',  '-png')

[h,p] = ttest2(summary_data_ko,summary_data_wt)
%% MGB analysis

%% Let's analyze the Rate Level function for single units
%% Load dataset
clear
path = pwd;
RLF_analysis_preprocess(path, 'MGB')
%% summarize WT the data
wt_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121020\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC36\121120\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122120\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC38\122220\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122720\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC40\122820\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC43\021821\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC43\021921\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC44\031721\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC44\031821\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC46\032421\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC46\032521_1\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\KeC46\032521_2\RateLevel'};
data_wt = [];
for i = 1:length(wt_file)
    cd(wt_file{i})
    load('summary_rlf_responses_v2.mat')
    badUnits = [];
    for j = 1:length(spikedata)
        if isempty(spikedata(j).resp)
            badUnits = [badUnits, j];
        end
    end
    spikedata(badUnits) = [];
    data_wt = [data_wt, spikedata];
end

%% Sumarize the KO data
ko_file = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122420\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC39\122520\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123020\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC41\123120\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC45\031721\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC45\031821\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC49\040921\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC49\041021\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC50\040921\RateLevel';...
    'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\KeC50\041021\RateLevel'};
data_ko = [];
for i = 1:length(ko_file)
    cd(ko_file{i})
    load('summary_rlf_responses_v2.mat')
    badUnits = [];
    for j = 1:length(spikedata)
        if isempty(spikedata(j).resp)
            badUnits = [badUnits, j];
        end
    end
    spikedata(badUnits) = [];
    data_ko = [data_ko, spikedata];
end
%%
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\RLF_analysis')
save('summaryData_raw_MGB.mat', 'data_wt', 'data_ko', 'ko_file', 'wt_file')
