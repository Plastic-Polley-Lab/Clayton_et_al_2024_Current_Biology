%%Step1: Pre-Process the data
clear
pair_pulse_PreProcess_Auto(50, 50, 0.01)
%% decoding with statistic inferences
load('summary.mat')
addpath('E:\Ke_Chen\MATLAB\Ephys-analysis\NeuralDecoder')
for i = 1: length(spikedata)
    if spikedata(i).keep ==1 && spikedata(i).resp == 1
        fprintf('Processing neuron # %d of total %d\n', i, length(spikedata))
        data = summary(i);
        try
        [results, options] = SU_decoding_PSTH(data);
        spikedata(i).decoding.results = results;
        spikedata(i).decoding.options = options;
        catch
            continue
        end
    end
end
save('summary_decoding_stats.mat', 'spikedata', 'summary', 'decoding_results')
%% summary the results across different sessions
clear
genotype = {'KO', 'WT'}
g = 1;
switch genotype{g}
    case 'KO'
        % These are the PTCHD1 KO data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091720\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092120\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092720\paired_pulse_noisebursts\'; ...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092820\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\paired_pulse_noisebursts'
            };
        
    case 'WT'
        % These are the PTCHD1 WT data
        files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091820\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\092020\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092420\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100720\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110620\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\paired_pulse_noisebursts\';...
            'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\paired_pulse_noisebursts\'
            };
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

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
save(['paired_pulse_results_',genotype{g},'.mat'], 'summary_data', 'files')
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

cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
save(['paired_pulse_results_', genotype{g}, '_stats.mat'], 'summary_data', 'files')
%% re-summarize data
summary_data = [];
for i = 1:length(files)
    cd(files{i})
    load('summary.mat')
    load('summary_decoding_stats.mat')
    clear decoding_results
    for j = 1:length(spikedata)
        
        if spikedata(j).resp == 1
            fprintf('Processing neuron # %d \n', j)
            options.binsize = 1;
            for z = 1:10 % 10 repetition of decoding
                temp = single_neuron_decoding_PSTH_v2(summary(j), options);
                decoding_run(z) = temp;
                clear temp
            end
            spikedata(j).confMatrix = reshape([decoding_run.confMatrix], size([decoding_run.confMatrix],1), size([decoding_run.confMatrix],1),[]);
            spikedata(j).performance = reshape([decoding_run.performance],size([decoding_run.confMatrix],1),[]);
            spikedata(j).run = decoding_run;
%             decoding_results = [decoding_results, decoding_sum];
            clear decoding_sum
        end
    end
    summary_data = [summary_data, spikedata];
%     summary_data = [summary_data, decoding_results];
end
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summarize the data into a small file for plots
clear
path = 'E:\Ke_Chen\Processed Data\PTCHD1-Project\Paired_pulse';
summaryData = summary_paired_pulse_analysis(path, 'MGB');

plots.wt.performance = reshape([summaryData.wt.data.performance_avg], 20, [])';
plots.ko.performance = reshape([summaryData.ko.data.performance_avg], 20, [])';
plots.wt.performance_stats = reshape([summaryData.wt.data_stats.performance_avg], 20, [])';
plots.ko.performance_stats = reshape([summaryData.ko.data_stats.performance_avg], 20, [])';

%% make some plots from the summary data

load('summary_data_all.mat')
% load('summary_data_all_MGB.mat.mat')
figure;
line_errorbar_drc(plots.wt.performance, plots.ko.performance)
% line_errorbar_drc(plots.wt.performance_stats, plots.ko.performance_stats)
xticks([2:2:20])
xticklabels({'1', '2', '4', '8', '16', '32', '64', '128', '256', '512'})
xlim([0, 20])
xlabel('Paired Pulse Interval (ms)')
ylabel('Performance')
set(gcf, 'Color', 'w')
% export_fig('Paired_pulse',  '-png')

% plot the spontanous firing rate
figure
[fig1, fig2] = ecdf_bar_plot([summaryData.wt.data.spont], [summaryData.wt.data.spont], summaryData.region)
set(get(fig1,'XLabel'), 'String', 'Spontaneous FR (HZ)');
set(get(fig2,'YLabel'), 'String', 'Spontaneous FR (HZ)');

%% calculate spontaneous activity
% clear
% genotype = {'KO', 'WT'}
% g = 2;
% switch genotype{g}
%     case 'KO'
%         % These are the PTCHD1 KO data
%         files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC22\091720\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092120\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092220\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC24\092320\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092620\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092720\paired_pulse_noisebursts\'; ...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Kec26\092820\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100820\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC28\100920\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC30\101920\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC31\102920\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\KeC34\111120\paired_pulse_noisebursts'
%             };
%         
%     case 'WT'
%         % These are the PTCHD1 WT data
%         files = {'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091820\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\091920\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC23\092020\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092420\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC25\092520\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100520\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100620\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC27\100720\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC29\101520\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC32\110220\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC33\110620\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111320\paired_pulse_noisebursts\';...
%             'E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\KeC35\111420\paired_pulse_noisebursts\'
%             };
% end
% 
% summary_data = [];
% for i = 1:length(files)
%     cd(files{i})
%     load('summary_decoding_stats.mat')
%     %         for j = 1:length(decoding_results)
%     %             decoding_results(j).files = files{i};
%     %
%     %         end
%     iter = 1;
%     spont_rate =[];
%     for j = 1:length(spikedata)
%         if spikedata(j).resp == 1 
%             spont = spikedata(j).clusterData.psth.scmatrix(:, 1:250);
%             spont_rate(iter) = mean(sum(spont,2))/0.25;
%             iter = iter + 1;
%         end
%     end
%     
%     summary_data = [summary_data, spont_rate];
% end

% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project')
% save(['paired_pulse_results_', genotype{g}, '_spont.mat'], 'summary_data', 'files')
%%
% load('paired_pulse_results_WT_spont.mat')
% summary_data_wt = summary_data;
% load('paired_pulse_results_KO_spont.mat')
% summary_data_ko = summary_data;
% [f_wt, x_wt] = ecdf(summary_data_wt);
% [f_ko, x_ko] = ecdf(summary_data_ko);
% figure;
% h1 = plot(x_wt, f_wt, '-k');
% hold on
% h2 = plot(x_ko, f_ko, '-r');
% xlabel('Spontaneous firing rate (Hz)')
% ylabel('Cumulative probablity')
% legend([h1, h2], {'WT', 'KO'})
% box off
% set(gca,'TickDir','out')
% set(gca,'fontsize',12)
% set(gca,'TickLengt', [0.015 0.015]);
% set(gca, 'LineWidth',1)
% set(gcf,'position',[100,200,500,500])
% 
% set(gcf, 'Color', 'w')
% export_fig('Spontanous_firing_ecdf',  '-png')
% 
% a=summary_data_wt';
% b = summary_data_ko';
% b(end+1:length(a)) =  NaN;
% figure
% bar_plot([a,b])
% set(gcf, 'Color', 'w')
% export_fig('Spontanous_firing_bar',  '-png')