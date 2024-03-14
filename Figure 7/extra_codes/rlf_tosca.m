%% Preprocess
% Step 1: assemble the tosca files with the spikes
path = {'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session2\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session3\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session4\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session2\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session3\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session4\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session5\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session6\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session2\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session3\ACtx-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session4\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session5\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC44\KeC44-Session5\MGB-Raw-Waveform';...
    'Z:\KeChen\Tosca Recording\KeC44\KeC44-Session4\MGB-Raw-Waveform';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC44\ACtx-Session3';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC44\ACtx-Session2';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43\ACtx_Session1\ACtx-Raw-Waveform';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43\ACtx_Session2';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43\MGB_Session3';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43\MGB-Session4'};

for i = 1:length(path)
    kilosort_spiketimes_Tosca_clusterData(path{i})
end

%% Preprocess
% Step2: re-organize the data structure
save_path = {'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session2.mat';...
    'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session3.mat';...
    'Z:\KeChen\Tosca Recording\KeC47\KeC47-Session4.mat';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session2.mat';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session3.mat';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session4.mat';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session5.mat';...
    'Z:\KeChen\Tosca Recording\KeC46\KeC46-Session6.mat';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session2.mat';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session3.mat';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session4.mat';...
    'Z:\KeChen\Tosca Recording\KeC45\KeC45-Session5.mat';...
    'Z:\KeChen\Tosca Recording\KeC44\KeC44-Session5.mat';...
    'Z:\KeChen\Tosca Recording\KeC44\KeC44-Session4.mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC44-Session3,mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC44-Session2.mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43-Session1.mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43-Session2.mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43-Session3.mat';...
    'F:\KeChen\RawData\Tosca_Ephys\KeC43-Session4.mat';...
    };

for i = 1:length(path)
    cd([path{i}, '_SpikeData'])
    files = dir('*.mat')
    for j = 1:length(files)
        fprintf('Processing Session %d Neuron %d\n', i, j)
        load(files(j).name)
        clusterData.run.name = files(j).name;
        spikedata(j) = clusterData.run;
        spikedata(j).stimData.trial =[]; % reduce the file size
    end
    save(save_path{i}, 'spikedata')
    clear spikedata
end
%% preprocess
% Step 3: raster the data
Tosca_raster


%% Preprocess
% Step 4: Manually check the spike quality
% check the GUI of the analysis_toolbox
% Let's analyze the Rate Level function for single units
%% Check the responses
clear; close all; clc
path = pwd;
RLF_analysis_preprocess_tosca(path, 'ACtx')
%% Load dataset
clear
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\ACtx\Tosca\rasterData')
% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-WT\MGB\Tosca\')

wt_file = dir('*_rlf_responses.mat');
data_wt =[];
for i = 1:length(wt_file)
    load(wt_file(i).name)
    badUnits = [];
    for j = 1:length(spikedata)
        if isempty(spikedata(j).resp)
            badUnits = [badUnits, j];
        end
    end
    spikedata(badUnits) = [];
    data_wt = [data_wt, spikedata];
end


cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\ACtx\Tosca\rasterData')
% cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\PTCHD1-KO\MGB\Tosca')

ko_file = dir('*_rlf_responses.mat');

data_ko = [];
for i = 1:length(ko_file)
    load(ko_file(i).name)
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
cd('E:\Ke_Chen\Processed Data\PTCHD1-Project\Tosca\RLF_analysis')
save('summaryData_raw.mat', 'data_wt', 'data_ko', 'ko_file', 'wt_file')
% save('summaryData_raw_MGB.mat', 'data_wt', 'data_ko', 'ko_file', 'wt_file')

%% Let's analyze excitatory neurons
clearvars -except data_wt data_ko
% load('summary_noise_ACtx_normalized.mat')
data_ko([data_ko.resp]==0) = [];
data_wt([data_wt.resp]==0) = [];

data_ko([data_ko.refra_violation_ratio]>0.005) =[];
data_wt([data_wt.refra_violation_ratio]>0.005) =[];

data = [data_ko, data_wt];
ko_indx = 1:length(data_ko);
wt_indx = length(data_ko)+1 : length(data);
%% summarize the data; let's first only look at excitatory responses
sign = 1; % extract exciation 
cluster1 = find([data.resp]==sign);
wt_cluster1 = intersect(cluster1, wt_indx);
ko_cluster1 = intersect(cluster1, ko_indx);
rlf =[];
%%
results.wt.exc = gain_threshold_rlf_tosca(data(wt_cluster1));
results.ko.exc = gain_threshold_rlf_tosca(data(ko_cluster1));
%%
figure
[fig1, fig2] = ecdf_bar_plot([results.wt.exc.MI], [results.ko.exc.MI], 'Monotonic Index');
set(get(fig1,'XLabel'), 'String', 'Monotonic Index');
set(get(fig2,'YLabel'), 'String', 'Monotonic Index');

%% organize the data for plots
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    for i = 1:length(results.(genotype{j}).exc)
%         plots.(genotype{j}).rlf(i,:) = results.(genotype{j})(i).rlf_avg;
        plots.(genotype{j}).exc.rlf(i,:) = results.(genotype{j}).exc(i).rlf_avg - ...
            results.(genotype{j}).exc(i).spont;

    end
%     plots.(genotype{j}).exc.gains       = [results.(genotype{j}).exc.gains_value];
%     plots.(genotype{j}).exc.threshold   = [results.(genotype{j}).exc.threshold];
    plots.(genotype{j}).exc.latency_p2t = [results.(genotype{j}).exc.latency_p2t];
    plots.(genotype{j}).exc.rlf_p_value = [results.(genotype{j}).exc.rlf_p_value];
    
%     nan_indx = find(isnan(plots.(genotype{j}).exc.threshold));
%     
%     plots.(genotype{j}).exc.gains(nan_indx)=[];
%     plots.(genotype{j}).exc.threshold(nan_indx) =[];
%     plots.(genotype{j}).exc.rlf(nan_indx,:) = [];
%     plots.(genotype{j}).exc.latency_p2t(nan_indx) = [];
end
%%
figure;
line_errorbar_drc(plots.wt.exc.rlf, plots.ko.exc.rlf)
xticks([1:1:7])
xticklabels({ '35', '45', '55', '65', '75', '85', '95'})
xlabel('Sound Intensity(dB)')
ylabel('Firing Rate (Hz)')
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
xlim([0, 8])
export_fig('Rate_Level_Function',  '-png', '-pdf')
