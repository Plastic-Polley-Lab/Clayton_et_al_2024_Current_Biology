%% extract waveform from Rach
%penNum - penetration number
%dir - directory into which you want to save the raw waveforms
%chans - read 1:128 with RSn1 and 128+(1:64) with RSn2
clearvars;
tankName = 'F:\KeChen\RawData\Rach_Ephys\KeC030920';
penNum = 8;
blockNums = [2:4];
chans = 1:16;
whichdigitizer = 'RSn1';
% chans = 128+(1:64);
% whichdigitizer = 'RSn2';
direc = 'F:\KeChen\RawData\Rach_Ephys\KeC030920\Raw_Waveform-1';
extractWaveforms_Rach(tankName,penNum,blockNums,chans,direc,whichdigitizer)
%% Offline sorting with kilosort
chanMapPath = 'F:\KeChen\MATLAB\chanMap32x2.mat';
% chanMapPath = 'F:\KeChen\MATLAB\chanMap16x1.mat'

%     chanMapPath = 'F:\KeChen\MATLAB\chanMap64.mat'
penetrationNumber = 8
rawDataPrefix = 'KeC030920-8'
rawDataPath = 'F:\KeChen\RawData\Rach_Ephys\KeC030920';
% goodchannel = chanMap(1:32);
NChan = 16;
% Chans = goodchannel;
Chans = 1:16;
region = 'MGB';
blockNums = [2:4];
process_pipeline_Rach(rawDataPath, rawDataPrefix, penetrationNumber, chanMapPath, NChan, Chans, region,blockNums )

%% assemble the clusters
analdir = 'F:\KeChen\RawData\Rach_Ephys\KeC022820\Sorted_Data\TRN-2\Penetration_1'
kilosort_spiketimes_to_clusterData(analdir,chanMapPath)
