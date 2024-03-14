function process_pipeline_Rach(rawDataPath, rawDataPrefix, penetrationNumber, chanMapPath, NChan, Chans, region,blockNums, raw_waveform_path, manual )
    %process_pipeline: Concatenate, sort, and pre-cluster arbitrary high-
    %                  density neural recording data.
    %Usage: 
    %       STATUS = process_pipeline(rawDataPath, rawDataPrefix, 
    %                                 penetrationNumber, processMode)
    %Options:
    %       rawDataPath       : The fully-qualified path to the directory
    %                           containing the unconcatenated MAT files of 
    %                           recordings
    %       rawDataPrefix     : A string prefix that is invariant among all
    %                           recordings
    %       penetrationNumber : Which penetration this is (per mouse)
    %       chanMapPath       : Which channel mapping to load
    
%     for eg. 
% 
%     clearvars;
%     chanMapPath = 'F:\KeChen\MATLAB\chanMap32x2.mat'
% %     chanMapPath = 'F:\KeChen\MATLAB\chanMap64.mat'
%     penetrationNumber = 2
%     rawDataPrefix = 'MMA012420-2'
%     rawDataPath = 'F:\EPHYS\Passive\MMA012420'
%     NChan = 64;
%     Chans = 64+(1:64);
%     region = 'MGB';
%     blockNums = 1:23;

%     clearvars;
%     chanMapPath = 'F:/EPHYS/Meenakshi_Code/Channel_Maps/chanMap_2trode_16each.mat'
%     penetrationNumber = 1
%     rawDataPrefix = 'MMA080218-1'
%     rawDataPath = 'F:\EPHYS\Practice_Awake\ACTXandMGB_Mouse3\MMA080218'
%     NChan = 32;

%     clearvars;
%     chanMapPath = 'F:/EPHYS/Meenakshi_Code/Channel_Maps/chanMap16_forSteve.mat'
%     penetrationNumber = 1
%     rawDataPrefix = 'MMA112018-2'
%     rawDataPath = 'F:\EPHYS\Awake_C57s\MGB_only\MMA112018'
%     NChan = 16;
    %% concatenate raw data files into one file to be sorted. 

    fpath = [rawDataPath '/Impale'];
    to_concat = findFilesToConcatenate(fpath, rawDataPrefix);
    to_concat = to_concat(blockNums); % add by Ke to only analyze data for some blocks
    concatenate_and_preprocess_Rach(rawDataPath, to_concat, penetrationNumber, rawDataPrefix,NChan,Chans,region,blockNums,raw_waveform_path);


    %% Do spike sorting on concatenated data
%     if manual ==0
%     penDir = fullfile(rawDataPath,sprintf('/Sorted_Data/%s/Penetration_%d/',region,penetrationNumber));   
% %     run_kilosort(penDir, chanMapPath, NChan);
%     run_master_kilosort2(penDir, chanMapPath, NChan) 
%     else
%         kilosort % add by Ke for manually run kilosort
%     end
end




