function extractWaveforms_Rach(tankName,penNum,blockNums,chans,direc,whichdigitizer)
%penNum - penetration number
%dir - directory into which you want to save the raw waveforms
%chans - read 1:128 with RSn1 and 128+(1:64) with RSn2
% clearvars;
% tankName = 'F:\EPHYS\Rach_Tanks\MMA012420';
% penNum = 2;
% blockNums = 1:23;
% chans = 1:128;
% whichdigitizer = 'RSn1';
% % chans = 128+(1:64);
% % whichdigitizer = 'RSn2';
% direc = 'F:\EPHYS\Passive\MMA012420\Raw_Waveform';

%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Remember to remove the (1)etc. from the blocknames
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tankNum = str2double(tankName(end-5:end));%6 digit dates
if exist(direc, 'dir')
    disp('Raw_Waveform directory already exists');
else
    mkdir(direc);
end

%Do all this stuff outside of primary loop, so as not grab multiple handles
%(seems to generate an OpenEx error)
TT = actxserver('TTank.X'); % Tank

if TT.ConnectServer('Local','Me') == 0 %test if can connect to TDT server
    error('Cannot connect to server.');
end
if TT.OpenTank(tankName, 'R') == 0 %test if can open tank ('R' means read only)
    error('Cannot open tank');
end
file = split(tankName, '\');
file = char(file(end));
for ii = 1:length(blockNums)  
    
    TT.ResetGlobals;
    
    blockName = sprintf('%s-%d-%d', file, penNum, blockNums(ii));%ensure the 6 digit date is zero padded
%     blockName = sprintf('CDQ%06d-%d-%d', tankNum, penNum, blockNums(ii));%ensure the 6 digit date is zero padded
%     blockName = sprintf('rlf_AT-%d-%d', penNum, blockNums(ii));%ensure the 6 digit date is zero padded
    TT.SelectBlock(blockName);
    TT.SetGlobals('MaxReturn=10000000');
    TT.SetGlobalV('WavesMemLimit', 1024^3);
    
    %grab the ticks
    numRead_ticks = TT.ReadEventsSimple('ImTr');
    ticks = TT.ParseEvInfoV(0, numRead_ticks, 6);
    
    %additional stuff to save
    numRead_spikes = TT.ReadEventsSimple('eNeu');
    
    spikes = TT.ParseEvInfoV(0, numRead_spikes, 6);
    sortcodes = TT.ParseEvInfoV(0, numRead_spikes,5);
    channels = TT.ParseEvInfoV(0, numRead_spikes, 4);
    snippets = TT.ParseEvV(0, numRead_spikes);
    
    results = [];
    results.info.tankName = tankName;
    results.info.blockNum = blockNums(ii);
    results.ticks = ticks;
    results.spikes = spikes;
    results.sortcodes = sortcodes;
    results.channels = channels;
    results.snippets = snippets;
    
    savename = [direc sprintf('/data_tank_%s_pen_%d_block_%d.mat',file,penNum,blockNums(ii))];
    save(savename,'results');
    
    disp(sprintf('Extracting data in tank %s: pen %d: block %d',file,penNum,blockNums(ii)));
    
%     data_raw = SEV2mat([tankName, '\', blockName]); % add by Ke, load all data simultaneously, very slow

    for i=1:length(chans)
        disp(sprintf('...extracting and saving channel %d',chans(i)));
        
        %initialize the results struct that will be saved later
        results=[];
%         data = data_raw.(whichdigitizer).data(i,:); % add by Ke

        %read the i^th channel
        TT.SetGlobalV('Channel', i); % edit by Ke 
        data = TT.ReadWavesV(whichdigitizer); % commented by Ke
%           data_raw = SEV2mat([tankName, '\', blockName], 'CHANNEL', i, 'DEVICE', whichdigitizer, 'TANK', tankName, 'BLOCK',blockName  ); % add by Ke
%           data = data_raw.(whichdigitizer).data; % add by Ke
        if isnan(data)
            tdt_struct = TDT2mat(sprintf('%s',tankName),sprintf('Pen-%d',penNum),sprintf('Block-%d',blockNums(ii)), 'STORE', 'Wave', 'CHANNEL', i);       
            if ~isempty(tdt_struct.streams)
                data = tdt_struct.streams.Wave.data';
            end
        end
        
        %create the data struct to save
        results.info.tankName = tankName;
        results.info.blockNum = blockNums(ii);
        results.info.channel = chans(i);
        
        %at this point, check again
        if isnan(data)
            disp('DATA RETURNED FROM BLOCK IS NAN, saving details and moving on...');
        end
        results.waveform = data'; % edit by Ke
%         results.waveform = data;
        savename = [direc sprintf('/waveform_tank_%s_pen_%d_block_%d_channel_%d.mat',file,penNum,blockNums(ii),chans(i))];
        save(savename,'results');
        clear data
    end
    clear data_raw
end

%close tank and release server
TT.CloseTank;
TT.ReleaseServer;


