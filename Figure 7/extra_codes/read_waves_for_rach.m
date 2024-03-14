function Y = read_waves(blockName, channelNum, trialNum, tankName, TT)
% READ_WAVES -- read continuous waveform from TDT data tank.
% Usage: Y = read_waves(blockNum, channelNum, trialNum, tankName)
%
% Reads waveform for one channel for a range of trials (i.e. stimulus
% repetitions).
%
% blockNum   : integer block number
% channelNum : channel number
% trialNum   : can be a vector [trial1 trial2] to extract waveforms for a
%              range of trials. Otherwise, scalar trial number.
% tankName   : optional, defaults to 'DemoTank2'
%

if nargin < 4, tankName = 'Tank2'; end

% TT = actxserver('TTank.X'); % Tank
% if TT.ConnectServer('Local', 'Impale') == 0,
%     error('Cannot connect to server.');
% end
if TT.OpenTank(tankName, 'R') == 0,
    error('Cannot open tank.');
    return;
end

% blockName = sprintf('Block-%d', blockNum);
if TT.SelectBlock(blockName) == 0,
    error('Cannot open block.');
end
    
% Find the trial markers (Ticks) to get the trial start times 
maxRead = 5000000;
numRead = maxRead;
numNewEvents = 0;
while numRead == maxRead,
    numRead = TT.ReadEventsV(maxRead, 'ImTr',0,0,0,0,'NEW');
    numNewEvents = numNewEvents + numRead;
end


ev = TT.ParseEvInfoV(0, numNewEvents, 6);

% Set global variables
TT.SetGlobalV('Channel', channelNum);

% Set waveform start and stop times corresponding to user-specified trials
TT.SetGlobalV('T1', ev(trialNum(1)));
TT.SetGlobalV('T2', ev(trialNum(end)+1));

% Retrieve the data
% y = TT.ReadWavesV('Wav1');
 y = TT.ReadWavesV('RSn1'); % add by Ke

% TT.ReleaseServer;

if nargout,
   Y = y;
   return;
end

% figure;
% plot(y);


%--------------------------------------------------------------------------
% END OF READ_WAVES.M
%--------------------------------------------------------------------------


