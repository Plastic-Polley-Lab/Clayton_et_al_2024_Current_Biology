function ev = read_ticks(blockName,tankName, TT)


% TT = actxserver('TTank.X'); % Tank
% if TT.ConnectServer('Local', 'Impale') == 0,
%     error('Cannot connect to server.');
% end
% if TT.OpenTank(tankName, 'R') == 0,
if TT.OpenTank(tankName, 'R') == 0,
    error('Cannot open tank.');
    return;
end
% blockName = sprintf('Block-%d', blockNum);
if TT.SelectBlock(blockName) == 0,
    error('Cannot open block.');
    ev = 0;
    return;
end
maxRead = 5000000;
numRead = maxRead;
numNewEvents = 0;
while numRead == maxRead,
    numRead = TT.ReadEventsV(maxRead, 'ImTr',0,0,0,0,'NEW');
    numNewEvents = numNewEvents + numRead;
end
ev = TT.ParseEvInfoV(0, numNewEvents, 6);

