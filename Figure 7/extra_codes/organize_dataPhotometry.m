%% Get data
%TOSCA
function data = organize_dataPhotometry(folder,mouse)
folder_tosca = [folder, mouse '\randreg_tosca.mat'];
data_tosca=load(folder_tosca);
dTosca=data_tosca.randreg_tosca;
%         dToscaTot{iMouse,iFiber}=dTosca;
disp('Tosca data loaded')

%FIBER PHOTOMETRY
folder_fibPhot=[folder, mouse '\randreg_fibPhot.mat'];
data_fibPhot=load(folder_fibPhot);
dFibPhot=data_fibPhot.randreg_fibPhot;
disp('Fiber photometry data loaded')

% Assemble the fiber photometry data with the stimulus data
for i = 1 : length(dTosca.set)
    for j = 1:length(dTosca.set(i).trials)
        dFibPhot.trial(dTosca.set(i).trials(j)).sets = i;
        dFibPhot.trial(dTosca.set(i).trials(j)).cycTab = dTosca.set(i).cycTab(j);
        dFibPhot.trial(dTosca.set(i).trials(j)).tOnset = dTosca.set(i).tOnset(j);
        dFibPhot.trial(dTosca.set(i).trials(j)).tOffset = dTosca.set(i).tOffset(j);
        dFibPhot.trial(dTosca.set(i).trials(j)).type = dTosca.set(i).type;
        dFibPhot.trial(dTosca.set(i).trials(j)).stimuli=getTimeBursts(i,dTosca.set(i).cycTab(j));
    end
end
%%
clearvars -except dFibPhot mouse 
%Build time
for i = 1:length(dFibPhot.trial)
    nSamples=length(dFibPhot.trial(i).fiber1.filtered);
    dFibPhot.trial(i).time=[0:nSamples-1]/dFibPhot.fs;
end
%% Build REG and RAND
for trial = 1:length(dFibPhot.trial)
    cyc = dFibPhot.trial(trial).cycTab;
    stimuli = dFibPhot.trial(trial).stimuli + dFibPhot.trial(trial).tOnset;
    time = dFibPhot.trial(trial).time;
    RAND1_start = stimuli(11); % add by Ke
    REG_start = stimuli(51);
    RAND_start   = stimuli(50+ cyc* 25 + 1);
    idxRAND1_Start = find(time>RAND1_start,1,'first'); % add by Ke
    idxREG_Start=find(time>REG_start,1,'first');
    idxRAND_Start=find(time>RAND_start,1,'first');
    duration = idxRAND_Start- idxREG_Start;
    
    dFibPhot.trial(trial).RAND1.time=time(idxRAND1_Start:idxREG_Start);
    dFibPhot.trial(trial).RAND1.stimuli = stimuli(11:51);
    dFibPhot.trial(trial).RAND1.fiber1 = dFibPhot.trial(trial).fiber1.signal(idxRAND1_Start:idxREG_Start);
    dFibPhot.trial(trial).RAND1.fiber2 = dFibPhot.trial(trial).fiber2.signal(idxRAND1_Start:idxREG_Start);
    
    
    dFibPhot.trial(trial).REG.time=time(idxREG_Start:idxRAND_Start);
    dFibPhot.trial(trial).REG.stimuli = stimuli(51:50+ cyc* 25 + 1);
    dFibPhot.trial(trial).REG.fiber1 = dFibPhot.trial(trial).fiber1.signal(idxREG_Start:idxRAND_Start);
    dFibPhot.trial(trial).REG.fiber2 = dFibPhot.trial(trial).fiber2.signal(idxREG_Start:idxRAND_Start);
    
    dFibPhot.trial(trial).RAND.time=time(idxRAND_Start:duration+idxRAND_Start);
    dFibPhot.trial(trial).RAND.stimuli = stimuli(50+ cyc* 25 + 1:end);
    dFibPhot.trial(trial).RAND.fiber1 = dFibPhot.trial(trial).fiber1.signal(idxRAND_Start:duration+idxRAND_Start);
    dFibPhot.trial(trial).RAND.fiber2 = dFibPhot.trial(trial).fiber2.signal(idxRAND_Start:duration+idxRAND_Start);
end
%%
% fiberTab={'fiber1','fiber2'};
% location = {'Rostral', 'Caudal'};
% 
% sets = 1
% cyc  = 4
% 
% trial = intersect(find([dFibPhot.trial.cycTab] == cyc),find([dFibPhot.trial.sets] == sets))
% figure(1);
% for iFiber = 1:2
%     subplot(2, 1, iFiber)
%     fiber=fiberTab{iFiber};
%     plot(dFibPhot.trial(trial).time, dFibPhot.trial(trial).(fiber).signal)
%     stimuli = dFibPhot.trial(trial).stimuli + dFibPhot.trial(trial).tOnset;
%     hold on
%     scatter(stimuli, max(dFibPhot.trial(trial).fiber1.signal) * ones(size(stimuli)), 1,'o')
%     title([location{iFiber}, ' Cyc-', num2str(dFibPhot.trial(trial).cycTab)])
%     set(gcf,'position',[100,200,1600,600])
% end
%% Calculate the mean value
% for i = 1:length(dFibPhot.trial)
%     dFibPhot.trial(i).REG.avg_fiber1 = mean(dFibPhot.trial(i).REG.fiber1);
%     dFibPhot.trial(i).REG.avg_fiber2 = mean(dFibPhot.trial(i).REG.fiber2);
%     dFibPhot.trial(i).RAND.avg_fiber1 = mean(dFibPhot.trial(i).RAND.fiber1);
%     dFibPhot.trial(i).RAND.avg_fiber2 = mean(dFibPhot.trial(i).RAND.fiber2);
% end
% % color =cbrewer('div', 'RdYlBu',4);
% fiberTab={'fiber1','fiber2'};
% avg_signal.cyc = sort(unique([dFibPhot.trial.cycTab]));
% avg_signal.sets = sort(unique([dFibPhot.trial.sets]));
% for j = 1:length(fiberTab)
%     avg_signal.(fiberTab{j}) = zeros(length(avg_signal.cyc), length(avg_signal.sets));
% %     figure;
%     for i = 1:length(dFibPhot.trial)
%         cyc = dFibPhot.trial(i).cycTab;
%         cyc_idx = find(avg_signal.cyc == cyc);
%         sets = dFibPhot.trial(i).sets;
%         avg_signal.(fiberTab{j})(cyc_idx, sets) = dFibPhot.trial(i).REG.(['avg_',fiberTab{j}])
% %         hold on
% %         scatter(cyc, dFibPhot.trial(i).REG.(['avg_',fiberTab{j}]), 'o', 'filled', 'MarkerFaceColor',color(sets,:))
% %         xlim([2, 14])
%     end
%     
% %     title(fiberTab{j})
% end
data.animal = mouse;
data.dFibPhot = dFibPhot;
% data.avg_signal = avg_signal;
% figure
% for 
% for i = 1:length(dFibPhot.trial)
%     cyc = dFibPhot.trial(i).cycTab;
%     sets = dFibPhot.trial(i).sets
%     hold on
%     scatter(cyc, dFibPhot.trial(i).REG.avg_fiber2, 'o', 'filled', 'MarkerFaceColor',color(sets,:))
%     xlim([2, 14])
% end
%% Calculate the PSD for REG and RAND



