% RAND-REG ANALYSIS SCRIPT   ---   RANDREG_ANALYZE_MIEC
% preliminary analysis of REG-RAND data, specifically looking at
% differences in autocorrelation and FFT between the REG and RAND epochs

% Blaise Robert, 04/28/2020; modified by Ke, 04/29/2020
%% Step 1:
% process the raw data and sort the behavioral events; there is for
% jittered data
clear
mouseTab={'DFI045','DFI047','DFI048'};
for i = 1:length(mouseTab)
    animalID = mouseTab{i};
    path = ['Z:\File Transfer\Blaise\RANDREG_forKe\Saved', '\', animalID, '\', 'Session 6'];
    cd(path)
    runFile = [animalID, '-Session6-Run2.txt']
    struct_tosca=sortFilesRandRegKe(runFile);
    cd(['Z:\File Transfer\Blaise\RANDREG_forKe\Saved', '\', animalID])
    save('randreg2_tosca.mat', 'struct_tosca')
end

%%
clear all
close all

%% Initialization
mouseTab={'DFI041','DFI042','DFI043','DFI044','DFI045','DFI047','DFI048'};
folder='Z:\File Transfer\Blaise\RANDREG_forKe\Saved\' ;
cd(folder)
mkdir('Pre-Processed')
cd('Pre-Processed')
%% for normal data
for i = 1: length(mouseTab)
    mouse = mouseTab{i};
    data = organize_dataPhotometry(folder,mouse);
    save([mouse, '_day1.mat'], 'data')
end
%% for jittered data
mouseTab={'DFI045','DFI047','DFI048'};
for i = 1: length(mouseTab)
    mouse = mouseTab{i};
    data = organize_dataPhotometry_v2(folder,mouse);
    save([mouse, '_jitter.mat'], 'data')
end
% Initialization for 2nd day
% mouseTab={'DFI043' 'DFI044'};
% fiberTab={'fiber1','fiber2'};
% folder='Z:\File Transfer\Blaise\RANDREG_forKe\Saved\' ;
% for i = 1: length(mouseTab)
%     mouse = mouseTab{i};
%     data(i) = organize_dataPhotometry(folder,mouse);
% end
%% calculate and plot the average
data = avg_signal(data);

% color =cbrewer('div', 'RdYlBu',4);
% for i = 1:6
%     figure
%     for j = 1 :length(data(i).REG.avg_signal.sets)
%         h{j} = scatter(data(i).REG.avg_signal.cyc, data(i).REG.avg_signal.fiber2(:,j) - data(i).RAND.avg_signal.fiber2(:,j) , 'o', 'filled', 'MarkerFaceColor',color(j,:))
%         hold on
%     end
%     xlim([2,14])
%     legend([h{1}, h{2}, h{3}, h{4}],{'Set1', 'Set2', 'Set3', 'Rand'})
% end
% xlim([2,14])
% legend([h{1}, h{2}, h{3}, h{4}],{'Set1', 'Set2', 'Set3', 'Rand'})
% xlabel('Cycle Size')
% ylabel('Avaerge dF/F')
%% Calculate the power spectral density and plot
for i = 1:length(data)
    temp = data(i);
    temp_new(i) = psd_signal(temp);
end
data = temp_new;
clear temp_new
%% let's see autocorrelation
for i = 1:length(data)
    temp = data(i);
    temp_new(i) = autoCor_signal(temp);
end
data = temp_new;
clear temp_new
%% let's see the cross-correlation between stimuli and signal
for i = 1:length(data)
    temp = data(i);
    temp_new(i) = cross_signal(temp);
end
data = temp_new;
clear temp_new
%% sumamry across mice
analyze_sig = 'cross_signal';
type = {'REG', 'RAND'}
fiber = {'fiber1', 'fiber2'};
fiberid = 1
clear data_ana
for i = 1:6
    for kk = 1:2
        data_ana(i).(type{kk}).(analyze_sig) = zeros(size(data(i).(type{kk}).(analyze_sig).(fiber{fiberid})));
        for j = 1:size(data_ana(i).(type{kk}).(analyze_sig), 1)
            for z = 1:size(data_ana(i).(type{kk}).(analyze_sig),2)
                if strcmp(analyze_sig, 'avg_signal')
                data_ana(i).(type{kk}).(analyze_sig)(j,z) = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid})(j, z));
                else
                data_ana(i).(type{kk}).(analyze_sig)(j,z) = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid}){j, z});
                end
                
            end
        end
    end
end
%%
figure(100)
across_sets_plot(data_ana, fiberid)
%% plot
figure(1);
subplot(1,2,1)
location = {'NDB', 'SI/GP'}
cyc = 12
color =cbrewer('div', 'RdYlBu',4);

for i = 1:6
    switch cyc
        case 4
            idx = 1
            cyclesize = '4'
        case 6
            idx = 2
            cyclesize = '6'
        case 8
            idx = 3
            cyclesize = '8'
        case 12
            idx = 5
            cyclesize = '12'
    end
    reg(i,:) = data_ana(i).REG.(analyze_sig)(idx,:);
    rand_d(i,:) = data_ana(i).RAND.(analyze_sig)(idx,:);
    start = 1;
    % add summary across sets
        reg_avg(i,1) = reg(i,1);
        reg_avg(i,2) = reg(i,4);
        rand_avg(i,1)  = rand_d(i,1);
        rand_avg(i,2)  = rand_d(i,4);

%     reg_avg(i,1) = mean(reg(i,1:3));
%     reg_avg(i,2) = mean(reg(i,4));
%     rand_avg(i,1)  = mean(rand_d(i,1:3));
%     rand_avg(i,2)  = mean(rand_d(i,4));
    for j = 1 :size(reg,2)
        h{j} = plot([start, start+1], [reg(i,j), rand_d(i,j)],'-o', 'color' ,color(j,:), 'MarkerFaceColor',color(j,:), 'LineWidth', 1)
        hold on
        start = start +3;
    end
    
end
xlim([0, 12])
ylim([0,0.4])
xticks([1 2 4 5 7 8 10,11])
xticklabels({'reg','rand','reg','rand','reg','rand','rand','rand'})
legend([h{1}, h{2}, h{3}, h{4}],{'Set1','Set2', 'Set3','RAND-Control'})
ylabel('Peak Cross-correlation')
title([location{fiberid}, ' Cycle Size ', cyclesize])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,600])

x = [0.7 2.3 3.7 5.3 6.7 8.3 9.7,11.3];
hold on
for i = 1:size(reg, 2)
    %     scatter(x(i, i+1), [mean(reg(:,i)), mean(rand_d(:,i))], 'd','MarkerFaceColor', color(i,:))
    errorbar(x(2*i-1:i*2),[mean(reg(:,i)), mean(rand_d(:,i))],([std(reg(:,i)), std(rand_d(:,i))])/sqrt(size(reg,1)), 'd','color',[0,0,0], 'MarkerFaceColor', [0,0,0], 'LineWidth', 1, 'CapSize',12)
end
%%
figure(2);
subplot(1,2,2)
color =cbrewer('div', 'RdYlBu',4);
for i = 1:size(reg_avg,1)
        start = 1;
        for j = 1 :size(reg_avg,2)
        h{j} = plot([start, start+1], [reg_avg(i,j), rand_avg(i,j)],'-o', 'color' ,color(j,:), 'MarkerFaceColor',color(j,:), 'LineWidth', 1)
        hold on
        start = start +3;
    end
end
xlim([0,6])
ylim([-0.1,0.4])
xticks([1 2 4 5])
xticklabels({'reg','rand','rand','rand'})
legend([h{1}, h{2}],{'REG-RAND','RAND-Control'})
ylabel('Peak Cross-correlation')
title([location{fiberid}, ' Cycle Size ', cyclesize])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,600])
x = [0.7 2.3 3.7 5.3];
hold on
for i = 1:size(reg, 2)
    %     scatter(x(i, i+1), [mean(reg(:,i)), mean(rand_d(:,i))], 'd','MarkerFaceColor', color(i,:))
    errorbar(x(2*i-1:i*2),[mean(reg_avg(:,i)), mean(rand_avg(:,i))],([std(reg_avg(:,i)), std(rand_avg(:,i))])/sqrt(size(reg_avg,1)), 'd','color',[0,0,0], 'MarkerFaceColor', [0,0,0], 'LineWidth', 1, 'CapSize',12)
end
%% get the dynamic change of cross-correlation;
for i = 1:length(data)
    temp = data(i);
    temp_new(i) = cross_signal_dynamic(temp);
end
data = temp_new;
clear temp_new

%% sumamry across mice
analyze_sig = 'cross_signal';
type = {'REG', 'RAND'}
fiber = {'fiber1', 'fiber2'};
fiberid = 2
clear data_ana
for i = 1:6
    for kk = 1:2
        numFrame = size(data(1).REG.cross_signal.fiber2{1},1);
        data_ana(i).(type{kk}).(analyze_sig) = zeros([size(data(i).(type{kk}).(analyze_sig).(fiber{fiberid})),numFrame]);
        for j = 1:size(data_ana(i).(type{kk}).(analyze_sig), 1)
            for z = 1:size(data_ana(i).(type{kk}).(analyze_sig),2)
                if strcmp(analyze_sig, 'avg_signal')
                data_ana(i).(type{kk}).(analyze_sig)(j,z) = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid})(j, z));
                else
                 range_idx = 800:1200;
                [data_ana(i).(type{kk}).(analyze_sig)(j,z,:), I] = max(data(i).(type{kk}).(analyze_sig).(fiber{fiberid}){j, z}(:,range_idx),[], 2);
                I
                end
                
            end
        end
    end
end
%%
%% plot
figure(1);
cyc = 12
switch cyc
    case 4
        idx = 1;
    case 12
        idx = 2;
end
subplot(1,2,idx)
% subplot(1,2,2)
location = {'NDB', 'SI/GP'}
color =cbrewer('div', 'RdYlBu',4);
clear reg rand_d
for i = 1:6
    switch cyc
        case 4
            idx = 1
            cyclesize = '4'
        case 6
            idx = 2
            cyclesize = '6'
        case 8
            idx = 3
            cyclesize = '8'
        case 12
            idx = 5
            cyclesize = '12'
    end
    reg(i,:,:) = squeeze(data_ana(i).REG.(analyze_sig)(idx,:,:));
    rand_d(i,:,:) = squeeze(data_ana(i).RAND.(analyze_sig)(idx,:,:));
    start = 0;
    % add summary across sets
    reg_avg(i,1,:) = mean(reg(i,1:3,:), 2);
    reg_avg(i,2,:) = mean(reg(i,4,:),2);
    rand_avg(i,1,:)  = mean(rand_d(i,1:3,:), 2);
    rand_avg(i,2,:)  = mean(rand_d(i,4,:),2);
    for j = 1 :size(reg,2)
        dynamic_sig = [squeeze(reg(i,j,:)); squeeze(rand_d(i,j,:))];
        h{j} = plot((1:10)+start, dynamic_sig,'-o', 'color' ,color(j,:), 'MarkerFaceColor',color(j,:), 'LineWidth', 1)
        hold on
        start = start+11;
    end
    
end
% xlim([0, 12])
% ylim([0,0.4])
% xticks([1 2 4 5 7 8 10,11])
% xticklabels({'reg','rand','reg','rand','reg','rand','rand','rand'})
legend([h{1}, h{2}, h{3}, h{4}],{'Set1','Set2', 'Set3','RAND-Control'})
ylabel('Peak Cross-correlation')
title([location{fiberid}, ' Cycle Size ', cyclesize])
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,1200,600])
% x = [0.7 2.3 3.7 5.3 6.7 8.3 9.7,11.3];
% hold on
% for i = 1:size(reg, 2)
%     %     scatter(x(i, i+1), [mean(reg(:,i)), mean(rand_d(:,i))], 'd','MarkerFaceColor', color(i,:))
%     errorbar(x(2*i-1:i*2),[mean(reg(:,i)), mean(rand_d(:,i))],([std(reg(:,i)), std(rand_d(:,i))])/sqrt(size(reg,1)), 'd','color',[0,0,0], 'MarkerFaceColor', [0,0,0], 'LineWidth', 1, 'CapSize',12)
% end
%
%
figure(2);
switch cyc
    case 4
        idx = 1;
    case 12
        idx = 2;
end
subplot(1,2,idx)
color =cbrewer('div', 'RdYlBu',4);
rectangle('Position',[0.5,-0.1, 5, 0.5], 'FaceColor',[color(3,:),0.5], 'EdgeColor',[color(3,:),0.5])
hold on
rectangle('Position',[5.5,-0.1, 5, 0.5], 'FaceColor',[0.5,0.5,0.5,0.5], 'EdgeColor',[0.5,0.5,0.5,0.5])
rectangle('Position',[11.5,-0.1, 5, 0.5], 'FaceColor',[color(3,:),0.5], 'EdgeColor',[color(3,:),0.5])
rectangle('Position',[16.5,-0.1, 5, 0.5], 'FaceColor',[0.5,0.5,0.5,0.5], 'EdgeColor',[0.5,0.5,0.5,0.5])

for i = 1:size(reg_avg,1)
    start = 0;
    for j = 1 :size(reg_avg,2)
        if j ==1
            h{j} = plot((1:10)+start, [squeeze(reg_avg(i,j,:)); squeeze(rand_avg(i,j,:))],'-o', 'color' ,color(j,:), 'MarkerFaceColor',color(j,:), 'LineWidth', 1)
        else
            h{j} = plot((1:10)+start, [squeeze(reg_avg(i,j,:)); squeeze(rand_avg(i,j,:))],'-o', 'color' ,color(4,:), 'MarkerFaceColor',color(4,:), 'LineWidth', 1)
            
        end
        hold on
        start = start +11;
    end
end
xlim([0,23])
ylim([-0.1,0.5])
xticks([1:10,12:21])
xticklabels({'1','2','3','4','5','1','2','3','4','5','1','2','3','4','5','1','2','3','4','5'})
xlabel('Frames (5 cycles)')
ylabel('Peak Cross-correlation')
text(2,0.42, 'REG', 'FontSize', 12)
text(7,0.42, 'RAND', 'FontSize', 12)
text(13,0.42, 'RAND', 'FontSize', 12)
text(18,0.42, 'RAND', 'FontSize', 12)
title([location{fiberid} 'Cycle size ', num2str(cyc)])
% legend([h{1}, h{2}],{'REG-RAND','RAND-Control'})
% ylabel('Peak Cross-correlation')
% title([location{fiberid}, ' Cycle Size ', cyclesize])
% box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,1200,600])


figure(3);
subplot(1,2,idx)
color =cbrewer('div', 'RdYlBu',4);
rectangle('Position',[0.5,-0.1, 5, 0.5], 'FaceColor',[color(3,:),0.5], 'EdgeColor',[color(3,:),0.5])
hold on
rectangle('Position',[5.5,-0.1, 5, 0.5], 'FaceColor',[0.5,0.5,0.5,0.5], 'EdgeColor',[0.5,0.5,0.5,0.5])
rectangle('Position',[11.5,-0.1, 5, 0.5], 'FaceColor',[color(3,:),0.5], 'EdgeColor',[color(3,:),0.5])
rectangle('Position',[16.5,-0.1, 5, 0.5], 'FaceColor',[0.5,0.5,0.5,0.5], 'EdgeColor',[0.5,0.5,0.5,0.5])
x = [1:10,12:21];
% hold on
for i = 1:size(reg_avg, 2)
    %     scatter(x(i, i+1), [mean(reg(:,i)), mean(rand_d(:,i))], 'd','MarkerFaceColor', color(i,:))
    mean_value = [squeeze(mean(reg_avg(:,i,:),1)); squeeze(mean(rand_avg(:,i,:),1))];
    sem_value = ([squeeze(std(reg_avg(:,i,:),0,1)); squeeze(std(rand_avg(:,i,:),0,1))])/sqrt(size(reg_avg,1))
    if i ==1
        errorbar(x(10*(i-1)+1:i*10),mean_value,sem_value, '-d','color',color(i,:), 'MarkerFaceColor', color(i,:), 'LineWidth', 1, 'CapSize',12)
    else
        errorbar(x(10*(i-1)+1:i*10),mean_value,sem_value, '-d','color',color(4,:), 'MarkerFaceColor', color(4,:), 'LineWidth', 1, 'CapSize',12)
        
    end
    hold on
end
xlim([0,23])
ylim([-0.1,0.5])
xticks([1:10,12:21])
xticklabels({'1','2','3','4','5','1','2','3','4','5','1','2','3','4','5','1','2','3','4','5'})
xlabel('Frames (5 cycles)')
ylabel('Peak Cross-correlation')
text(2,0.42, 'REG', 'FontSize', 12)
text(7,0.42, 'RAND', 'FontSize', 12)
text(13,0.42, 'RAND', 'FontSize', 12)
text(18,0.42, 'RAND', 'FontSize', 12)
title([location{fiberid} 'Cycle Size ', num2str(cyc)])
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,1200,600])
%% plot raw cross-correlation
% temp = data(3);
% type = {'REG', 'RAND'}
% for sets = 1:4
%     cyc = [4, 6, 10, 12];
%     color =cbrewer('div', 'RdYlBu',length(cyc));
%     figure
%     for i = 1:length(cyc)
%         for j = 1:2
%             subplot(length(cyc),2, (i-1) *2+ j)
%             trial = intersect(find([temp.dFibPhot.trial.sets] == sets), find([temp.dFibPhot.trial.cycTab] == cyc(i)))
%             dF = temp.dFibPhot.trial(trial).(type{j}).fiber2;
%             stimulus = temp.dFibPhot.trial(trial).(type{j}).stimuli;
%             Fs = temp.dFibPhot.fs;
%             signal = writeaudio_signal(stimulus, Fs);
%             [acor,lag] = xcorr(dF(1:length(signal)),signal,'coeff', 800);
%             plot(lag/Fs, acor, 'color', color(i,:))
%             ylim([-0.3, 0.3])
%             ylabel([' Cyc ', num2str(cyc(i))])
%             xlabel('Lags (s)')
%             if i == 1
%                 title([type{j}])
%             end
%         end
%     end
%     suptitle([temp.animal, ' ', 'Set ', num2str(sets), ' Caudal ', 'XCorr(dF,sound)'])
%     set(gcf, 'Position',[200, 200, 600, 800])
% end
% [pxx, f] = periodogram(signal,[], length(signal),Fs);
% [normalizedACF, lags] =  xcorr(signal+ 0.1* normrnd(0, 1, length(signal),1),signal,'coeff', 5000);
% % figure;
% % plot(f, pxx)
% % xlim([0.4, 6])
% figure;
% plot(lags/Fs, normalizedACF)

%% plot cycle by cycle
temp = data(3);
type = 'REG'

sets = 1
cyc = 12
trial = intersect(find([temp.dFibPhot.trial.sets] == sets), find([temp.dFibPhot.trial.cycTab] == cyc))
dF = temp.dFibPhot.trial(trial).(type).fiber2;
stimulus = temp.dFibPhot.trial(trial).(type).stimuli;
% plot_patternSeq(diff(stimulus), cyc)
time = temp.dFibPhot.trial(trial).(type).time;
signal_psth = [];
for i = 1: 25*cyc %length(stimulus)
    idx = find(time > stimulus(i), 1,'first');
    signal_psth(i,:) = dF(idx:idx + 100);
end

harmonic = 1;
color =cbrewer('div', 'RdYlBu',cyc*harmonic);
figure;
for i = 1:25*cyc %length(stimulus)
    rem = mod(i,harmonic*cyc);
    if rem == 0
        subplot(cyc/2* harmonic, 2, cyc*harmonic)
        plot((signal_psth(i,:)-min(signal_psth(i,:)))./max((signal_psth(i,:)-min(signal_psth(i,:)))), 'color',color(cyc,:))
        title(['Within Cycle ', num2str(cyc)])
        %         plot(signal_psth(i,:), 'color',color(4,:))
        xlim([0, 100])
        xlabel('Time (ms)')
        ylabel('Normalized dF')
    else
        subplot(cyc/2 * harmonic, 2, rem)
        plot((signal_psth(i,:)-min(signal_psth(i,:)))./max((signal_psth(i,:)-min(signal_psth(i,:)))), 'color',color(rem,:))
        title(['Within Cycle ', num2str(rem)])
        %         plot(signal_psth(i,:), 'color',color(rem,:))
        xlim([0, 100])
        xlabel('Time (ms)')
        ylabel('Normalized dF')
        
    end
    hold on
end
xlim([0, 100])
suptitle([temp.animal, ' ', 'Set ', num2str(sets), ' Caudal ', type])
%%
avg_signal_psth = mean(signal_psth,1);
hold on
plot(avg_signal_psth, 'color',[0, 0, 0])

%%
Fs = temp.dFibPhot.fs;
signal = writeaudio_signal(stimulus, Fs);

[acor,lag] = xcorr(dF(1:end-1),signal,'coeff', 5000);
figure
plot(lag/Fs, acor)
% [pxx, f] = periodogram(signal,[], length(signal),Fs);
[normalizedACF, lags] =  xcorr(signal,signal,'coeff', 5000);
% figure;
% plot(f, pxx)
% xlim([0.4, 6])
figure;
plot(lags/Fs, normalizedACF)

%% new analysis based on jittered stimulus
% Add that folder plus all subfolders to the path.
folder = 'F:\KeChen\MATLAB\Tosca\MATLAB'
addpath(genpath(folder));
runFile = 'Z:\File Transfer\Blaise\RANDREG_forKe\Saved\DFI045\Session 6\DFI045-Session6-Run2.txt';
struct_tosca=sortFilesRandRegKe(runFile)