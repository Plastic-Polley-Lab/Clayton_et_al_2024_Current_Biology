%% Summarize basic data for the FM at different rates 

clear all 
close all 

%CHANGE THIS PATH 
files = {'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_151_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_152_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_153_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_154_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_155_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_156_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_157_FM_rate.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\FM_rate_mod\FM_rate_modChAT_150_FM_rate.mat'}


files = unique(files) %ensure no duplicates

%% Load in data 

sprate_mat = [ ] 
for m = 1:length(files)
    clear dat
   load(files{m})
   sprate_mat(m,:,:,:,:) = dat.sorted_face_stim; 
    
end
    
%% just plot out means to get a sense of the data 
mean_srm = squeeze(mean(sprate_mat,4)); %Average over trials

rep_rate = 4:-0.5:3; 
for ts = 1:3 
    time_string(ts,:) = 151:150/rep_rate(ts):(151+150/rep_rate(ts)*6)
end

figure;
cnt = 1;
for i = 1:2
    for j = 3:-1:1
        subplot(2,3,cnt)
        imagesc(squeeze(mean_srm(:,i,j,:)))
        caxis([0 2])
        cnt = cnt + 1;
        xlim([0 1500])
        xlim([76 600])
        xticks([151:150:1076])
        xticklabels([0:6])
        xlabel('Time re: stimulus onset')
        colormap('hot')
        xline(151,'w','LineWidth',2)
        hold on
        for t = 1:6
           xline(time_string(j,t),'w','LineWidth',2) 
        end
        if i == 1;
            ylabel('Down-sweeps')
        end
        
        if i == 4;
            ylabel('Up-sweeps')
        end
        
        if i<3
            title([num2str(rep_rate(j)) ' Hz' ])
        end
        
    end
end



figure;
cnt = 1;

for j = 1:7;
    subplot(7,1,j)
    shadedErrorBar(1:4500,squeeze(mean(mean(mean_srm(:,:,8-j,:),2))),...
        squeeze(std(mean(mean_srm(:,:,8-j,:),2)))./sqrt(8),'LineProps','r')
    ylim([-.5 2])
    box off
    xlim([76 1051])
    xticks([151:150:1076])
    xticklabels(string([0:6]))
    hold on
    if j == 7
        cnt = 0;
    end
end
xlabel('Time re sound onset')
ylabel('Facial mvmt amp.')
set(gca,'xcolor','k','ycolor','k')
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 2\';
set(gcf,'Position',[162   149   560   945])
print(gcf,[savedir 'FM_mean_resp.pdf'],'-dpdf')


%% Get the synchronization rates for each individual stimuli and mouse
windows = floor(6./[4:-0.5:1]*150); %Set the windows for FFT analysis
%windows = 6*ones(1,7)*150 

for m = 1:8
    cnt = 1;
    for i = 1:2
        for j = 1:7 
        clear P2 P1 f L t f 
        Fs = 150;
        T = 1/Fs;
        L = windows(j);
        t = (0:L-1)*T;
        
        Y = fft(squeeze(mean_srm(m,i,j,151:151+windows(j)-1)));
        
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        fourier_sync{i,j}.psd(m,:) = P1;
        fourier_sync{i,j}.f = Fs*(0:(L/2))/L;
        
        clear Y

    end
    end
end



figure;
cnt = 1;
freqs = 1:0.5:4;
for i  = 1:7
    subplot(7,1,cnt)
    clear temp_dat
    temp_dat = (fourier_sync{1,8-i}.psd+fourier_sync{1,8-i}.psd)./2
    plot(fourier_sync{1,8-i}.f,temp_dat,'Color',[.5 .5 .5])
    hold on
    plot(fourier_sync{1,8-i}.f,mean(temp_dat),'LineWidth',2)
    xlim([0 10])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    ylim([0 0.5])
    cnt = cnt +1;
    
    if i <= 7 
        title([num2str(freqs(i)) ' Hz'])
    end
end
set(gcf,'Position',[162   149   560   945])
%print(gcf,[savedir 'FM_rate_fft.pdf'],'-dpdf')



cnt = 1 %Account for multiple 
shoulder = [1]
for i = 1:2
    for j = 1:7 
    
    sig_amp(i,j,:) =(fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j))))-...
        mean([fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j))-shoulder),...
        fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j))+shoulder)],2);
        
    dB_SNR(i,j,:) =10*log10(fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j)))./...
        mean([fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j))-shoulder),...
        fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==freqs(8-j))+shoulder)],2)); 
  
    if cnt == 7 
        cnt = 0; 
    end
    
    cnt = cnt + 1;

    end
end



dB_SNR_all = squeeze(mean(dB_SNR,1));
sig_amp_all = squeeze(mean(sig_amp,1)); 



figure
errorbar(freqs,flipud(mean(sig_amp_all,2)),flipud(std(sig_amp_all,[],2))./sqrt(8),'k','LineWidth',2)
box off 
hold on
plot(flipud(freqs),flipud(sig_amp_all),'Color',[0.5 0.5 0.5])
xlabel('Stimulus Frequency (hz)')
xlim([0.75 4.25]) 
ylabel('Phase-locking @ stimulus frequency') 
title('FM rate following response') 



figure
errorbar(freqs,flipud(mean(dB_SNR_all,2)),flipud(std(dB_SNR_all,[],2))./sqrt(8),'k','LineWidth',2)
box off 
hold on
plot(flipud(freqs),flipud(dB_SNR_all),'Color',[0.5 0.5 0.5])
xlabel('Stimulus Frequency (hz)')
xlim([0.75 4.25]) 
ylabel('Phase-locking @ stimulus frequency (dB SNR)') 
title('FM rate following response') 

%print(gcf,[savedir 'FM_rate_fft_dB_SNR.pdf'],'-dpdf')

%% Statistics 



allValues = [flipud(dB_SNR_all)'];
%Now we want to format this into a table with the correct labels
for i = 1:7
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
factorNames = {'Temporal_modulation_rate'};
mvmtLabels = [1:7]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V7~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Temporal_modulation_rate');
rAnovaResults










