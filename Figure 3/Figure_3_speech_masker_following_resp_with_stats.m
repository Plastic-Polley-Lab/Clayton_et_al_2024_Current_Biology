clear all
close all
%CHANGE THESE PATHS 
files = {'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_151_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_152_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_153_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_154_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_155_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_156_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_157_speech_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\speech_in_noise\ChAT_150_speech_in_noise.mat'}

files = unique(files); %ensure no duplicates

%% Load in data

sprate_mat = [ ];
for m = [1:8]% length(files)
    clear dat
    load(files{m})
    sprate_mat(m,:,:,:,:) = dat.sorted_face_stim;
end


%% Average across trials 
mean_srm = mean(sprate_mat,4)-mean(sprate_mat(:,:,:,1:800),4);

figure;
cnt = 1; 
for i = 1:2
    for j = 1:6
        subplot(2,6,cnt)
       imagesc(squeeze(mean_srm(:,i,j,:)))
       caxis([0 2])
        cnt = cnt + 1;
        xlim([800 1750])
    end
end

masker_levs = [-inf 10 20 30 40 50]; 

figure;
cnt = 1; 
for i = 1
    for j = [1:6]; 
subplot(6,1,cnt)
h1 = shadedErrorBar(1:1801,mean(squeeze(mean_srm(:,1,j,:))),...
   std(squeeze(mean_srm(:,1,j,:)))./sqrt(8),'lineProps','r')
hold on
h2 = shadedErrorBar(1:1801,mean(squeeze(mean_srm(:,2,j,:))),...
   std(squeeze(mean_srm(:,2,j,:)))./sqrt(8),'lineProps','b')
ylim([-.5 3]) 
xlim([651 1750]) 
xticks([651:150:1800]) 
xticklabels([-1:6]) 
hold on
if j == 6
    cnt = 0;
    ylabel('Facial mvmt amplitude (z-score)') 
end
cnt = cnt + 1;  
    title([num2str(masker_levs(j)) 'dB Masker Level'])
    end
    
    
end
xlabel('Time re: sentence onset') 
legend([h1.mainLine h2.mainLine],{'Gee','Ha'})
% 
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\f3_elements\';
%print(gcf,[savedir 'Speech_in_noise_examples.pdf'],'-dpdf')



%% Frequency analysis 
for m = 1:8
    cnt = 1;
    for i = 1:2
        for j = 1:6
        clear P2 P1 f
        Fs = 150;
        T = 1/Fs;
        L = 900 ;
        t = (0:L-1)*T;
        
        Y = fft(squeeze(mean_srm(m,i,j,800:1700)));
        
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        fourier_sync{i,j}.psd(m,:) = P1;
        f = Fs*(0:(L/2))/L;
        
        fourier_sync{i,j}.f = f; 
        clear Y
        end
    end
end


%% Plot ffts 
figure;
cnt = 1;
%masker_level
masker_levs = [-inf 10 20 30 40 50]; 
for i  = 1:2
    for j = 1:6
    subplot(2,6,cnt)
    plot(f,squeeze(fourier_sync{i,j}.psd'),'Color',[.5 .5 .5])
    hold on
    plot(f,mean(squeeze(fourier_sync{i,j}.psd)),'k','LineWidth',2)
    xlim([0 10])
    xlabel('Frequency (Hz)')
    ylabel('Power')
    cnt = cnt + 1;
    ylim([0 0.6]) 
    xline(1,'r--')
    
    if i == 1
        title([num2str(masker_levs(j)) ' dB masker'])
    end
end
end



%% Compute db SNR 

cnt = 1 %Account for multiple 
shoulder = [1]
for i = 1:2
    for j = 1:6
    
    sig_amp(i,j,:) =(fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1)))-...
        mean([fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1)-shoulder),...
        fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1)+shoulder)],2);
        
    dB_SNR(i,j,:) =10*log10(fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1))./...
        mean([fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1)-shoulder),...
        fourier_sync{i,j}.psd(:,find(round(fourier_sync{i,j}.f,2)==1)+shoulder)],2)); 
  
    if cnt == 7 
        cnt = 0; 
    end
    
    cnt = cnt + 1;

    end
end



%% Plot entrainment values for both speech tokens
figure; 
subplot(1,2,1)
errorbar(1:6,mean(dB_SNR(1,:,:),3),std(dB_SNR(1,:,:),[],3)./sqrt(8),'k','LineWidth',2)
hold on 
plot(squeeze(dB_SNR(1,:,:)),'Color',[0.5 0.5 0.5])
xticks([1:6])
xticklabels(masker_levs) 
xlim([0.75 6.25]) 
title('Gee') 
box off 
ylabel('Speech following response, dB SNR') 
ylim([-5 15]) 

subplot(1,2,2) 
errorbar(1:6,mean(dB_SNR(2,:,:),3),std(dB_SNR(2,:,:),[],3)./sqrt(8),'k','LineWidth',2)
hold on 
plot(squeeze(dB_SNR(2,:,:)),'Color',[0.5 0.5 0.5])
xticks([1:6])
xticklabels(masker_levs) 
xlim([0.75 6.25]) 
title('Ha') 
box off 
ylabel('Speech following response, dB SNR') 
ylim([-5 15]) 
%print(gcf,[savedir 'Speech_following_response.pdf'],'-dpdf')


%% Do statistics 

%Now do stats on above
allValues = [squeeze(dB_SNR(1,:,:))'; squeeze(dB_SNR(2,:,:))'];
%Now we want to format this into a table with the correct labels
for i = 1:6
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:16
    if ismember(i,[1:8])
        expGroup{i,1} = 'Gee';
    else
        expGroup{i,1} = 'Ha';
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Noise_level'};
freqLabels = [1:6]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V6~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Noise_level');
rAnovaResults
%Also do post-hoc tests for each frequency. Do as a two sample t test
%comparing noise to sham. One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:6 %For 6 frequencies
    nums1 = allValues(1:7,i); nums2 = allValues(8:end,i);
    [h,p] = ttest2(nums1,nums2);
    disp(freqLabels{i}); disp(p); disp('');
end
