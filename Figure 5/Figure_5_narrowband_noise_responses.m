%% Analyze the noise exposure data 
clear all 
close all

%CHANGE THIS PATH 
 load('\\apollo\Polley_Lab\Mouse Videography Audiogram\Data\Data Summary\sham_data_unfilt_z.mat')
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 5\';% 'C:\Users\claytok\Dropbox\File_Transfer\Dan - Kameron\Staff scientist\Conference Presentations\ARO23\Facial_videography\Figure5\Facial_videography\'
%Average all data across trials
sham_data = squeeze(nanmean(nbn1_facial_z_unfilt,5));
sham_data_tbyt = nbn1_facial_z_unfilt;

%Generate RLFs 
time_win = 151:200; 
for m = 1:size(sham_data,1)
    for d = 1:5
        for i = 1:size(sham_data,4) 
            
           clear temp_ind
    % do for 8 kHz 
           [dum temp_ind] = max(sham_data(m,d,1,i,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        sham_8kHz_rlf(m,d,i) = nanmean(sham_data(m,d,1,i,temp_ind-2:temp_ind+2))%-...
           %mean(sham_data(m,d,1,i,1:150)); 
        
        %do for 32 kHz 
        clear temp_ind
        [dum temp_ind] = max(sham_data(m,d,2,i,time_win));
        temp_ind = temp_ind +150; % Correct for ROI
       sham_32kHz_rlf(m,d,i) = nanmean(sham_data(m,d,2,i,temp_ind-2:temp_ind+2))%-...
        %mean(sham_data(m,d,2,i,1:150))
        end 
    end
end

% Get baseline functions 
baseline_rlf_8_kHz_sham = squeeze(nanmean(sham_8kHz_rlf(:,1:2,:),2)); 
baseline_rlf_32_kHz_sham = squeeze(nanmean(sham_32kHz_rlf(:,1:2,:),2)); 


%% Now load the exposed data and do the same thing 
load('\\apollo\Polley_Lab\Mouse Videography Audiogram\Data\Data Summary\exposed_data_unfilt_z.mat')
exposed_data = squeeze(nanmean(nbn1_facial_z_unfilt_exposed,5));
exposed_data_tbyt = nbn1_facial_z_unfilt_exposed;


%Sham scatter plots 
%turn all data into rlfs at 8 kHz
time_win = 151:300; 
for m = 1:size(exposed_data,1)
    for d = 1:5
        for i = 1:size(exposed_data,4) 
            
           clear temp_ind
        % do for 8 kHz 
           [dum temp_ind] = max(exposed_data(m,d,1,i,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        exposed_8kHz_rlf(m,d,i) = mean(exposed_data(m,d,1,i,temp_ind-2:temp_ind+2))%-...
           %mean(exposed_data(m,d,1,i,1:150)); 
        
        %do for 32 kHz 
        clear temp_ind
        [dum temp_ind] = max(exposed_data(m,d,2,i,time_win));
        temp_ind = temp_ind +150; % Correct for ROI
        exposed_32kHz_rlf(m,d,i) = mean(exposed_data(m,d,2,i,temp_ind-2:temp_ind+2))%-...
       % mean(exposed_data(m,d,2,i,1:150))
     
        end 
    end
end

baseline_exp_rlf_8_kHz = squeeze(nanmean(exposed_8kHz_rlf(:,1:2,:),2)); 
baseline_exp_rlf_32_kHz = squeeze(nanmean(exposed_32kHz_rlf(:,1:2,:),2)); 

%% figure A) Avg. across baseline days and frequencies to show growth function

%Get heat map for both group 
 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)

sham_time_series_32 = squeeze(nanmean(nanmean(sham_data(:,2,:,1:6,:),3),2));
exposed_time_series_32 = squeeze(nanmean(nanmean(exposed_data(:,2,:,1:6,:),3),2))

octave_band_resp = cat(1,sham_time_series_32,exposed_time_series_32) 

%Get the tones and noise time series data 
load('\\Apollo\Polley_lab\Mouse Videography Audiogram\Data\Data Summary\z_spectral_stim_unfilt.mat')
% Get mean responses
wn_mean_resp = squeeze(nanmean(z_spectral_stim.data.wn,3)); 
tones_mean_resp = squeeze(nanmean(z_spectral_stim.data.tones,4)); 
tones_32_resp = squeeze(tones_mean_resp(:,2,:,:)); 

%70 dB 
figure;
h1 = shadedErrorBar(1:size(octave_band_resp,3),squeeze(mean(wn_mean_resp(:,3,:))),...
squeeze(std(wn_mean_resp(:,3,:)))./sqrt(8))
hold on
h2 = shadedErrorBar(1:size(octave_band_resp,3),squeeze(mean(octave_band_resp(:,5,:))),...
squeeze(std(octave_band_resp(:,5,:)))./sqrt(20),'r')
h3 = shadedErrorBar(1:size(octave_band_resp,3),squeeze(mean(tones_32_resp(:,3,:))),...
squeeze(std(tones_32_resp(:,3,:)))./sqrt(8),'b')
xlim([75 301])
xticks([76:75:451])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('TME (z-score)') 
box off 
legend([h1.mainLine h2.mainLine h3.mainLine],{'WN','1 oct noise','Tone'}) 
set(gca,'xcolor','k','ycolor','k')


%Get avg response for exposed group in baseline
baseline_exp_group = (baseline_exp_rlf_8_kHz + baseline_exp_rlf_32_kHz)./2;
%Get avg response for sham group in baseline 
baseline_sham_group = (baseline_rlf_8_kHz_sham + baseline_rlf_32_kHz_sham)./2; 

avg_response = [baseline_exp_group; baseline_sham_group]; 

figure;
shadedErrorBar(30:10:80,mean(avg_response(:,1:6)), std(avg_response(:,1:6))./sqrt(size(avg_response,2)))
hold on 
title('Intensity growth function for 1 oct NBN evoked facial movements')
xticks([30:10:80])
xlim([25 85]) 
xlabel('dB SPL')
ylabel('Facial motion energy (z-score)') 
box off 
set(gca,'xcolor','k','ycolor','k')


%Repeated measures ANOVA with main effect of sound level 
%Now do stats on above
allValues = [avg_response] 

%allValues([1 3 11 13],:) = []; 
%Now we want to format this into a table with the correct labels
for i = 1:size(allValues,2)
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Noise_level'};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
freqLabels = [1:size(allValues,2)]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V8~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Noise_level');
rAnovaResults

 
%% Figure 5C) for sham mice, plot each days
figure;
plot(1:5,squeeze(sham_8kHz_rlf(:,:,5))','Color',[0.5 0.5 1])
hold on
plot(1:5,squeeze(sham_32kHz_rlf(:,:,5))','Color',[1 0.5 0.5])
h1 = errorbar(1:5,mean(mean(sham_8kHz_rlf(:,:,5),3)),std(mean(sham_8kHz_rlf(:,:,5),3))./sqrt(8),...
    'LineWidth',2,'Color','b');
h2 = errorbar(1:5,mean(mean(sham_32kHz_rlf(:,:,5),3)),std(mean(sham_32kHz_rlf(:,:,5),3))./sqrt(8),...
    'LineWidth',2,'Color','r');
ylim([0 3.5]) 
box off 
ylabel('Facial mvmt response at 70 dB SPL (z-score)') 
set(gca,'xcolor','k','ycolor','k')
xlim([0.5 5.5]) 
xticks([1:5]) 
xticklabels({'D1','D2','D3','D8','D14'})
xlabel('Measurement day') 
legend([h1 h2], {'8 kHz','32 kHz'}) 

%Now perform a repeated measures ANOVA with day as factor 
clear allValues
clear varNames
allValues = [squeeze(sham_8kHz_rlf(:,:,5)) squeeze(sham_32kHz_rlf(:,:,5))]

%allValues([1 3 11 13],:) = []; 
%Now we want to format this into a table with the correct labels
for i = 1:size(allValues,2)
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable
dataTable.Properties.VariableNames =["lo_d1","lo_d2","lo_d3","lo_d4","lo_d5",...
    "hi_d1","hi_d2","hi_d3","hi_d4","hi_d5"]; 

within = table([1 1 1 1 1 2 2 2 2 2]',[1 2 3 4 5 1 2 3 4 5]','VariableNames',{'Freq' 'Day'}) % within model
within.Freq = categorical(within.Freq);
within.Day = categorical(within.Day);

temp_rm = fitrm(dataTable,'lo_d1-hi_d5~1','WithinDesign',within,'WithinModel','Freq*Day');

[ranovatbl,A,C,D] = ranova(temp_rm,'WithinModel','Freq*Day')
