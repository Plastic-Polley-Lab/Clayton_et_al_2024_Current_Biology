%% Analyze the noise exposure data 
clear all 
close all

%CHANGE THIS PATH 
 load('\\apollo\Polley_Lab\Mouse Videography Audiogram\Data\Data Summary\sham_data_unfilt_z.mat')
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 5\';% 'C:\Users\claytok\Dropbox\File_Transfer\Dan - Kameron\Staff scientist\Conference Presentations\ARO23\Facial_videography\Figure5\Facial_videography\'
%Average all data across trials
sham_data = squeeze(nanmean(nbn1_facial_z_unfilt,5));
sham_data_tbyt = nbn1_facial_z_unfilt;

%% Get response thresholds and latencies 
for f = 1:2 
    for d = 1:5
        h = getResponseStats(squeeze(sham_data_tbyt(:,d,f,:,:,:)))
        threshold_sham(d,f,:) = h.threshold; 
        latency_sham(d,f,:,:) = h.latency; 
        disp('here')
    end
end
%%
%Sham scatter plots 
%turn all data into rlfs at 8 kHz
time_win = 151:300; 
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


%% Plot fold change 

%Calculate fold change at each time point for each mouse 

for d = 1:5 
    clear temp_data
    temp_data = squeeze(sham_8kHz_rlf(:,d,:));
    
    fold_change_sham(d,:,:) = temp_data./baseline_rlf_8_kHz_sham;
   
end 

for d = 1:5 
    clear temp_data
    temp_data = squeeze(sham_32kHz_rlf(:,d,:));
    
    fold_change_32sham(d,:,:) = temp_data./baseline_rlf_32_kHz_sham;
   
end 


%% Now load the exposed data and do the same thing 
load('\\apollo\Polley_Lab\Mouse Videography Audiogram\Data\Data Summary\exposed_data_unfilt_z.mat')
exposed_data = squeeze(nanmean(nbn1_facial_z_unfilt_exposed,5));

exposed_data_tbyt = nbn1_facial_z_unfilt_exposed;

%Now eliminate mice that don't have 2 wk data 
exposed_data  = exposed_data([1:3 6:12],:,:,:,:); 
exposed_data_tbyt  = exposed_data_tbyt([1:3 6:12],:,:,:,:,:); 

for f = 1:2 
    for d = 1:5
        h = getResponseStats(squeeze(exposed_data_tbyt(:,d,f,:,:,:)))
        threshold_exposed(d,f,:) = h.threshold; 
        latency_exposed(d,f,:,:) = h.latency; 
        disp('here')
    end
end


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
 cmap(cmap<0) = 0 ; 

sham_time_series = squeeze(nanmean(nanmean(sham_data(:,1:2,:,1:6,:),3),2));
exposed_time_series = squeeze(nanmean(nanmean(exposed_data(:,1:2,:,1:6,:),3),2))

octave_band_resp = cat(1,sham_time_series,exposed_time_series) 
imagesc(flipud(squeeze(mean(octave_band_resp))))
xlim([75 301])
xticks([76:75:451])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticks([1:7]) 
yticklabels(string([80:-10:30]))
colormap(cmap)
caxis([-2 2 ]) 
box off 
colorbar
set(gca,'xcolor','k','ycolor','k')
print(gcf,[savedir 'Time_series_response_Oct_band_noise.pdf'],'-dpdf')



%Get avg response for exposed group 
baseline_exp_group = (baseline_exp_rlf_8_kHz + baseline_exp_rlf_32_kHz)./2;
%Get avg response for 
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
%legend([h(:).mainLine],{'WN','Tone','0.5 Oct NBN','1 Oct NBN'},'Location','Northwest')
box off 
set(gca,'xcolor','k','ycolor','k')
%print(gcf,[savedir 'Oct_band_noise_rlf.pdf'],'-dpdf')

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

 
%% Figure B) for sham mice, plot each level
%Sham_8kHz_rlf 
sham_rlf = (sham_8kHz_rlf + sham_32kHz_rlf)./2;

figure;
plot(1:5,squeeze(sham_rlf(:,:,5))','Color',[0.5 0.5 0.5])
hold on
errorbar(1:5,mean(mean(sham_rlf(:,:,5),3)),std(mean(sham_rlf(:,:,5),3))./sqrt(8),...
    'LineWidth',2,'Color','k')
ylim([0 3.5]) 
box off 
ylabel('Facial mvmt response at 70 dB SPL (z-score)') 
set(gca,'xcolor','k','ycolor','k')
xlim([0.5 5.5]) 
xticks([1:5]) 
xticklabels({'D1','D2','D3','D8','D14'})
xlabel('Measurement day') 
print(gcf,[savedir 'SEFM_over_days.pdf'],'-dpdf')

%Now perform a repeated measures ANOVA with day as factor 
clear allValues
clear varNames
allValues = [squeeze(sham_rlf(:,:,5))] 

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
factorNames = {'Day'};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
freqLabels = [1:size(allValues,2)]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V5~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Day');
rAnovaResults


%% Now plot mean response heat maps 
figure
days = {'Baseline Day 1','Baseline day 2','2 hr','5 days','2 wk'}; 
for d = 1:5
    subplot(1,5,d)
       imagesc(flipud(squeeze(nanmean(exposed_data(:,d,1,1:6,:)))))
       hold on 
       xlim([75 301]) 
       xticks([76 151 151+75]) 
       xticklabels({'-0.5','0','0.5','1.0'}) 
       colorbar
       caxis([-3 3]) 
       %ylim([-1 7]) 
       colormap(cmap)
       yticks([1:6]) 
       yticklabels(string([80:-10:30])) 
       title(days{d}) 
       if d == 1
           ylabel('Intensity (dB SPL)') 
           xlabel('Time re: sound onset') 
       end
end   
sgtitle('Exposed data over days')
set(gca,'xcolor','k','ycolor','k')
 set(gcf,'Position',[54         848        1548         250])
print(gcf,[savedir '8 kHz heat_map_exposed.pdf'],'-dpdf')


figure
days = {'Baseline Day 1','Baseline day 2','2 hr','5 days','2 wk'}; 
for d = 1:5
    subplot(1,5,d)
       imagesc(flipud(squeeze(nanmean(exposed_data(:,d,2,1:6,:)))))
       hold on 
       xlim([75 301]) 
       xticks([76 151 151+75]) 
       xticklabels({'-0.5','0','0.5','1.0'}) 
       colorbar
       caxis([-3 3]) 
       %ylim([-1 7]) 
       colormap(cmap)
       yticks([1:6]) 
       yticklabels(string([80:-10:30])) 
       title(days{d}) 
       if d == 1
           ylabel('Intensity (dB SPL)') 
           xlabel('Time re: sound onset') 
       end
end   
sgtitle('Exposed data over days')
 set(gcf,'Position',[54         848        1548         250])
 set(gca,'xcolor','k','ycolor','k')
print(gcf,[savedir '32 kHz heat_map_exposed.pdf'],'-dpdf')


%% Now plot fold change for each mouse 
%Calculate fold change at each time point for each mouse 

for d = 1:5 
    clear temp_data
    temp_data = squeeze(exposed_8kHz_rlf(:,d,:));
    fold_change(d,:,:) = temp_data./baseline_exp_rlf_8_kHz;
   
end 

for d = 1:5 
    clear temp_data
    temp_data = squeeze(exposed_32kHz_rlf(:,d,:));
    fold_change_32(d,:,:) = temp_data./baseline_exp_rlf_32_kHz;   
end 


%% Plot change in threshold 
threshold_exposed(threshold_exposed == 125) = 100; 
threshold_sham(threshold_sham == 125) = 100; 

% Now plot threshold shifts 
threshold_shift_sham = threshold_sham - nanmean(threshold_sham(1:2,:,:),1)

% Now plot threshold shifts 
threshold_shift_ex = threshold_exposed - nanmean(threshold_exposed(1:2,:,:),1)

%Now do stats on the 8 kHz responses 

%Now, do stats on the pre vs. 2 wk post thresholds 
clear allValues varNames dataTable
noise_thr = squeeze(threshold_shift_ex([1 2 3 5],1,:));
sham_thr = squeeze(threshold_shift_sham([1 2 3 5],1,:)); 

allValues = [noise_thr sham_thr]; 
%Now we want to format this into a table with the correct labels
for i = 1:4
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues', 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:18
    if ismember(i,[1:10])
        expGroup{i,1} = 'Noise';
    else
        expGroup{i,1} = 'Sham';
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Day'};
freqLabels = [1:4]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V4~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Day');
rAnovaResults


%% Repeat for 32 kHz 

%Now, do stats on the pre vs. 2 wk post thresholds 
clear allValues varNames dataTable
noise_thr = squeeze(threshold_shift_ex([1 2 3 5],2,:));
sham_thr = squeeze(threshold_shift_sham([1 2 3 5],2,:)); 

allValues = [noise_thr sham_thr]; 
%Now we want to format this into a table with the correct labels
for i = 1:4
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues', 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:18
    if ismember(i,[1:10])
        expGroup{i,1} = 'Noise';
    else
        expGroup{i,1} = 'Sham';
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Day'};
freqLabels = [1:4]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V4~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Day');
rAnovaResults


%% Now express gain changes relative to threshold for 8 kHz

%Does it clean up gain change result to show threshold
levels = 30:10:100; 
fold_change_8k_re_thresh = nan(5,10,8)
for m = 1:10
    for d = 1:5
        ind = find(threshold_exposed(d,1,m)==levels)
        fold_change_8k_re_thresh(d,m,1:length(ind:8)) = fold_change(d,m,ind:end)
        
    end
end


levels = 30:10:100; 
fold_change_8k_re_thresh_sham = nan(5,8 ,8)
for m = 1:8
    for d = 1:5
        ind = find(threshold_sham(d,1,m)==levels)
        fold_change_8k_re_thresh_sham(d,m,1:length(ind:8)) = fold_change_sham(d,m,ind:end)
        %latency_8k_re_thresh(d,m,1:length(ind:8)) = latency_exposed(d,1,m,ind:end); 
        
    end
end


levels = 30:10:100; 
fold_change_32k_re_thresh = nan(5,10,8)
for m = 1:10
    for d = 1:5
        ind = find(threshold_exposed(d,2,m)==levels)
        fold_change_32k_re_thresh(d,m,1:length(ind:8)) = fold_change_32(d,m,ind:end)
        
    end
end


levels = 30:10:100; 
fold_change_32k_re_thresh_sham = nan(5,8 ,8)
for m = 1:8
    for d = 1:5
        ind = find(threshold_sham(d,1,m)==levels)
        fold_change_32k_re_thresh_sham(d,m,1:length(ind:8)) = fold_change_32sham(d,m,ind:end)
        %latency_8k_re_thresh(d,m,1:length(ind:8)) = latency_exposed(d,1,m,ind:end); 
        
    end
end


figure; 
hold on 
tps = {'0 dB SL','10 dB SL','20 dB SL'};
for l = 1:3
subplot(1,3,l)

hold on 
for m = 1:10
    ind = find(~isnan(fold_change_8k_re_thresh(:,m,l)));
    plot(ind,squeeze(fold_change_8k_re_thresh(ind,m,l)),'Color',[0.5 0.5 0.5])
    hold on
end 
h2 = errorbar(1:5,squeeze(nanmean(fold_change_8k_re_thresh(:,:,l),2)),...
    squeeze(nanstd(fold_change_8k_re_thresh(:,:,1),[],2)./sqrt(8)),'k','LineWidth',2)
ylim([0 6])
xticks([1:5]) 
xlim([0.5 5.5]) 
xticklabels({'-2','-1','2 hr','5 days','2 wk'}); 
yline(1,'--')
ylabel('Fold change re: threshold')
xlabel('Time point') 
%legend([h2],{'8 kHz'}) 
box off
 set(gca,'xcolor','k','ycolor','k')
title(tps{l}) 
end
sgtitle('8 kHz responses re: threshold, exposed') 
set(gcf,'Position',[129         680        1026         324])
print(gcf,[savedir 'exposed_fold_change_re_thresh.pdf'],'-dpdf')
 %set(gcf,'Position',[54         848        1548         250])

figure; 
hold on 
tps = {'0 dB SL','10 dB SL','20 dB SL'};
for l = 1:3
subplot(1,3,l)

hold on 
for m = 1:10
    ind = find(~isnan(fold_change_32k_re_thresh(:,m,l)));
    plot(ind,squeeze(fold_change_32k_re_thresh(ind,m,l)),'Color',[0.5 0.5 0.5])
    hold on
end 
h2 = errorbar(1:5,squeeze(nanmean(fold_change_32k_re_thresh(:,:,l),2)),...
    squeeze(nanstd(fold_change_32k_re_thresh(:,:,1),[],2)./sqrt(8)),'k','LineWidth',2)
ylim([0 6])
xticks([1:5]) 
xticklabels({'-2','-1','2 hr','5 days','2 wk'}); 
yline(1,'--')
xlim([0.5 5.5]) 
ylabel('Fold change re: threshold')
xlabel('Time point') 
%legend([h2],{'8 kHz'}) 
box off
title(tps{l}) 
 set(gca,'xcolor','k','ycolor','k')
end
sgtitle('32 kHz responses re: threshold, exposed') 
set(gcf,'Position',[129         680        1026         324])
print(gcf,[savedir '32_kHz_exposed_fold_change_re_thresh.pdf'],'-dpdf')
 
 %% Now do 1 sample t-tests with population mean of 1
 for d = 1:5 
     for l = 1:3
  [t p_32(d,l)] = ttest(fold_change_32k_re_thresh(d,:,l),1); 
     end
 end
 
  for d = 1:5 
     for l = 1:3
  [t p_8(d,l)] = ttest(fold_change_8k_re_thresh(d,:,l),1); 
     end
  end

%% Now plot scatters 
intens_range = [1:2]; 
%8 kHz sham and 
figure
for m = 1:10 
plot([0.875+rand(1,1)*0.25 1.875+rand(1,1)*0.25],[nanmean(fold_change_8k_re_thresh(5,m,intens_range),3)...
    nanmean(fold_change_32k_re_thresh(5,m,intens_range),3)],'-r.','MarkerSize',25)
hold on
end
plot([0.875 1.125],nanmean(nanmean(fold_change_8k_re_thresh(5,:,intens_range),3)).*ones(1,2),'k',...
    'LineWidth',2)
plot([1.875 2.125],nanmean(nanmean(fold_change_32k_re_thresh(5,:,intens_range),3)).*ones(1,2),'k',...
    'LineWidth',2)

for m = 1:8 
plot([3.875+rand(1,1)*0.25  4.875+rand(1,1)*0.25],[nanmean(fold_change_8k_re_thresh_sham(5,m,intens_range),3)...
    nanmean(fold_change_32k_re_thresh_sham(5,m,intens_range),3)],'-k.','MarkerSize',25)
end

plot([3.875 4.125],nanmean(nanmean(fold_change_8k_re_thresh_sham(5,:,intens_range),3)).*ones(1,2),'k',...
    'LineWidth',2)
plot([4.875 5.125],nanmean(nanmean(fold_change_32k_re_thresh_sham(5,:,intens_range),3)).*ones(1,2),'k',...
    'LineWidth',2)
yline([1],'k--') 
box off 
ylim([0 4]) 
xticks([1 2 4 5]) 
xticklabels({'8','32','8','32'}) 
set(gca,'xcolor','k','ycolor','k')
ylabel('Fold change near threshold, post-exposure') 
print(gcf,[savedir 'Effect_summary.pdf'],'-dpdf')
% Now do a repeated measures ANOVA 

%Now, do stats on the pre vs. 2 wk post thresholds 
clear allValues varNames dataTable
 sham_fc = [nanmean(fold_change_8k_re_thresh_sham(5,:,intens_range),3);...
     nanmean(fold_change_32k_re_thresh_sham(5,:,intens_range),3)]';
noise_fc = [nanmean(fold_change_8k_re_thresh(5,:,intens_range),3);...
    nanmean(fold_change_32k_re_thresh(5,:,intens_range),3)]';

allValues = [noise_fc;sham_fc] 
%Now we want to format this into a table with the correct labels
for i = 1:2
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:18
    if ismember(i,[1:10])
        expGroup{i,1} = 'Noise';
    else
        expGroup{i,1} = 'Sham';
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Frequency'};
freqLabels = [1:2]';
freqLabels = arrayfun(@num2str, freqLabels, 'UniformOutput', 0);
withinTable = table(freqLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V2~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Frequency');
rAnovaResults
%Also do post-hoc tests for each frequency. Do as a two sample t test
%comparing noise to sham. One-sided test, Bonferonni corrected, so net no
%effect.
allValues1 = allValues
for i = 1:2 %For 6 frequencies
    nums1 = allValues1(1:10,i); nums2 = allValues1(11:end,i);
    [h,p] = ttest2(nums1,nums2);
    disp(freqLabels{i}); disp(p); disp('');
end




%% Correlations
face_an_index = [1:3 6 8:10]; %Mice which have 2 wk data 
abr_ind_2wk = [1:3 5:8]; %Mice which have 2 wk data 

load('C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 5\abr_summary_values.mat')

%Resp amplitude re: thresh 
face_amp_re_thresh =  nanmean(fold_change_8k_re_thresh(5,:,1:2),3);
figure;
h = plot(cc.fold_change_w4_amp_nthr(abr_ind_2wk),face_amp_re_thresh(face_an_index),'ko')
box off
xlabel('Fold change W4 resp near thresh')
ylabel('Fold change facial mvmt amp near thresh') 
xlim([-.1 3.5])
ylim([-.1 3.5])
title('Near threshold amplitudes') 
%Plot best fit line
coefficients = polyfit(cc.fold_change_w4_amp_nthr(abr_ind_2wk),...
    face_amp_re_thresh(face_an_index), 1);
numFitPoints = 1000; 
xFit = linspace(min(cc.fold_change_w4_amp_nthr(abr_ind_2wk)),...
     max(cc.fold_change_w4_amp_nthr(abr_ind_2wk)), numFitPoints);
yFit = polyval(coefficients, xFit);
hold on;
plot(xFit, yFit, 'k', 'LineWidth', 1);
print(gcf,[savedir 'W4_resp_amp_re_thresh_correlation.pdf'],'-dpdf')
[r p] = corr(cc.fold_change_w4_amp_nthr(abr_ind_2wk)',face_amp_re_thresh(face_an_index)')



face_amp_re_thresh =  nanmean(fold_change_8k_re_thresh(5,:,1:2),3);
figure;
h = plot(cc.fold_change_w1_amp_nthr(abr_ind_2wk),face_amp_re_thresh(face_an_index),'ko')
box off
xlabel('Fold change w1 resp near thresh')
ylabel('Fold change facial mvmt amp near thresh') 
xlim([-.1 3.5])
ylim([-.1 3.5])
title('Near threshold amplitudes') 
%Plot best fit line
coefficients = polyfit(cc.fold_change_w1_amp_nthr(abr_ind_2wk),...
    face_amp_re_thresh(face_an_index), 1);
numFitPoints = 1000; 
xFit = linspace(min(cc.fold_change_w1_amp_nthr(abr_ind_2wk)),...
     max(cc.fold_change_w1_amp_nthr(abr_ind_2wk)), numFitPoints);
yFit = polyval(coefficients, xFit);
hold on;
plot(xFit, yFit, 'k', 'LineWidth', 1);
print(gcf,[savedir 'w1_resp_amp_re_thresh_correlation.pdf'],'-dpdf')
[r p] = corr(cc.fold_change_w1_amp_nthr(abr_ind_2wk)',face_amp_re_thresh(face_an_index)')



face_amp_re_thresh =  nanmean(fold_change_8k_re_thresh(5,:,1:2),3);
figure;
h = plot(cc.fold_change_w2_amp_nthr(abr_ind_2wk),face_amp_re_thresh(face_an_index),'ko')
box off
xlabel('Fold change w2 resp near thresh')
ylabel('Fold change facial mvmt amp near thresh') 
xlim([-.1 3.5])
ylim([-.1 3.5])
title('Near threshold amplitudes') 
%Plot best fit line
coefficients = polyfit(cc.fold_change_w2_amp_nthr(abr_ind_2wk),...
    face_amp_re_thresh(face_an_index), 1);
numFitPoints = 1000; 
xFit = linspace(min(cc.fold_change_w2_amp_nthr(abr_ind_2wk)),...
     max(cc.fold_change_w2_amp_nthr(abr_ind_2wk)), numFitPoints);
yFit = polyval(coefficients, xFit);
hold on;
plot(xFit, yFit, 'k', 'LineWidth', 1);
print(gcf,[savedir 'w2_resp_amp_re_thresh_correlation.pdf'],'-dpdf')
[r p] = corr(cc.fold_change_w2_amp_nthr(abr_ind_2wk)',face_amp_re_thresh(face_an_index)')



%% Helper functions
function h = getResponseStats(mvmt_resp)
level_vec = 30:10:100;
for m = 1:size(mvmt_resp,1)
    for i = 1:size(mvmt_resp,2)
        m 
        i
            clear temp_mat;
        temp_mat = squeeze(mvmt_resp(m,i,:,:));
     
        if sum(~isnan(temp_mat))==0
            latency(m,i) = NaN;
            sig_resp(m,i) = NaN;
            fano(m,i) = NaN;
            
        else
    
            
            %Interpolate between any NaNs
            if sum(isnan(temp_mat(:,1:500)),'all')>0
                for t = 1:size(temp_mat,1)
                    clear nan_ind
                    nan_ind = find(isnan(temp_mat(t,1:700)));
                    for n = 1:length(nan_ind)
                        if nan_ind(n) == 1
                            temp_mat(t,nan_ind(n)) = temp_mat(t,nan_ind(n)+1);
                        else
                            temp_mat(t,nan_ind(n)) = nanmean(temp_mat(t,nan_ind(n)-...
                                1:nan_ind(n)+1));
                        end
                    end
                end
            end
            
            %Find peak time
            clear peak_time peak_val resp_vec
            [peak_val peak_time] = max(nanmean(temp_mat(:,151:300)));
            peak_time = peak_time + 150;
            resp_vec = [nanmean(temp_mat(:,1:150),2) nanmean(temp_mat(:,peak_time-2:peak_time+2),2)];
            
            %Test for responses using paired t-test
            [sig_resp(m,i)] = ttest(resp_vec(:,1),(resp_vec(:,2)),'alpha',.05)
            resp_vec_mat(m,i,:,:) = resp_vec;
            fano(m,i) = nanstd(resp_vec(:,2)-resp_vec(:,1))./nanmean(resp_vec(:,2)-resp_vec(:,1))
            
            
            %Now determine time to half max
            clear xx lat_ind
            xx = nanmean(temp_mat(:,151:300)); %Sound comes on at 151
            [maxi maxi_ind] = max(xx);
            %Work backward from i and find half max  (is half max really best?)
            lat_ind = find(fliplr(xx(1:maxi_ind))<maxi/2,1); %Find first point from max below half max, note sample 1 = 0
            if sig_resp(m,i) == 1 && ~isempty(lat_ind)
                latency(m,i) = (maxi_ind-lat_ind+1) *1/150;
            elseif ~isempty(lat_ind)
                latency(m,i) = (maxi_ind-lat_ind+1) *1/150;
                sig_resp(m,i) = 0;
                fano(m,i) = NaN;
            else
                latency(m,i) = NaN;
                sig_resp(m,i) = 0;
                fano(m,i) = NaN;
            end
            resp(m,i,:) = xx;
        end
    end
    %Get threshold and eliminate spurious responses below threshold
    if sum(isnan(sig_resp(m,:)))==0
        
        if sum(sig_resp(m,:))>0 & sum(isnan(sig_resp(m,:)))==0
            %This is the simple version
            %threshold(m) = level_vec(find(sig_resp(m,:)==1,1));
            %More robust is find highest level that doesn't elicit a response
            if sum(sig_resp(m,:))==length(level_vec)
                threshold(m) = level_vec(1);
            else
                if find(sig_resp(m,:)==0,1,'last')<=8
                    if find(sig_resp(m,:)==0,1,'last')+1<=8
                        threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last')+1);
                    else
                        threshold(m) = 100;
                    end
                    
                else
                    threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last'));
                end
                sig_resp(m,1:find(sig_resp(m,:)==0,1,'last')) = 0;
                %latency(m,1:find(sig_resp(m,:)==0,1,'last')) = NaN;
                fano(m,1:find(sig_resp(m,:)==0,1,'last')) = NaN;
            end
        else
            threshold(m) = 125; %Distinguish between no response observed and a response
        end
        
    else
        threshold(m) = NaN;
    end
    
end

h.threshold = threshold;
%Threshold latency by sig_responses
h.latency = latency;

end



