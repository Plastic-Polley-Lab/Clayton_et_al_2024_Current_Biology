%% Analyze ABR Data 
clear all
close all

%% Plot thersholds 
%CHANGE THESE PATHS
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 5\'
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\ABR\'))
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\ephys\'))

mouse_struct{1}.thresholds.pre = [35 40 35 45 55 65]; 
mouse_struct{1}.thresholds.post = [35 45 60 95 95 100]; 
mouse_struct{1}.id = 'ChAT 82'; 
mouse_struct{1}.sessions = [2 4]; 


mouse_struct{2}.thresholds.pre = [35 35 40 45 40 75]; 
mouse_struct{2}.thresholds.post = [40 35 50 95 90 100]; 
mouse_struct{2}.id = 'ChAT 84'; 
mouse_struct{2}.sessions = [2 4]; 


mouse_struct{3}.thresholds.pre = [40 50 45 55 55 80]; 
mouse_struct{3}.thresholds.post = [35 40 60 100 90 95]; 
mouse_struct{3}.id = 'ChAT 85'; 
mouse_struct{3}.sessions = [2 4]; 


mouse_struct{4}.thresholds.pre = [30 35 35 50 45 75]; 
mouse_struct{4}.thresholds.post = [45 40 50 100 95 100]; 
mouse_struct{4}.id = 'ChAT 95'; 
mouse_struct{4}.sessions = [2 4]; 


mouse_struct{5}.thresholds.pre = [35 30 30 40 45 65];
mouse_struct{5}.thresholds.post = [40 40 65 95 95 90]; 
mouse_struct{5}.id = 'ChAT 111'; 
mouse_struct{5}.sessions = [2 4];

mouse_struct{6}.thresholds.pre = [35 35 45 40 80 95];
mouse_struct{6}.thresholds.post = [40 40 40 100 95 100]; 
mouse_struct{6}.id = 'Rbp4 3'; 
mouse_struct{6}.sessions = [2 4];

mouse_struct{7}.thresholds.pre = [30 30 40 50 50 75];
mouse_struct{7}.thresholds.post = [40 45 85 100 95 95]; 
mouse_struct{7}.id = 'ChAT 140'; 
mouse_struct{7}.sessions = [2 4];

mouse_struct{8}.thresholds.pre = [35 30 40 45 45 60];
mouse_struct{8}.thresholds.post = [40 45 80 100 95 95]; 
mouse_struct{8}.id = 'ChAT 141'; 
mouse_struct{8}.sessions = [2 4];



%% Plot threshold shift
cnt = 1;
for m = 1:length(mouse_struct)
    if ~isempty(mouse_struct{m}.thresholds.pre)
        
        threshold_shift(cnt,:) = mouse_struct{m}.thresholds.post-mouse_struct{m}.thresholds.pre
        thresholds(cnt,1,:) = mouse_struct{m}.thresholds.pre
        thresholds(cnt,2,:) = mouse_struct{m}.thresholds.post
        cnt = cnt + 1;
        
    end
end


figure;
h1 = shadedErrorBar(1:6,nanmean(thresholds(:,1,:)),nanstd(thresholds(:,1,:))./sqrt(8))
hold on
h2 = shadedErrorBar(1:6,nanmean(thresholds(:,2,:)),nanstd(thresholds(:,2,:))./sqrt(8),'r')
ylim([0 100])
xticks([1:6])
xticklabels(string([8 11.3 16 22.6 32]))
xlabel('Frequency')
ylabel('Threshold, dB SPL')
box off
xlim([1 5])
legend([h1.mainLine,h2.mainLine],{'pre','post'})
set(gca,'xcolor','k','ycolor','k')
%print(gcf,[savedir 'thresholds.pdf'],'-dpdf')

%% Compute repeated measures anova 
clear allValues dataTable 
% Summarize the results in a scatter plot %Now do stats on above
thr = reshape(thresholds(:,:,1:5),[16 5]) 
allValues = [thr];
%Now we want to format this into a table with the correct labels
for i = 1:5
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:16
    if ismember(i,[1:8])
        expGroup{i,1} = 'Pre'
    else
         expGroup{i,1} = 'Post'
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Frequency'};
mvmtLabels = [1:5]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V5~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Frequency');
rAnovaResults;


%% Get the amplitudes of waves 1-4 for 8 and 32 kHz 
 f8k_amps_w1 = nan(9,17,2);  f8k_amps_w2 = nan(9,17,2);
  f8k_amps_w3 = nan(9,17,2);  f8k_amps_w4 = nan(9,17,2);
    f8k_amps_w5 = nan(9,17,2); % f8k_amps_w4 = nan(10,5,2);
 
     f32k_amps_w1 = nan(9,17,2);  f32k_amps_w2 = nan(9,17,2);
     f32k_amps_w3 = nan(9,17,2);  f32k_amps_w4 = nan(9,17,2);
     f32k_amps_w5 = nan(9,17,2);  f8k_threshold = nan(12,2); f32k_threshold = nan(12,2); 
     %Get all data for this mouse (for all days)
     for m = 1:length(mouse_struct)
         mouse = mouse_struct{m}.id;
         sessions = mouse_struct{m}.sessions;
         
         for i = sessions
             finalPath = abrPath(mouse, i); %function abrPath needs to be set to the location of the ABR data 
             allData(i).day = i;
             allData(i).data = abrGetSessionData(finalPath);
         end
         
         %Extract the peak to peak amplitude data
         ind = find([allData(2).data.freq] == 8); %Find frequency
         levels = [20:5:100];
         cnt = 1
         for j = sessions %Number of days measured
             ind = find([allData(j).data.freq] == 8); %Find frequency
             %inds = logical(sum(allData(j).data(ind).levels == levels,2));
             %More efficient is the following
             [inds pos] = ismember(allData(j).data(ind).levels',levels);
             pos (pos==0) = [];
             f8k_amps_w1(m,pos,cnt) = abs(allData(j).data(ind).pk2pk_amp(inds));
             f8k_amps_w2(m,pos,cnt) = abs(allData(j).data(ind).w2_amp(inds));
             f8k_amps_w3(m,pos,cnt) = abs(allData(j).data(ind).w3_amp(inds));
             f8k_amps_w4(m,pos,cnt) = abs(allData(j).data(ind).w4_amp(inds));
             f8k_amps_w5(m,pos,cnt) = abs(allData(j).data(ind).w5_amp(inds));
             f8k_threshold(m,cnt) = allData(j).data(ind).threshold; 
             cnt = cnt+1;
         end
         
         %Extract the peak to peak amplitude data
         ind = find([allData(2).data.freq] == 32); %Find frequency
         levels = [20:5:100];
         cnt = 1
         for j = sessions %Number of days measured
             ind = find([allData(j).data.freq] == 32); %Find frequency
             %inds = logical(sum(allData(j).data(ind).levels == levels,2));
             %More efficient is the following
             [inds pos] = ismember(allData(j).data(ind).levels',levels);
             pos (pos==0) = [];
             if sum(inds) ==4
                 f32k_amps_w1(m,pos,cnt) = allData(j).data(ind).pk2pk_amp(inds);
                 f32k_amps_w2(m,pos,cnt) = allData(j).data(ind).w2_amp(inds);
                 f32k_amps_w3(m,pos,cnt) = allData(j).data(ind).w3_amp(inds);
                 f32k_amps_w4(m,pos,cnt) = allData(j).data(ind).w4_amp(inds);
                 f32k_amps_w5(m,pos,cnt) = allData(j).data(ind).w5_amp(inds);
                 f32k_threshold(m,cnt) = allData(j).data(ind).threshold;
                 cnt = cnt+1;
             else
                 
                 f32k_amps_w1(m,pos,cnt) = allData(j).data(ind).pk2pk_amp(inds);
                 f32k_amps_w2(m,pos,cnt) = allData(j).data(ind).w2_amp(inds);
                 f32k_amps_w3(m,pos,cnt) = allData(j).data(ind).w3_amp(inds);
                 f32k_amps_w4(m,pos,cnt) = allData(j).data(ind).w4_amp(inds);
                 f32k_amps_w5(m,pos,cnt) = allData(j).data(ind).w5_amp(inds);
                 f32k_threshold(m,cnt) = allData(j).data(ind).threshold;
                 cnt = cnt+1;
             end
         end
         
     end

   
 %% Plot narrowband noise thresholds 
 
figure; 
subplot(1,2,1) 
plot(f8k_threshold','Color',[0.5 0.5 0.5])
xticks([1 2])
xticklabels({'Pre','Post'})
box off 
xlim([0.5 2.5]) 
ylabel('Threshold') 
hold on
set(gca,'xcolor','k','ycolor','k')
plot([0.85 1],ones(1,2)*nanmean(f8k_threshold(:,1)),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(f8k_threshold(:,2)),'k','LineWidth',2)
ylim([20 100]) 
title('8 kHz NBN thresholds') 

subplot(1,2,2)
plot(f32k_threshold','Color',[0.5 0.5 0.5])
xticks([1 2])
xticklabels({'Pre','Post'})
box off 
xlim([0.5 2.5]) 
ylabel('Threshold') 
hold on
set(gca,'xcolor','k','ycolor','k')
plot([0.85 1],ones(1,2)*nanmean(f32k_threshold(:,1)),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(f32k_threshold(:,2)),'k','LineWidth',2)
ylim([20 100]) 
title('32 kHz NBN thresholds') 
%print(gcf,[savedir 'NBN_thresholds.pdf'],'-dpdf')

%Test thresholds
[h p] = ttest(f8k_threshold(:,1),f8k_threshold(:,2))
[h p] = ttest(f32k_threshold(:,1),f32k_threshold(:,2))


clear allValues dataTable varNames
% Summarize the results in a scatter plot %Now do stats on above
allValues = [f8k_threshold(1:8,1),f32k_threshold(1:8,1);...
    f8k_threshold(1:8,2),f32k_threshold(1:8,2);];
%Now we want to format this into a table with the correct labels
for i = 1:2
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable
for i = 1:16
    if ismember(i,[1:8])
        expGroup{i,1} = 'Pre'
    else
         expGroup{i,1} = 'Post'
    end
end
dataTable.expGroup = expGroup;
%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Frequency'};
mvmtLabels = [1:2]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V2~expGroup','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Frequency');
rAnovaResults;

%% Now do the same thing, but with threshold corrected responses 

% %% Now take the threshold approach 
w4_8k_re_thresh = nan(9,17,2); 
for m = 1:9
    for i = 1:2
        ind = find(levels == f8k_threshold(m,i))
        w4_8k_re_thresh(m,1:(17-ind+1),i) = ... 
            f8k_amps_w4(m,ind:end,i); 
        
    end
end

w2_8k_re_thresh = nan(9,17,2); 
for m = 1:9
    for i = 1:2
        ind = find(levels == f8k_threshold(m,i))
        w2_8k_re_thresh(m,1:(17-ind+1),i) = ... 
            f8k_amps_w2(m,ind:end,i); 
        
    end
end


w3_8k_re_thresh = nan(9,17,2); 
for m = 1:9
    for i = 1:2
        ind = find(levels == f8k_threshold(m,i))
        w3_8k_re_thresh(m,1:(17-ind+1),i) = ... 
            f8k_amps_w3(m,ind:end,i); 
        
    end
end



%Do the same thing for w1
% %% Now take the threshold approach 
w1_8k_re_thresh = nan(9,17,2); 
for m = 1:9
    for i = 1:2
        ind = find(levels == f8k_threshold(m,i))
        w1_8k_re_thresh(m,1:(17-ind+1),i) = ... 
            f8k_amps_w1(m,ind:end,i); 
        
    end
end


%% Try mean normalization approach 

w4_8k_re_thresh_norm = w4_8k_re_thresh./nanmean(w4_8k_re_thresh(:,:,1),2)

w1_8k_re_thresh_norm = w1_8k_re_thresh./nanmean(w1_8k_re_thresh(:,:,1),2)


w3_8k_re_thresh_norm = w3_8k_re_thresh./nanmean(w3_8k_re_thresh(:,:,1),2)


w2_8k_re_thresh_norm = w2_8k_re_thresh./nanmean(w2_8k_re_thresh(:,:,1),2)

%% First plot out all waves 

figure;
shadedErrorBar(0:5:30,nanmean(w4_8k_re_thresh_norm(:,1:7,1),1),...
    nanstd(w4_8k_re_thresh_norm(:,1:7,1),[],1)./sqrt(8))
hold on 
shadedErrorBar(0:5:30,nanmean(w4_8k_re_thresh_norm(:,1:7,2),1),...
    nanstd(w4_8k_re_thresh_norm(:,1:7,2),[],1)./sqrt(8),'r')
box off 
title('WIV') 
xlabel('dB SL') 
ylabel('Amplitude normalized to mean baseline amp (uV)')
print(gcf,[savedir 'WIV_amplitude.pdf'],'-dpdf')

figure;
shadedErrorBar(0:5:30,nanmean(w1_8k_re_thresh_norm(:,1:7,1),1),...
    nanstd(w1_8k_re_thresh_norm(:,1:7,1),[],1)./sqrt(8))
hold on 
shadedErrorBar(0:5:30,nanmean(w1_8k_re_thresh_norm(:,1:7,2),1),...
    nanstd(w1_8k_re_thresh_norm(:,1:7,2),[],1)./sqrt(8),'r')
title('WI') 
box off 
xlabel('dB SL') 
ylim([0 3.5]) 
ylabel('Amplitude normalized to mean baseline amp (uV)')
%print(gcf,[savedir 'WI_amplitude.pdf'],'-dpdf')


figure;
shadedErrorBar(0:5:30,nanmean(w2_8k_re_thresh_norm(:,1:7,1),1),...
    nanstd(w2_8k_re_thresh_norm(:,1:7,1),[],1)./sqrt(8))
hold on 
shadedErrorBar(0:5:30,nanmean(w2_8k_re_thresh_norm(:,1:7,2),1),...
    nanstd(w2_8k_re_thresh_norm(:,1:7,2),[],1)./sqrt(8),'r')
title('WII') 
box off 
xlabel('dB SL') 
ylim([0 3.5]) 
ylabel('Amplitude normalized to mean baseline amp (uV)')
%print(gcf,[savedir 'WII_amplitude.pdf'],'-dpdf')



%% Now do quantifaction 

figure
plot([nanmean(w4_8k_re_thresh_norm(:,5:7,1),2) nanmean(w4_8k_re_thresh_norm(:,5:7,2),2)]','Color',[0.5 0.5 0.5])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'pre','post'}) 
hold on
plot([0.85 1],ones(1,2)*nanmean(nanmean(w4_8k_re_thresh_norm(:,5:7,1),2)),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(nanmean(w4_8k_re_thresh_norm(:,5:7,2),2)),'k','LineWidth',2)
title('Change in WIV amplitude ([p = 0.02))') 
ylabel('Normalized  WIV amp') 
set(gca,'xcolor','k','ycolor','k')
box off 
%print(gcf,[savedir 'W4_amplitude_line_plot.pdf'],'-dpdf')

mean_w4_pre = nanmean(nanmean(w4_8k_re_thresh_norm(:,5:7,1),2))
mean_w4_post = nanmean(nanmean(w4_8k_re_thresh_norm(:,5:7,2),2))


figure
plot([nanmean(w1_8k_re_thresh_norm(:,5:7,1),2) nanmean(w1_8k_re_thresh_norm(:,5:7,2),2)]','Color',[0.5 0.5 0.5])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'pre','post'}) 
hold on
set(gca,'xcolor','k','ycolor','k')
plot([0.85 1],ones(1,2)*nanmean(nanmean(w1_8k_re_thresh_norm(:,5:7,1),2)),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(nanmean(w1_8k_re_thresh_norm(:,5:7,2),2)),'k','LineWidth',2)
ylim([0 4.5])
ylabel('Normalized  W1 amp') 
box off 
title('Normalized WI amplitude') 
print(gcf,[savedir 'W1_amplitude_line_plot.pdf'],'-dpdf')



mean_w1_pre = nanmean(nanmean(w1_8k_re_thresh_norm(:,5:7,1),2))
mean_w1_post = nanmean(nanmean(w1_8k_re_thresh_norm(:,5:7,2),2))

se_w1_pre = nanstd(nanmean(w1_8k_re_thresh_norm(:,5:7,1),2))./sqrt(8)
se_w1_post = nanstd(nanmean(w1_8k_re_thresh_norm(:,5:7,2),2))./sqrt(8)


figure
plot([nanmean(w2_8k_re_thresh_norm(:,5:7,1),2) nanmean(w2_8k_re_thresh_norm(:,5:7,2),2)]','Color',[0.5 0.5 0.5])
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'pre','post'}) 
hold on
set(gca,'xcolor','k','ycolor','k')
plot([0.85 1],ones(1,2)*nanmean(nanmean(w2_8k_re_thresh_norm(:,5:7,1),2)),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(nanmean(w2_8k_re_thresh_norm(:,5:7,2),2)),'k','LineWidth',2)
ylim([0 4.5])
ylabel('Normalized  WII amp') 
box off 
title('Normalized WII amplitude') 
%print(gcf,[savedir 'WII_amplitude_line_plot.pdf'],'-dpdf')


mean_w2_pre = nanmean(nanmean(w2_8k_re_thresh_norm(:,5:7,1),2))
mean_w2_post = nanmean(nanmean(w2_8k_re_thresh_norm(:,5:7,2),2))

se_w2_pre = nanstd(nanmean(w2_8k_re_thresh_norm(:,5:7,1),2))./sqrt(8)
se_w2_post = nanstd(nanmean(w2_8k_re_thresh_norm(:,5:7,2),2))./sqrt(8)



%% Summary metric to correlate 

cc.fold_change_w1_amp_nthr = (nanmean(w1_8k_re_thresh_norm(:,5:7,2),2)./nanmean(w1_8k_re_thresh_norm(:,5:7,1),2))';
cc.fold_change_w4_amp_nthr = (nanmean(w4_8k_re_thresh_norm(:,5:7,2),2)./nanmean(w4_8k_re_thresh_norm(:,5:7,1),2))';
cc.fold_change_w2_amp_nthr = (nanmean(w2_8k_re_thresh_norm(:,5:7,2),2)./nanmean(w2_8k_re_thresh_norm(:,5:7,1),2))';

disp('here') 

%save([savedir 'abr_summary_values.mat'],'cc')
