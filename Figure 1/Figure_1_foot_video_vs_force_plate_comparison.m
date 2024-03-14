%% Goals here are:
%1. Plot out foot data vs. startle pad data 
%2. Get thresholds for both 
%3. Get latencies for both 
%4. Correlate both trial vy trial 

clear all 
close all 

load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\z_face_v_dlcfoot.mat')
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\Videography\figure_making'))

foot_resp = z_face_v_dlcfoot.matrix_data.foot;
startle_resp = z_face_v_dlcfoot.matrix_data.startle;

%% Plot startle measured with force plate
 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
 cmap(cmap<0) = 0;
mean_startle_resp = squeeze(nanmean(z_face_v_dlcfoot.matrix_data.startle,3));

%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_startle_resp))))
xlim([5001 20001])
xticks([5001:5000:20001])
xticklabels(string([-0.5:0.5:1]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-.250 .250]) 
box off 
colorbar
title('Startle mvmt')
print(gcf,[savedir 'Startle_raw.pdf'],'-dpdf')

%% Plot video of foot 

mean_foot_resp = squeeze(nanmean(z_face_v_dlcfoot.matrix_data.foot,3));
%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_foot_resp))))
xlim([75 301])
xticks([76:75:451])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-3 3]) 
box off 
colorbar
title('Foot mvmt')



%% Get thresholds
foot_resp_dat.stats = getResponseStats(foot_resp) 
startle_resp_dat = getResponseStatsStartle(startle_resp);

corrmat(1,:,:,:) = foot_resp_dat.stats.resp_strength(:,:,:); 
corrmat(2,:,:,:) = startle_resp_dat.stats.resp_strength(:,:,:); 

%Load in the face data from the other session 
load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\face_dat.mat')

%% Plot thresholds/latencies 
resp_str = {'foot_resp','startle_resp'}; 

figure;
codec = [1 2 ]; 
for mo = 1:2
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,' resp_str{codec(mo)} '_dat.stats.threshold,400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(' resp_str{codec(mo)} '_dat.stats.threshold).*ones(1,2),''k-'',''LineWidth'',2)'])
end
box off 
xlim([0.5 2.5])
xticks([1:2])
xticklabels({'Foot mvmt','Startle'})
ylabel('Threshold (dB SPL)') 
set(gca,'xcolor','k','ycolor','k')
box off 
ylim([20 140])



figure;
codec = [1 2];  
for mo = 1:2
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2)).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 2.5])
xticks([1:2])
set(gca,'xcolor','k','ycolor','k')
box off 
ylim([0.01 0.06]) 
xticklabels({'Foot mvmt','Startle'})
ylabel('Latency (time to half max)') 
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 1\';
print(gcf,[savedir 'Latency.pdf'],'-dpdf')

figure;
codec = [1 2]; 
for mo = 1:2
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.fano(:,7:11),2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.fano(:,7:11),2)).*ones(1,2),''k-'',''LineWidth'',2)'])
end
box off 
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Foot mvmt','Startle'})
ylabel('Intrasubject variance')
ylim([0.1 1]) 



%% Get intrasubject variance 

for mo = 1:2
    for i = 1:11
        for t = 1:20 
            cov_inter(mo,i,t) = nanstd(corrmat(mo,:,i,t))./nanmean(corrmat(mo,:,i,t));
        end
    end
end

%Plot for high stimulus levels where we know there is a response 
figure;
codec = [1 2]; 
for mo = 1:2
    
   scatter(ones(1,20)*mo+.25*rand(1,20)-0.125,squeeze(nanmean(cov_inter(codec(mo),7:11,:),2)),400,'b.')
    hold on
   plot([mo-0.125 mo+0.125],nanmean(squeeze(nanmean(cov_inter(codec(mo),7:11,:),2)))*ones(1,2),'k-','LineWidth',2)
    
end
box off 
ylabel('Coefficient of variation (intersubject)') 
box off 
xlim([0.5 2.5])
xticks(1:2)
xticklabels({'Foot mvmt','Startle'})
ylabel('Intersubject variance (time to half max)') 
ylim([0 3.5]) 


%% Now do stats 

%Do a one way anova to test 
%threshold_stat_mat(:,1) = face_resp_dat.stats.threshold;
threshold_stat_mat(:,1) = foot_resp_dat.stats.threshold;
threshold_stat_mat(:,2) = startle_resp_dat.stats.threshold; 

%Test for difference in threshold 
[h p] = ttest(threshold_stat_mat(:,1),threshold_stat_mat(:,2))


latency_stat_mat(:,1) = nanmean(foot_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,2) = nanmean(startle_resp_dat.stats.latency(:,7:11),2); 

[h p] = ttest(latency_stat_mat(:,1),latency_stat_mat(:,2))


% Summarize the results in a scatter plot %Now do stats on above
allValues = [threshold_stat_mat];
%Now we want to format this into a table with the correct labels
for i = 1:3
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Movement_type'};
mvmtLabels = [1:3]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V3~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Movement_type');
rAnovaResults;
%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:2 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p); disp('');
end

%Compute holm-bonferroni correction
[dum p_ind] = sort(p); 
for i = 1:2
    p(p_ind(i)) = p(p_ind(i))*(2-i+1)
end

clear p p_ind; 

%% Helper functions 
function h = getResponseStats(mvmt_resp)
level_vec = 15:10:115;
for m = 1:size(mvmt_resp,1)
    for i = 1:size(mvmt_resp,2)
        clear temp_mat;
        temp_mat = squeeze(mvmt_resp(m,i,:,:));
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
        [sig_resp(m,i)] = ttest(resp_vec(:,1),(resp_vec(:,2)),'alpha',.002)
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
        else
            latency(m,i) = NaN;
            sig_resp(m,i) = 0; 
            fano(m,i) = NaN; 
        end
        resp(m,i,:) = xx;  
    end
    %Get threshold and eliminate spurious responses below threshold
    if sum(sig_resp(m,:))>0
        %This is the simple version
        %threshold(m) = level_vec(find(sig_resp(m,:)==1,1));
        %More robust is find highest level that doesn't elicit a response
        if sum(sig_resp(m,:))==length(level_vec)
            threshold(m) = level_vec(1);
        else
            if find(sig_resp(m,:)==0,1,'last')<=10
                threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last')+1);
            else
                threshold(m) = NaN; %level_vec(find(sig_resp(m,:)==0,1,'last'));
            end
            sig_resp(m,1:find(sig_resp(m,:)==0,1,'last')) = 0;
            latency(m,1:find(sig_resp(m,:)==0,1,'last')) = NaN;
            fano(m,1:find(sig_resp(m,:)==0,1,'last')) = NaN;
        end
    else
        threshold(m) = 125;
    end
   
end

h.threshold = threshold;
%Threshold latency by sig_responses 
h.latency = latency;
h.resp = resp; 
h.fano = fano; 
h.resp_vec = resp_vec_mat; 
h.resp_strength = resp_vec_mat(:,:,:,2)-resp_vec_mat(:,:,:,1);
end

% Helper function 2 get startle threshold 
function startle_resp_dat = getResponseStatsStartle(startle_resp)

level_vec = 15:10:115; 

% Threshold
% startle_p2t_noise = max(startle_resp(:,:,:,9800:10000),[],4)-... %First estimate the noise floor 
%     min(startle_resp(:,:,:,(9800:10000)),[],4);
% startle_p2t = max(startle_resp(:,:,:,10000:10200),[],4)-... %Then estimate the 
%     min(startle_resp(:,:,:,(10000:10200)),[],4);

startle_p2t_noise = max(startle_resp(:,:,:,1:10000),[],4)-... %First estimate the noise floor 
    min(startle_resp(:,:,:,(1:10000)),[],4);
startle_p2t = max(startle_resp(:,:,:,10001:20000),[],4)-... %Then estimate the 
    min(startle_resp(:,:,:,(10001:20000)),[],4);



disp('here')
for m = 1:size(startle_resp,1)
    for i = 1:size(startle_resp,2)
        [h p] = ttest(startle_p2t(m,i,:),startle_p2t_noise(m,i,:),'alpha',0.002);
        sig_resp(m,i) = h;
    end
    
    if sum(sig_resp(m,:))==length(level_vec)
        threshold(m) = level_vec(1);
    else
        if find(sig_resp(m,:)==0,1,'last')<=10
            threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last')+1);
        else
            threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last'));
        end
        sig_resp(m,1:find(sig_resp(m,:)==0,1,'last')) = 0;
    end
end


%Now get latencies (time to half maximum) 
for m = 1:size(startle_resp,1)
    for i = 1:size(startle_resp,2)
        if sig_resp(m,i)==0 
            latency(m,i) = NaN; 
        else
            clear temp_vec
            clear xx lat_ind
            xx = squeeze(nanmean(startle_resp(m,i,:,10000:10200),3)); %Sound comes on at sample 10000 
            [maxi maxi_ind] = max(xx);
            %Work backward from i and find half max  (is half max really best?)
            lat_ind = find(flipud(xx(1:maxi_ind))<maxi/2,1);
            if ~isempty(lat_ind)
                latency(m,i) = (maxi_ind-lat_ind+1) *1/10000;
            else
                disp('here')
                %Note, in some cases the peak to trough amplitude is
                %dominated by the minimum. In these cases, the latency
                %index will be empty. 
                  [mini mini_ind] = min(xx);
                  lat_ind = find(flipud(xx(1:mini_ind))>mini/2,1);
                  latency(m,i) = (mini_ind-lat_ind+1) *1/10000;
            end
        end
    end
end

%Finally get fano factors 
for m = 1:size(startle_resp,1)
    for i = 1:size(startle_resp,2)
        if sig_resp(m,i)==0 
            fano(m,i) = NaN; 
        else
            fano(m,i) = nanstd(startle_p2t(m,i,:))./nanmean(startle_p2t(m,i,:))
        end
    end
end

startle_resp_dat.stats.resp_strength = startle_p2t %- startle_p2t_noise; 
startle_resp_dat.stats.latency = latency;
startle_resp_dat.stats.threshold = threshold; 
startle_resp_dat.stats.fano = fano; 
end

