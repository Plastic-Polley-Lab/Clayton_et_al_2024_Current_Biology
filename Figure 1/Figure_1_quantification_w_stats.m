%% Now get rate level functions for each 
clear all
close all
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\Videography'))


%Load all data
load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\z_face_v_asr.mat') %CHANGE THIS PATH 


face_resp = z_face_v_asr.matrix_data.facial;
pupil_resp = z_face_v_asr.matrix_data.pupil;
jaw_resp = z_face_v_asr.matrix_data.jaw;
nose_resp =z_face_v_asr.matrix_data.nose;
ear_resp = z_face_v_asr.matrix_data.ear;
startle_resp = z_face_v_asr.matrix_data.startle;
eye_resp =abs(z_face_v_asr.matrix_data.eyelid); %Rectify eye data 


%% Plot mouse by mouse startle 
startle_resp_dat = getResponseStatsStartle(startle_resp);

%all other forms of response 
resp_str = {'face_resp','pupil_resp','jaw_resp','nose_resp','ear_resp',...
    'eye_resp','startle_resp'}; 
%% Get thresholds, fano factor trial by trial variance, and latencies by mouse
%Leave out face for now
for mvmt = 1:length(resp_str)-1
    eval([resp_str{mvmt} '_dat.stats = getResponseStats(' resp_str{mvmt} ')'])
    eval(['corrmat(mvmt,:,:,:) = ' resp_str{mvmt} '_dat.stats.resp_strength(:,:,:)'])
end

%Add corrmat 
corrmat(7,:,:,:) = startle_resp_dat.stats.resp_strength; 


% Save the face data 
%% Now plot for each on the same plot 
figure;
codec = [1 5 4 2 6 3 7]; 
for mo = 1:7
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,' resp_str{codec(mo)} '_dat.stats.threshold,400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],mean(' resp_str{codec(mo)} '_dat.stats.threshold).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 7.5])
xticks([1:7])
xticklabels({'Face','Ear','Nose','Pupil','Eye','Jaw','Startle'})
ylabel('Threshold (dB SPL)') 
%print(gcf,[savedir 'threshold_quantification.pdf'],'-dpdf')




%Do a one way anova to test 
threshold_stat_mat(:,1) = face_resp_dat.stats.threshold;
threshold_stat_mat(:,2) = ear_resp_dat.stats.threshold;
threshold_stat_mat(:,3) = nose_resp_dat.stats.threshold; 
threshold_stat_mat(:,4) = pupil_resp_dat.stats.threshold; 
threshold_stat_mat(:,5) = eye_resp_dat.stats.threshold; 
threshold_stat_mat(:,6) = jaw_resp_dat.stats.threshold; 
threshold_stat_mat(:,7) = startle_resp_dat.stats.threshold; 

clear allValues dataTable 
% Summarize the results in a scatter plot %Now do stats on above
allValues = [threshold_stat_mat];
%Now we want to format this into a table with the correct labels
for i = 1:7
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Movement_type'};
mvmtLabels = [1:7]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V7~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Movement_type');
rAnovaResults;
%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:6 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p); disp('');
end

%Compute holm-bonferroni correction
[dum p_ind] = sort(p); 
for i = 1:6 
    p(p_ind(i)) = p(p_ind(i))*(6-i+1)
end

clear p p_ind; 

%% Now plot latency for high intensities only 
figure;
codec = [1 5 4 2 6 3 7]; 
for mo = 1:7
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2)).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 7.5])
xticks([1:7])
xticklabels({'Face','Ear','Nose','Pupil','Eye','Jaw','Startle'})
ylabel('Latency (time to half max)') 

figure;
codec = [1 5 4 2 6 3 7]; 
for mo = [1:3 5:7]
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.latency(:,7:11),2)).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 7.5])
xticks([1:7])
xticklabels({'Face','Ear','Nose','Pupil','Eye','Jaw','Startle'})
ylabel('Latency (time to half max)') 
%print(gcf,[savedir 'latency_quantification.pdf'],'-dpdf')

latency_stat_mat(:,1) = nanmean(face_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,2) = nanmean(ear_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,3) = nanmean(nose_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,4) = nanmean(pupil_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,5) = nanmean(eye_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,6) =nanmean( jaw_resp_dat.stats.latency(:,7:11),2);
latency_stat_mat(:,7) = nanmean(startle_resp_dat.stats.latency(:,7:11),2);

clear allValues dataTable 
% Summarize the results in a scatter plot %Now do stats on above
allValues = [latency_stat_mat];
%Now we want to format this into a table with the correct labels
for i = 1:7
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Movement_type'};
mvmtLabels = [1:7]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V7~1','WithinDesign',withinTable);
[rAnovaResults_latency] = ranova(rmIHC, 'WithinModel','Movement_type');
%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:6 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p); disp('');
end

%Compute holm-bonferroni correction
[dum p_ind] = sort(p); 
for i = 1:6 
    p(p_ind(i)) = p(p_ind(i))*(6-i+1)
end

clear p p_ind; 

%% Now plot the coefficient of variation for high intensities only 
figure;
codec = [1 5 4 2 6 3 7]; 
for mo = 1:7
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.fano(:,7:11),2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.fano(:,7:11),2)).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 7.5])
xticks([1:7])
xticklabels({'Face','Ear','Nose','Pupil','Eye','Jaw','Startle'})
ylabel('Coefficient of variation (intrasubject)') 
%print(gcf,[savedir 'intrasubject_var.pdf'],'-dpdf')

fano_stat_mat(:,1) = nanmean(face_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,2) = nanmean(ear_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,3) = nanmean(nose_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,4) = nanmean(pupil_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,5) = nanmean(eye_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,6) =nanmean( jaw_resp_dat.stats.fano(:,7:11),2);
fano_stat_mat(:,7) = nanmean(startle_resp_dat.stats.fano(:,7:11),2);

clear allValues dataTable 
% Summarize the results in a scatter plot %Now do stats on above
allValues = [fano_stat_mat];
%Now we want to format this into a table with the correct labels
for i = 1:7
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Movement_type'};
mvmtLabels = [1:7]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V7~1','WithinDesign',withinTable);
[rAnovaResults_cov] = ranova(rmIHC, 'WithinModel','Movement_type');

%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:6 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p_cov(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p_cov); disp('');
end


%Compute holm-bonferroni correction
[dum p_ind] = sort(p_cov); 
for i = 1:6 
    p_cov(p_ind(i)) = p_cov(p_ind(i))*(6-i+1)
end

clear p_cov p_ind; 


%% Plot the intersubject variability using CoV 
%Strategy 1, take the trial by trial responses and get the std and mean 
for mo = 1:7
    for i = 1:11
        for t = 1:20 
            cov_inter(mo,i,t) = nanstd(corrmat(mo,:,i,t))./nanmean(corrmat(mo,:,i,t));
        end
    end
end

%Plot for high stimulus levels where we know there is a response 
figure;
codec = [1 5 4 2 6 3 7]; 
for mo = 1:7
    
   scatter(ones(1,20)*mo+.25*rand(1,20)-0.125,squeeze(nanmean(cov_inter(codec(mo),7:11,:),2)),400,'b.')
    hold on
   plot([mo-0.125 mo+0.125],nanmean(squeeze(nanmean(cov_inter(codec(mo),7:11,:),2)))*ones(1,2),'k-','LineWidth',2)
    
end
box off 
xlim([0.5 7.5])
xticks([1:7])
xticklabels({'Face','Ear','Nose','Pupil','Eye','Jaw','Startle'})
ylabel('Coefficient of variation (intersubject)') 
%print(gcf,[savedir 'intersubject_var.pdf'],'-dpdf')
%Now do stats 
cov_stat_mat = squeeze(mean(cov_inter(codec,7:11,:),2))';

clear allValues dataTable 
% Summarize the results in a scatter plot %Now do stats on above
allValues = [cov_stat_mat];
%Now we want to format this into a table with the correct labels
for i = 1:7
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
%Add in the exposure group variable

%Finally create a table that reflects the within subject factors, which
%here are frequency
factorNames = {'Movement_type'};
mvmtLabels = [1:7]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V7~1','WithinDesign',withinTable);
[rAnovaResults_cov] = ranova(rmIHC, 'WithinModel','Movement_type');
rAnovaResults
%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:6 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p_cov(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p_cov); disp('');
end

%Compute holm-bonferroni correction
[dum p_ind] = sort(p_cov); 
for i = 1:6 
    p_cov(p_ind(i)) = p_cov(p_ind(i))*(6-i+1)
end




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
                threshold(m) = level_vec(find(sig_resp(m,:)==0,1,'last'));
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