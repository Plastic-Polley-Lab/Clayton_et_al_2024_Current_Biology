%% Gaps in noise final analysis 

clear all
close all
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\Videography\'))
files = {'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_150_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_151_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_152_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_153_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_154_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_155_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_156_gaps_in_noise.mat',...
    'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\gaps_in_noise\ChAT_157_gaps_in_noise.mat'};

files = unique(files); %ensure no duplicate files get analyzed 

%% Load in data

sprate_mat = [ ];
for m = [1:8]% length(files)
    clear dat
    load(files{m})
    gap_mat(m,:,:,:) = dat.sorted_face_stim;
end

%% Investigate trial by trial noise 

ggg = reshape(gap_mat,[8 13*20 300])
figure
cnt = 1; 
for m = 1:8
   subplot(2,4,cnt)
   histogram(std(squeeze(ggg(m,:,1:150))'),100)
   cnt = cnt + 1; 
   std_vals = std(squeeze(ggg(m,:,1:150))'); 
   ggg_star = squeeze(ggg(m,:,:)); 
   ggg_star(std_vals>1.5,:) = NaN;  
   ggg(m,:,:) = ggg_star;
end

gap_mat = reshape(ggg,[8 13 20 300]); 


%% Average over trials after zero centering 
mean_gaps = squeeze(nanmean(gap_mat,3)); 

%% Plot average response
 figure;
  cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
imagesc(flipud(squeeze(mean(mean_gaps))))
xlabel('Time re: gap termination (s)')
xticks([76:75:300])
xticklabels(string([-.5:.5:.5]))
yticks([1:13])
yticklabels(string(fliplr([0  30  40  50    60    70    80    90   100   150   200   250   500])))
ylabel('Gap Duration') 
colormap(cmap)
caxis([-1.5 1.5]) 
box off 
colorbar
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\f2_elements\';
% print(gcf,[savedir 'Gaps_mean_resp.pdf'],'-dpdf')


%Plot mouse by mouse responses 
for m = 1:8 
 figure;
  cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
imagesc(flipud(squeeze(mean_gaps(m,:,:))))
xlabel('Time re: gap termination (s)')
xticks([76:75:300])
xticklabels(string([-.5:.5:.5]))
yticks([1:13])
yticklabels(string(fliplr([0  30  40  50    60    70    80    90   100   150   200   250   500])))
ylabel('Gap Duration') 
colormap(cmap)
caxis([-1.5 1.5]) 
box off 
colorbar
title(string(m))
end 
%savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\f2_elements\';
% print(gcf,[savedir 'Gaps_mean_resp.pdf'],'-dpdf')


%% Now plot growth functions 

for m = 1:size(mean_gaps,1)
    for g = 1:size(mean_gaps,2)
        clear temp_vec
        temp_vec = squeeze(mean_gaps(m,g,:)); 
        [x i] = max(temp_vec(151:200)); 
        gap_resp(m,g) = mean(temp_vec(i+150-2:i+150+2))-mean(temp_vec(1:150)); 
    end
end 


figure
plot(mean(gap_resp),'k','LineWidth',2)
hold on 
plot(gap_resp','Color',[.5 .5 .5])
xticks([1:13])
ylim([-0.2 2])
xticklabels(string([0    30    40    50    60    70    80    90   100   150   200   250   500]))
xlabel('Gap Duration')
ylabel('mean resp') 
box('off') 
% print(gcf,[savedir 'Gaps_rlf.pdf'],'-dpdf')



allValues = [gap_resp];
%Now we want to format this into a table with the correct labels
for i = 1:13
    v = strcat('V',num2str(i));
    varNames{i,1} = v;
end
dataTable = array2table(allValues, 'VariableNames', varNames);
factorNames = {'Gap_duration'};
mvmtLabels = [1:13]';
mvmtLabels = arrayfun(@num2str, mvmtLabels, 'UniformOutput', 0);
withinTable = table(mvmtLabels,'VariableNames',factorNames);
%Now that we have the tables set up, run the anova model
rmIHC = fitrm(dataTable, 'V1-V13~1','WithinDesign',withinTable);
[rAnovaResults] = ranova(rmIHC, 'WithinModel','Gap_duration');
rAnovaResults
%Also do post-hoc tests for each movement compared to face.
%One-sided test, Bonferonni corrected, so net no
%effect.
for i = 1:12 %For 6 frequencies
    nums1 = allValues(1:8,1); nums2 = allValues(1:8,i+1);
    [h,p(i)] = ttest(nums1,nums2);
    disp(mvmtLabels{i}); disp(p); disp('');
end

%% Quantify gap detection thresholds and latencies 
gap_resp_dat.stats = getResponseStats(gap_mat);
resp_str = {'gap_resp'}; 

figure;
subplot(1,2,1)
codec = [1 2 3]; 
for mo = 1
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,' resp_str{codec(mo)} '_dat.stats.threshold,400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],mean(' resp_str{codec(mo)} '_dat.stats.threshold).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 3.5])
xticks([1])
xlim([0 2]) 
xticklabels({'Face'})
ylabel('Threshold (Gap duration in ms)') 
title('Threshold') 

%Get the mean and median gap response 




subplot(1,2,2)
codec = [1 2 3]; 
for mo = 1
    
    eval(['scatter(ones(1,8)*mo+.25*rand(1,8)-0.125,nanmean(' resp_str{codec(mo)} '_dat.stats.latency,2),400,''b.'')'])
    hold on
    eval(['plot([mo-0.125 mo+0.125],nanmean(nanmean(' resp_str{codec(mo)} '_dat.stats.latency,2)).*ones(1,2),''k-'',''LineWidth'',2)'])
    
end
box off 
xlim([0.5 3.5])
xticks([1])
xlim([0 2]) 
xticklabels({'Face'})
ylabel('Latency re: gap termination (s)') 
title('Latency')
ylim([0 0.5]) 
%print(gcf,[savedir 'Gaps_threshold.pdf'],'-dpdf')



%% Helper functions 
function h = getResponseStats(mvmt_resp)
level_vec = [0  30  40 50 60 70 80 90 100 150 200 250 500];
for m = 1:size(mvmt_resp,1)
    for i = 1:size(mvmt_resp,2)
        clear temp_mat;
        temp_mat = squeeze(mvmt_resp(m,i,:,:));
        
        %Find peak time
%         clear peak_time peak_val resp_vec
        [peak_val peak_time] = max(nanmean(temp_mat(:,151:200)));
        peak_time = peak_time + 150;
         resp_vec = [nanmean(temp_mat(:,1:150),2) nanmean(temp_mat(:,peak_time-2:peak_time+2),2)];
        
%          
%         %Test for responses using paired t-test 
%         [sig_resp(m,i)] = ttest(resp_vec(:,1),(resp_vec(:,2)),'alpha',.05,'tail','left')
%         resp_vec_mat(m,i,:,:) = resp_vec;
%         fano(m,i) = nanstd(resp_vec(:,2)-resp_vec(:,1))./nanmean(resp_vec(:,2)-resp_vec(:,1))
        
        
        [sig_resp(m,i) p_resp(m,i)] = ttest2(nanmean(temp_mat(:,121:150)),nanmean(temp_mat(:,peak_time:peak_time+30)),'alpha',.002)


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
        threshold(m) = 115;
    end

end

h.threshold = threshold;
%Threshold latency by sig_responses 
h.latency = latency;
h.resp = resp; 
%h.fano = fano; 
%h.resp_vec = resp_vec_mat; 
%h.resp_strength = resp_vec_mat(:,:,:,2)-resp_vec_mat(:,:,:,1);
end




