%% Analyze the ipsi vs. contra data 
clear all 
close all

mid = [142 163 164 165 166 167 169 170]; 
full_mat = nan(8,5,20,4500); 
for m = 1:8
    clear sorted_mat
    load(['C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\Visual_grating_data\visual_gratings_ChAT_' num2str(mid(m)) '.mat'])
    full_mat(m,:,1:size(sorted_mat,2),:) = sorted_mat; 
end


%% Now plot the data out 

trial_avg = nanmean(full_mat,3);
trial_avg = nanmean(trial_avg,2); 

figure;
h1 = shadedErrorBar(1:4500, squeeze(nanmean(trial_avg(:,1,:))),...
    squeeze(nanstd(trial_avg(:,1,:)))./sqrt(length(mid)),'lineProp','-k')
load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\aud_resp.mat') %Change this directory 
%Also plot the auditory data 
aud_avg = trial_avg; 
h1 = shadedErrorBar(1:4500, squeeze(nanmean(trial_avg(:,1,:))),...
    squeeze(nanstd(trial_avg(:,1,:)))./sqrt(length(mid)),'lineProp','-r')
xlim([0 400]) 
xticks([1:150:400])
xticklabels([-1:2])
ylim([-.5 4])
xlabel('Time re: grating onset (s)')
ylabel('Facial movement energy (z-score)')
%legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine h5.mainLine],{'-90 deg.','0','45' ,'90','180'})
title('Drifting grating responses') 
set(gca,'xcolor','k','ycolor','k')
set(gca,'FontName','Arial')


%% Now calculate the peak response for both spatial positions 
trial_avg = squeeze(trial_avg)
%For each condition, estimate the strength of the response
for m = 1:length(mid)
    [peak_val peak_time] = max(squeeze(trial_avg(m,150:300)));
    peak_time = peak_time +149;
    resp_strength(m) = nanmean(trial_avg(m,peak_time-2:peak_time+2));

end

figure; 
%plot(1*rand(1,8)*0.25+.875,resp_strength','.k','MarkerSize',24)
hold on 
ylim([-0.5 4]); 
plot([0.75 1],[nanmean(resp_strength) nanmean(resp_strength)],'k','LineWidth',2) 
xlim([0.5 2.5]) 
box off 
%ylabel('Facial movement amplitude')
ylabel('Facial motion energy (z-score)')
xticks([1 2])
xticklabels({'Grating','60 dB Noise'}) 
%Aud response strength:
aud_resp =  [2.1922  2.6731 3.5176 2.0511 3.3482 2.1017 1.8933 1.9059];
%plot(1*rand(1,8)*0.25+1.875,aud_resp','.k','MarkerSize',24)
plot([1 2], [resp_strength' aud_resp'],'Color',[0.5 0.5 0.5])
hold on 
ylim([-0.5 4]); 
plot([2.0 2.25],[nanmean(aud_resp) nanmean(aud_resp)],'k','LineWidth',2) 
set(gca,'xcolor','k','ycolor','k')
set(gca,'FontName','Arial')





