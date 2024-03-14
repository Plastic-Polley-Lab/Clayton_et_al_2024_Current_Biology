%% Videography figure 1: Movement plots

clear all
close all

%Load all data
%Set path here 
load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\z_face_v_asr.mat') %CHANGE THIS PATH 

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
%% Face

%Get mean and standard error of the mean 
mean_face_resp = squeeze(nanmean(z_face_v_asr.matrix_data.facial,3));

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
figure; 
imagesc(flipud(squeeze(mean(mean_face_resp))))
xlim([75 601])
xticks([76:75:601])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-5 5]) 
title('Facial Movement energy')
box off
colorbar



%% Pupil 
%Get mean and standard error of the mean 

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
mean_pupil_resp = squeeze(nanmean(z_face_v_asr.matrix_data.pupil,3));
%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_pupil_resp))))
xlim([75 601])
xticks([76:75:601])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-.5 .5]) 
box off 
colorbar
title('Pupil diameter')


%% Nose 
 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
mean_nose_resp = squeeze(nanmean(z_face_v_asr.matrix_data.nose,3));
%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_nose_resp))))
xlim([75 601])
xticks([76:75:601])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-3 3]) 
box off 
colorbar
title('Nose mvmt')


%% Ear 

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
mean_nose_resp = squeeze(nanmean(z_face_v_asr.matrix_data.ear,3));
%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_nose_resp))))
xlim([75 601])
xticks([76:75:601])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-5 5]) 
box off 
colorbar
title('Ear mvmt')

%% Jaw

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
mean_jaw_resp = squeeze(nanmean(z_face_v_asr.matrix_data.jaw,3));
%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_jaw_resp))))
xlim([75 601])
xticks([76:75:601])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-2 2]) 
box off 
colorbar
title('Jaw mvmt')


%% Startle  

 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
mean_startle_resp = squeeze(nanmean(z_face_v_asr.matrix_data.startle,3));

%Z-score startle data 
%z-score startle data
mean_startle_resp_z = (mean_startle_resp - mean(squeeze(mean(mean_startle_resp(:,:,1:10000),2)),2))...
     ./ std(squeeze(mean(mean_startle_resp(:,:,1:10000),2)),[],2); 

%Plot heat map 
figure; 
imagesc(flipud(squeeze(mean(mean_startle_resp_z))))
xlim([5001 20001])
xticks([5001:5000:20001])
xticklabels(string([-0.5:0.5:1]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticklabels(string([115:-10:15]))
colormap(cmap)
caxis([-250 250]) 
box off 
colorbar
title('Startle mvmt')










