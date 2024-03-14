%% Figure 4C 
% Plots response heat maps for pure tones vs. broadband noise 

clear all
close all
 cmap=cbrewer('div','RdBu',80);
 cmap = flip(cmap)
 cmap(cmap<0) = 0; 

%% Load data 
%CHANGE THIS PATH 
load('\\Apollo\Polley_lab\Mouse Videography Audiogram\Data\Data Summary\z_spectral_stim_unfilt.mat')

% Get mean responses
wn_mean_resp = squeeze(nanmean(z_spectral_stim.data.wn,3)); 

tones_mean_resp = squeeze(nanmean(z_spectral_stim.data.tones,4)); 

savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\summary_data\figure 5\';
%% For each stimulus plot heat map of response 

figure; 
imagesc(flipud(squeeze(mean(wn_mean_resp))))
xlim([75 301])
xticks([76:75:451])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticks([1:4]) 
yticklabels(string([90:-20:30]))
colormap(cmap)
caxis([-3 3]) 
box off 
colorbar
title('White noise')

%% Tones 
%8 kHz tones response 
figure; 
imagesc(flipud(squeeze(mean(mean((tones_mean_resp(:,:,:,:)))))))
xlim([75 301])
xticks([76:75:451])
xticklabels(string([-0.5:0.5:3]))
xlabel('Time re: sound onset')
ylabel('dB SPL')
yticks([1:4]) 
yticklabels(string([90:-20:30]))
colormap(cmap)
caxis([-3 3]) 
box off 
colorbar
title('Tone') 


