clear all 
close all


addpath(genpath('C:\Users\kx776\Dropbox\codeBase\Ephys\utils\'))
%% Db
penetration{1} = 'Ash116_080223_Penetration_5';
penetration{2} = 'Ash116_080223_Penetration_7';
penetration{3} = 'Ash116_080223_Penetration_9';
penetration{4} = 'Ash116_080323_Penetration_2';
penetration{5} = 'Ash116_080323_Penetration_4';
penetration{6} = 'Ash116_080323_Penetration_6';
penetration{7} = 'Ash116_080423_Penetration_2'; 
penetration{8} = 'Ash116_080423_Penetration_4';
penetration{9} = 'Ash117_080823_Penetration_2';
penetration{10} = 'Ash117_080923_Penetration_2';
penetration{11} = 'Ash117_080923_Penetration_4'; 
penetration{12} = 'Ash117_081023_Penetration_2';
penetration{13} = 'Ash117_081023_Penetration_4'; 
penetration{14} = 'Ash118_081423_Penetration_2';
penetration{15} = 'Ash118_081423_Penetration_4'; 
penetration{16} = 'Ash118_081523_Penetration_2'; 
penetration{17} = 'Ash118_081523_Penetration_5';
penetration{18} = 'Ash119_090723_Penetration_2';
penetration{19} = 'Ash119_090723_Penetration_4'; 

penetration = unique(penetration)

%CHANGE THIS PATH 
file_dir = 'C:\Users\kx776\Dropbox\codeBase\Videography\summary_data\PV_stim_ephys\';
%% Load in data 
all_tone = []; all_noise = []; psth_x_mw = []; p2t = []; sr = []; depth = [];
for s = [1:19];
   load([file_dir penetration{s} '.mat']) 
   all_noise = cat(1,all_noise,session_data.tone_data); 
   sr = [sr; session_data.tone_sr']; 
   psth_x_mw = cat(1,psth_x_mw,session_data.psth_x_mw); 
   p2t = [p2t; session_data.p2t']; 
   depth = [depth; session_data.depth'];
end
sr(isnan(sr)) = 0; 

rs_ind = find(p2t>.6); 

imagesc(squeeze(mean(psth_x_mw(p2t>.6,10,:,:),3)))



%Exclude units with laser contamination. 
psth_x_mw([533 525 362 361 356 273 242 237 143 142 141 140 130 122 99 18 ],:,:,:) = []; 
p2t([533 525 362 361 356 273 242 237 143 142 141 140 130 122 99 18 ]) = []; 
depth([533 525 362 361 356 273 242 237 143 142 141 140 130 122 99 18 ]) = []; 

rs_ind = find(p2t>.6);
psth_x_mw(rs_ind([249 259 260 205 96 82]),:,:,:) = []; 
p2t(rs_ind([249 259 260 205 96 82])) = []; 
depth(rs_ind([249 259 260 205 96 82])) = []; 

rs_ind = find(p2t>.6);
psth_x_mw(rs_ind([146]),:,:,:) = []; 
p2t(rs_ind([146])) = []; 
depth(rs_ind([146])) = []; 

%% Get the condition of interest 
deep_layer_rs = squeeze(mean(psth_x_mw(p2t>.6 & depth<0,10,:,:),3));
deep_layer_rs_ds  = reshape(deep_layer_rs,[373 10 1200/10])
deep_layer_rs_ds = squeeze(sum(deep_layer_rs_ds,2))*100

upper_layer_rs = squeeze(mean(psth_x_mw(p2t>.6 & depth>0,10,:,:),3));
upper_layer_rs_ds  = reshape(upper_layer_rs,[size(upper_layer_rs,1) 10 1200/10])
upper_layer_rs_ds = squeeze(sum(upper_layer_rs_ds,2))*100

fs_all_layers = squeeze(mean(psth_x_mw(p2t<.5,10,:,:),3));
fs_all_layers_ds  = reshape(fs_all_layers,[size(fs_all_layers,1) 10 1200/10])
fs_all_layers_ds = squeeze(sum(fs_all_layers_ds,2))*100


%% Plot data 
figure
subplot(3,1,2)
shadedErrorBar(1:120,mean(upper_layer_rs_ds),std(upper_layer_rs_ds)./sqrt(size(upper_layer_rs_ds,1)),'b')
ylim([0 10])
title('Upper layer RS units')
box off
xticks([0:25:120])
xticklabels(string([-500:250:500]))
subplot(3,1,3)
shadedErrorBar(1:120,mean(deep_layer_rs_ds),std(deep_layer_rs_ds)./sqrt(size(deep_layer_rs_ds,1)),'c')
ylim([0 10])
title('Deep layer RS units')
ylabel('sp/s')
xticks([0:25:120])
xticklabels(string([-500:250:500]))
xlabel('Time re: laser onset')
box off
hold on
subplot(3,1,1)
shadedErrorBar(1:120,mean(fs_all_layers_ds),std(fs_all_layers_ds)./sqrt(size(fs_all_layers_ds,1)),'r')
title('FS units')
box off
xticks([0:25:120])
xticklabels(string([-500:250:500]))
hold on
%Also plot out the laser stim 
stim_mat = zeros(25,48);
stim_mat(:,21:2:44,:) = 1;
stim_vec = reshape(stim_mat,[1200 1]);
stim_vec(stim_vec==0) = NaN; 
plot((0:1199)/10,stim_vec*40,'LineWidth',3,'Color','b')
set(gca,'fontname','Arial')
set(gca,'xcolor','k','ycolor','k')
set(gcf,'Position',[211   499   360   663])
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\Revision_figures\';
print('-vector','-dpdf',[savedir  'PV_stimulation_split_by_deep_and_upper_layers.pdf '])

