%% Now get rate level functions for each 
clear all
close all

%Load all data
load('C:\Users\kx776\Dropbox\codeBase\Videography\Data\Figure 1\z_face_v_asr.mat') %CHANGE THIS PATH

cmap=cbrewer('div','RdBu',80);
cmap = flip(cmap)


mean_face_resp = squeeze(nanmean(z_face_v_asr.matrix_data.facial,3));
mean_pupil_resp = squeeze(nanmean(z_face_v_asr.matrix_data.pupil,3))
mean_jaw_resp = squeeze(nanmean(z_face_v_asr.matrix_data.jaw,3));
mean_nose_resp = squeeze(nanmean(z_face_v_asr.matrix_data.nose,3));
mean_ear_resp = squeeze(nanmean(z_face_v_asr.matrix_data.ear,3));
mean_startle_resp = squeeze(nanmean(z_face_v_asr.matrix_data.startle,3));
mean_eye_resp = squeeze(nanmean(abs(z_face_v_asr.matrix_data.eyelid),3));

%% Plot RLF for Face 

time_win = 151:300;
%For face, get RLF
for c = 1:8
    for int = 1:size(mean_face_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_face_resp(c,int,151:200));
        temp_ind = temp_ind +150; % Correct for ROI 
        face_rlf(c,int) = mean(mean_face_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

figure; 
subplot(2,3,1)
shadedErrorBar(15:10:115,mean(face_rlf),std(face_rlf)./sqrt(8))
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
title('Facial motion energy')
ylim([0 9]) 
hold on 
plot(15:10:115,mean(face_rlf),'k.')


%% Plot RLF for pupil
for c = 1:8
    for int = 1:size(mean_pupil_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_pupil_resp(c,int,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        pupil_rlf(c,int) = mean(mean_pupil_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,2)
shadedErrorBar(15:10:115,mean(pupil_rlf),std(pupil_rlf)./sqrt(8))
xlim([10 120])
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
title('Pupil') 
ylim([0 9]) 
hold on
plot(15:10:115,mean(pupil_rlf),'k.')

%% Plot RLF for jaw
%For jaw, get RLF
for c = 1:8
    for int = 1:size(mean_jaw_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_jaw_resp(c,int,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        jaw_rlf(c,int) = mean(mean_jaw_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,3) 
shadedErrorBar(15:10:115,mean(jaw_rlf),std(jaw_rlf)./sqrt(8))
hold on
plot(15:10:115,mean(jaw_rlf),'k.')
xlim([10 120])
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
title('Jaw') 
ylim([0 9]) 

%% Plot RLF for nose
%For nose, get RLF
for c = 1:8
    for int = 1:size(mean_nose_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_nose_resp(c,int,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        nose_rlf(c,int) = mean(mean_nose_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,4) 
shadedErrorBar(15:10:115,mean(nose_rlf),std(nose_rlf)./sqrt(8))
xlim([10 120])
xticks(15:20:115)
hold on
plot(15:10:115,mean(nose_rlf),'k.')
box off
xlabel('Intensity')
ylabel('Peak response') 
title('nose') 
ylim([0 9]) 


%% Plot RLF for ear 
time_win = 151:300;
%For ear, get RLF
for c = 1:8
    for int = 1:size(mean_ear_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_ear_resp(c,int,time_win));
        temp_ind = temp_ind +150; % Correct for ROI 
        ear_rlf(c,int) = mean(mean_ear_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,5)
shadedErrorBar(15:10:115,mean(ear_rlf),std(ear_rlf)./sqrt(8))
hold on
plot(15:10:115,mean(ear_rlf),'k.')
xlim([10 120])
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
title('Ear') 
ylim([0 9]) 

%% Plot RLF for startle
time_win = 10001:20000;
%For ear, get RLF
for c = 1:8
    for int = 1:size(mean_startle_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_startle_resp(c,int,time_win));
        temp_ind = temp_ind +10000; % Correct for ROI 
        startle_rlf(c,int) = mean(mean_startle_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,6)
shadedErrorBar(15:10:115,mean(startle_rlf),std(startle_rlf)./sqrt(8))
xlim([10 120])
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
%ylim([0 9]) 
title('Startle') 

%Z-scored startle 
% mean_startle_resp_z = (mean_startle_resp - mean(squeeze(mean(mean_startle_resp(:,:,1:10000),2)),2))...
%      ./ std(squeeze(mean(mean_startle_resp(:,:,1:10000),2)),[],2); 
time_win = 10001:20000;
%For ear, get RLF
for c = 1:8
    for int = 1:size(mean_startle_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_startle_resp(c,int,time_win));
        temp_ind = temp_ind +10000; % Correct for ROI 
        startle_rlf(c,int) = mean(mean_startle_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end

subplot(2,3,6)
shadedErrorBar(15:10:115,mean(startle_rlf),std(startle_rlf)./sqrt(8))
hold on
plot(15:10:115,mean(startle_rlf),'k.')
xlim([10 120])
xticks(15:20:115)
box off
xlabel('Intensity')
ylabel('Peak response') 
%ylim([0 9]) 
title('Startle') 


%% Make RLF for the eyelid
time_win = 151:300;
%For face, get RLF
for c = 1:8
    for int = 1:size(mean_eye_resp,2)
        clear temp_ind
        [dum temp_ind] = max(mean_eye_resp(c,int,151:200));
        temp_ind = temp_ind +150; % Correct for ROI 
        eye_rlf(c,int) = mean(mean_eye_resp(c,int,temp_ind-2:temp_ind+2)); 
    end
end


%% Plot all the data together 
%Order  face pupil jaw nose ear startle
%savedir = 'C:\Users\claytok\Dropbox\File_Transfer\Dan - Kameron\Staff scientist\Conference Presentations\ARO23\Facial_videography\Figure1\RLFs\'
figure;
h(1) = shadedErrorBar(15:10:115,mean(face_rlf),std(face_rlf)./sqrt(8),'LineProp','r')
hold on 
h(2) = shadedErrorBar(15:10:115,mean(pupil_rlf),std(pupil_rlf)./sqrt(8),'LineProp','b')
h(3) = shadedErrorBar(15:10:115,mean(jaw_rlf),std(jaw_rlf)./sqrt(8),'LineProp','g')
h(4) = shadedErrorBar(15:10:115,mean(nose_rlf),std(nose_rlf)./sqrt(8),'LineProp','c')
h(5) = shadedErrorBar(15:10:115,mean(ear_rlf),std(ear_rlf)./sqrt(8),'LineProp','y')
h(6) = shadedErrorBar(15:10:115,mean(eye_rlf),std(eye_rlf)./sqrt(8))
legend ([h(1:6).mainLine],{'face','pupil','jaw','nose','ear','eyelid'},'Location','Northwest')
box off 
xlabel('dB SPL')
savedir = 'C:\Users\kx776\Dropbox\codeBase\Videography\figure_making\Revision_figures\';
ylabel('Movement amplitude (z-scored)')
print(gcf,[savedir 'Video_rlf_V2.pdf'],'-dpdf')


figure;
h2 = shadedErrorBar(15:10:115,mean(startle_rlf),std(startle_rlf)./sqrt(8))
legend ([h2.mainLine],{'startle'},'Location','Northwest')
box off 
xlabel('dB SPL')
ylabel('Movement amplitude (force plate)')
%print(gcf,[savedir 'Startle_rlf_V2.pdf'],'-dpdf')





%% Now that we have everything, plot another way 

figure;
subplot(2,3,1)
imagesc(face_rlf)
colormap('hot')
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
ylabel('Mouse #') 
colorbar
title('Facial motion energy')

subplot(2,3,2)
imagesc(pupil_rlf)
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
colorbar
title('Pupil')
%caxis([0 9]) 

subplot(2,3,3)
imagesc(jaw_rlf)
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
title('Jaw')
colorbar
%caxis([0 9]) 

subplot(2,3,4)
imagesc(nose_rlf)
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
%caxis([0 9]) 
xlabel('Intensity (dB SPL)') 
ylabel('Mouse #') 
title('Nose')
colorbar

subplot(2,3,5)
imagesc(ear_rlf)
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
title('Ear')
colorbar
%caxis([0 9]) 

subplot(2,3,6)
imagesc(startle_rlf)
title('Startle')
%caxis([0 9]) 
xticks(1:2:11)
xticklabels(string(15:20:115))
colorbar
%caxis([0 9]) 






