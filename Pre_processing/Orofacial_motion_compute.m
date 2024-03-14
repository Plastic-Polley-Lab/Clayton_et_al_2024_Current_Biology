%%
%Updated 7/14/22 to include comments for takesian lab
%this script takes all of the videos in one run of one session, then
%calculates the absolute difference between every frame of every video,
%utilizing the roi from first video of data set, set by user 

function get_orofacial_ROIs_compute_Motion_ChatMice_startle(mouse_num, session_num, run_num)
%% load the files-- manually for every mouse
% clear; clc; close all

%Specify directories to use
%flipped:
% data_root = 'C:\MouseVideoAudiogram\Data';
% animalID = sprintf('Chat%02d',mouse_num);
% data_dir = fullfile(sprintf('C:/MouseVideoAudiogram/Data/Chat%02d/Session %d/flipped_avi',mouse_num,session_num)); %location of the data (videos)
% aviFolder = fullfile(sprintf('C:/MouseVideoAudiogram/Data/Chat%02d/Session %d/flipped_avi',mouse_num,session_num)); %location of videos
% savingDir = fullfile(data_root,sprintf('Chat%02d/proc_orofacial/Session %d/Run%d',mouse_num,session_num,run_num)); %saving location for processed videos
%chat
% data_root = '\\apollo\Polley_Lab\Mouse Videography Audiogram\Data';
% animalID = sprintf('Chat%02d',mouse_num);
% data_dir = fullfile(data_root,sprintf('Chat%02d/Session %d',mouse_num,session_num)); %location of the data (videos)
% aviFolder = fullfile(data_root,sprintf('Chat%02d/Session %d',mouse_num,session_num)); %location of videos
% savingDir = fullfile(data_root,sprintf('Chat%02d/proc_orofacial/Session %d/Run%d',mouse_num,session_num,run_num)); %saving location for processed videos
%KS
data_root = 'D:\Videography\Data';
animalID = sprintf('ChAT_%03d',mouse_num);
data_dir = fullfile(data_root,sprintf('ChAT_%03d/Session %d',mouse_num,session_num)); %location of the data (videos)
aviFolder = fullfile(data_root,sprintf('ChAT_%03d/Session %d',mouse_num,session_num)); %location of videos
savingDir = fullfile(data_root,sprintf('ChAT_%03d/proc_orofacial/Session %d/Run%d',mouse_num,session_num,run_num)); %saving location for processed videos

cd(aviFolder);
mkdir(savingDir);
%% Selecting the ROI-- manually for the first video, then apply to all

%Naming and reading the video file
file = dir('*.avi');
%Filter files further
flist = string(vertcat(file(:).name));
flist = flist(contains(flist,['Run' num2str(run_num)])) 

filename = flist(1);
K=VideoReader(filename);
im1 = (read(K,K.NumFrames));
figure
imshow(im1(:,:,1));
% get the first ROI, which is posterior to the whisker pad
% here you need to preload the facial_ROI.mat to reuse the same ROIs
% use h = imrect(gca); to generate new ROIs and drag
% and place the rectangle
% h=imrect(gca, crop1);
h=imrect(gca);

% get the coordinates
Crop=getPosition(h);
crop1(1,1)=Crop(1,1);
crop1(1,2)=Crop(1,2);
crop1(1,3)=Crop(1,3);
crop1(1,4)=Crop(1,4);

% % get the whisker pad
% % h1=imrect(gca, crop2);
% h1=imrect(gca);
%
% % get the coordinates
% Crop_new=getPosition(h1);
% crop2(1,1)=Crop_new(1,1);
% crop2(1,2)=Crop_new(1,2);
% crop2(1,3)=Crop_new(1,3);
% crop2(1,4)=Crop_new(1,4);
%% save the data if running for the first time-- only can save one ROI

cd(savingDir);
save('facial_ROI.mat','crop1');
% save('facial_ROI_fullface.mat','crop1');
figure
im1_gray = rgb2gray(im1);
imshow(im1_gray);
hold on
rectangle('Position',crop1,'EdgeColor','b','LineWidth',2,'LineStyle','--'); hold on; % add the roi
% rectangle('Position',crop2,'EdgeColor','r','LineWidth',2,'LineStyle','--'); hold on; % add the roi
% export_fig('ROI',  '-png')
% export_fig('ROI_fullface',  '-png')

%% Applying the crop (crop1) to the rest of the videos and their frames

cd(aviFolder);
files = dir('*.avi'); %all of the files in one session
%Filter files further
flist = string(vertcat(file(:).name));
flist = flist(contains(flist,['Run' num2str(run_num)])) 

flist(1) = []; %000 trial, do not include in data
% files(end) = [];
% files(error_trials) = [];


roi = cell(1,length(flist));

% test for bad videos :
% for i = 1:length(files)
%     filename = files(i).name; 
%     v = VideoReader(filename);
%     numframes = v.NumFrames;
%     disp(numframes)
%     i
% end

tic %starts elapsed timer

%parallelize reading the video files
parfor i = 1:length(flist)
    filename = flist{i};
    v = VideoReader(filename);
    numframes = v.NumFrames;
    for ii = 1:numframes
        im = read(v,ii);
        im_gray = rgb2gray(im);
        roi{i}(:,:,ii) = imcrop(im_gray, crop1); %cropping to the ROI manually created in the first run
        %roi_anterior{i}(:,:,ii) = imcrop(im_gray, crop2);
    end
    i
end

toc
% save('roi_face_raw.mat','roi');
%%
% roi_face = [];
% for i = 1:length(roi)
%     roi_face = cat(3,roi_face,roi{i});
%     i
% end
%
% move_face = abs(diff(double(roi_face), 1, 3));
% Mtrace_face = mean(squeeze(mean(move_face,1)),1);

%% compute diff between frames without concatenating everything
tic
mean_diff_face = cell(1,length(roi)); %initialize a space for mean diff

for i = 1:length(roi)
    diff_face = abs(diff(double(roi{i}), 1, 3));
    if (i < length(roi))
        left = double(roi{i}(:,:,end)); %last frame
        right = double(roi{i+1}(:,:,1)); %first frame of next video
        inbetween_diff = abs(right - left);
        diff_face = cat(3,diff_face,inbetween_diff); %concat with all differences from before
    end
    mean_diff_face{i} = mean(squeeze(mean(diff_face,1)),1); %taking the average of the concatonated differences
    i
end
toc

Mtrace_face = []; %initialize a space for mean diff for every file
for i = 1:length(roi)
    Mtrace_face = [Mtrace_face,mean_diff_face{i}];
end
%% Saving all of the calculated differences

Mtrace_face = [0 Mtrace_face]; %adding zero to the start of the array
Mtrace_face_split = cell(1,length(files));
framesinAvi = [];
for i = 1:length(files)
    framesinAvi = [framesinAvi size(roi{i},3)];
end
cumframesAvi = [0 cumsum(framesinAvi)];
for i = 1:length(files)
    Mtrace_face_split{i} = Mtrace_face(cumframesAvi(i)+(1:framesinAvi(i)));
end

% % %if appending to an already saved variable:
% cd(savingDir);
% load([animalID,'_Orofacial.mat']);

Orofacial.Mtrace_face = Mtrace_face;
Orofacial.Mtrace_face_split = Mtrace_face_split;
Orofacial.roi.face = crop1;

% Orofacial.Mtrace_face_full = Mtrace_face;
% Orofacial.Mtrace_face_full_split = Mtrace_face_split;
% Orofacial.roi.facefull = crop1;

cd(savingDir);
save([animalID,'_Run' num2str(run_num) '_Orofacial.mat'], 'Orofacial')
end
