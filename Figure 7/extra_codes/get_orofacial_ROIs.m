%% load the files
clear; clc; close all
load('Z:\KeChen\Analysis\KeC35\Session 1\Archives010421\facial_ROI.mat')
file = dir('*.avi');
filename = file(1).name;
K=VideoReader(filename);
im1 = (read(K,10));
figure
imshow(im1(:,:,1));
%% get the first ROI, which is posterior to the whisker pad
% here you need to preload the facial_ROI.mat to use the same ROIs
% otherwise, you need to use h = imrect(gca); to generate new ROIs
h=imrect(gca, crop1);
%% get the coordinates
Crop=getPosition(h);
crop1(1,1)=Crop(1,1);
crop1(1,2)=Crop(1,2);
crop1(1,3)=Crop(1,3);
crop1(1,4)=Crop(1,4);
%% get the whisker pad
h1=imrect(gca, crop2);
%% get the coordinates
Crop_new=getPosition(h1);
crop2(1,1)=Crop_new(1,1);
crop2(1,2)=Crop_new(1,2);
crop2(1,3)=Crop_new(1,3);
crop2(1,4)=Crop_new(1,4);
%% save the data
save('facial_ROI.mat','crop1','crop2'); 