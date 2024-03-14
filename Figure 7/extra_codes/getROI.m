function [crop1,crop2]=getROI(filename)            % RV


% This function will wait until the user draw 2 rectangle (Region Of
% Interest ROI), one of which have to be centered on the rat's mouth and
% the other in a corner of the original image for background correction.

%for clarification of the function ask Roberto.

% input  = filename; is the name of the .AVI file.
% output = times; is the timestamps of CinePlex File.
%          crop; coordinates of the ROI around the Rat's mouth.
%          crop_norm; coordinates of a ROI used for background correction.
% These outputs will feed the "Licks_No_Paw_Cineplex_Filter_Noise_Removal"
% function.

K=VideoReader(filename);
im1 = (read(K,10));
imshow(im1(:,:,1));

h=imrect(gca);Crop=getPosition(h);
 crop1(1,1)=Crop(1,1);
 crop1(1,2)=Crop(1,2);
 crop1(1,3)=Crop(1,3);
 crop1(1,4)=Crop(1,4);

h1=imrect(gca);Crop_new=getPosition(h1);
 crop2(1,1)=Crop_new(1,1);
 crop2(1,2)=Crop_new(1,2);
 crop2(1,3)=Crop_new(1,3);
 crop2(1,4)=Crop_new(1,4);
 
save('facial_ROI.mat','crop1','crop2'); %save the time-stamp of the event in a struct.



% 
% if isfield(expe{1},'Null')==1
%     clear('Null','Odor','Somato','Tone1','Left_LED');
% else
%     clear('Odor','Somato','Tone1','Left_LED');
% end
%     



