function [Times,Event,crop,crop_norm]=preLickNoPaw(filename)            % RV


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
im1 = (read(K,1));
imshow(im1(:,:,1));

h=imrect(gca);Crop=getPosition(h);
 crop(1,1)=Crop(1,1);
 crop(1,2)=Crop(1,2);
 crop(1,3)=Crop(1,3);
 crop(1,4)=Crop(1,4);

h1=imrect(gca);Crop1=getPosition(h1);
 crop_norm(1,1)=Crop1(1,1);
 crop_norm(1,2)=Crop1(1,2);
 crop_norm(1,3)=Crop1(1,3);
 crop_norm(1,4)=Crop1(1,4);
 
 % if the Times.text is not stored in the current path, change the
 % following line.
 Times = dlmread('Times.txt');
 % StimInfo is a file used by Roberto to stored all the info of the
 % experiment. Change this line to adapt it as the user prefer.
 load('StimInfo','expe');
 % info about the timestamps of the experiment from StimInfo. See below.
 % creating the Event. mat a structure that will be used for orofacial
 % movement detection.
 
 if isfield(expe{1},'Null')==1
     Event.Null=expe{1}.Null;
     Event.Odor=expe{1}.Odor.Odor.Event;
     Event.Somato=expe{1}.Somato.Somato.Event;
     Event.Left_LED=expe{1}.Left_LED.Left_LED.Event;
     Event.Tone1=expe{1}.Tone.Tone.Event;
 else
     Event.Odor=expe{1}.Odor.Odor.Event;
     Event.Somato=expe{1}.Somato.Somato.Event;
     Event.Left_LED=expe{1}.Left_LED.Left_LED.Event;
     Event.Tone1=expe{1}.Tone.Tone.Event;
 end
 
save('Event.mat','-struct','Event'); %save the time-stamp of the event in a struct.
save('Times','Times')              ; %save the time-stamp of the video in a mat file.

expe{1}.OroFacial.roi(1,:) = Crop;  %saving the coordinates of the ROI for the orofacial movement (pixel intensity changes) detection
expe{1}.OroFacial.roi(2,:) = Crop1; %saving the coordinates of the control ROI for the orofacial movement (pixel intensity changes) detection

save('StimInfo','expe');

if isfield(expe{1},'Null')==1
    clear('Null','Odor','Somato','Tone1','Left_LED');
else
    clear('Odor','Somato','Tone1','Left_LED');
end
    



