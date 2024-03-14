function Licks_No_Paw_Cineplex_Filter_Noise_Removal(filename,precue,postcue,pawfilter)
%This function was made 03/07/12 by Matt and was modified 12/08/14 by
%Roberto.

%%%%%edges is set so that edges(preframes + 1) is time zero%%%%%%
%this is the first frame following the event. When edges is used, time zero
%is the first bin looked at for a mouth movement.

%This function is an updated version of Licks_One_Tone_No_Paw_Cineplex in
%which the number of Events to analyze is not limited to just one event as
%in the One_Tone version. This script was further updated to
%Licks_No_Paw_Cineplex_Norm, which allows for a normalizing of changes in
%pixel intensity to the background. An additional input, crop_norm, has
%been added to this script. This input contains coordinates within the
%video file of an area with no motion during the recording. Any changes
%within this area are assumed to be due to changes in light or some other
%disturbance within the recording.

%This function opens a video file associated withss a recorded .plx file, and
%pulls video frames around each event time in the variable events. The
%script then determines the relative amount of rat mouth movements frame to frame by
%taking the difference in intensity of each pixel from two consecutive
%frames. 
%The crop variable provides dimensions of how to crop the video images so
%that only the mouth of the animal is visible. The times variable contains
%the timestamps of each image in the video. The Events variable is a
%structure containing the events to be analyzed. Each variable within the
%structure contains the timestamps of when the event occurs relative to the video timestamps
%in times. The filename variable contains the filename of the video to be
%analyzed

%Crop
[Times,Events,crop,crop_norm]=preLickNoPaw(filename);    
load('StimInfo','expe');
%Cineplex files have an image so crops have x's and y's swapped
%figure(2),imshow(im1(round(crop(1,2)):round(crop(1,2))+round(crop(1,4)),round(crop(1,1)):round(crop(1,1))+round(crop(1,3))))
%xright = round(crop(1,2))+round(crop(1,4));
xright = round(crop(1,2))+round(crop(1,4));
xleft = round(crop(1,2));
ytop = round(crop(1,1));
ybottom = round(crop(1,1))+round(crop(1,3));  %

%crop_norm is a crop of a blank area of the video so that background noise
%can be subtracted from the imagecrop.
xnright = round(crop_norm(1,2))+round(crop_norm(1,4));
xnleft = round(crop_norm(1,2));
yntop = round(crop_norm(1,1));
ynbottom = round(crop_norm(1,1))+round(crop_norm(1,3)); 


%pawfilter sets the value of the pixel which is considered white enough to
%pick up a paw in the image
%pawfilter = 110;
n = 0;

% save the coordinate used for crop and crop norm
expe{1}.OroFacial.crop      = [xright,xleft,ytop,ybottom];
expe{1}.OroFacial.crop_norm = [xnright,xnleft,yntop,ynbottom];

%The sampling rate is 30 hz. 

expe{1}.OroFacial.framerate = 30;

%duration is the number of seconds to do the analysis following the time
expe{1}.OroFacial.postCS = postcue;

%preduration is the time to analyze before the event
expe{1}.OroFacial.preCS = precue;

%edges provides the x-axis scale for creating a plot
expe{1}.OroFacial.edges = -1*expe{1}.OroFacial.preCS: (1/expe{1}.OroFacial.framerate): ...
    (expe{1}.OroFacial.postCS - 1/expe{1}.OroFacial.framerate);

%%%%%edges is set so that edges(preframes + 1) is time zero%%%%%%

frames    = round(expe{1}.OroFacial.postCS *expe{1}.OroFacial.framerate);
preframes = round(expe{1}.OroFacial.preCS*expe{1}.OroFacial.framerate)  ;

%This gives the total number of pixels
pixel_number = (xright - xleft)*(ybottom - ytop);  


%I is the video object which contains each frame in the video
I = VideoReader(filename);

%This section of code progresses through each Event in the event structure.
%It then cycles through each time of the individual event while 
%concurrently looping through each time in the video times
%variable. Flag1 keeps track of where to restart the times loop after
%an event is matched with the time variable. When the times loop surpasses a time
%in the events variable, the first video frame following the event can be
%pulled up at that time.

Event_Names = fieldnames(Events);

for i = 1: length(Event_Names)    

    %this checks whether the event has any trials within it
    if isempty(Events.(Event_Names{i}))
        
        %this sets Trial to an empty set and moves on to the next event
        Trial.(Event_Names{i}) = [];
        expe{1}.OroFacial.Trial.(Event_Names{i}) = [];
        Trial_no_out.(Event_Names{i}) = [];
        expe{1}.OroFacial.Trial_no_out.(Event_Names{i}) = [];
        continue    
    end
    
    
    %Flag1 is a place keeper for the loop through the frames of the video. 
    flag1 = 1;
    %This for loop iterates through each time in the events variable. 
for m = 1:length(Events.(Event_Names{i}))

    Break = 0;
    
    %This loop goes through each time in the video to find the first frame
    %following an event. Flag1 allows this loop to start where it left off
    %from the previous event. This makes the assumption that the event
    %variables are in chronological order.
    for n = flag1:length(Times)
        
        %this finds whether a video frame time is greater than the event
        %time
        if Events.(Event_Names{i})(m) <= Times(n)
            
            %This loop pulls up each frame from a time window around the
            %event. Since each movement data point is a difference in pixel
            %intensity from one frame to the next, the process starts on
            %the 2nd frame of the prefame window. This results in one less
            %datapoint than the frames + preframes total. This convention
            %was used to match the method used in the earlier analysis.
            
            %h is the size of the filter to be used. This takes an average of the 
            %values within the filter mask (size h) for an individual pixel
            h = ones(5,5)/25;
            
            %offset provides the time offset between the event and the timestamp of the
            %first image following the event
            Offset.(Event_Names{i})(m) = Times(n) - Events.(Event_Names{i})(m);
            expe{1}.OroFacial.Offset.(Event_Names{i})(m) = Times(n) - Events.(Event_Names{i})(m);
            
            for k = 1: (frames + preframes)
                
                %This loop runs through frames + preframes times (if 1 sec pre and 2 sec
                %post this is 90 times (30 hz sampling rate) through with the event occurring at 31)
                %The psth is defined by im2, So the first frame following the event should
                %be the preframes + 1 frame (or 31 in our example). This value is the first
                %frame following the event minus the frame just before the event. For this reason, im1 
                %should be taken preframes + 1 before the event. In our example where there
                %are 30 preframes and the event occurs just prior to frame 31, im1 should
                %be frame zero and im2 frame 1 giving 30 im2 frames prior to the event.
                im1 = (read(I,n - preframes - 2 + k));
                image = imfilter(double(im1(:,:,1)), h);
                %image = imrotate(image, -90);
                
                im2 = (read(I,(n - preframes - 1 + k)));
                image2 = imfilter(double(im2(:,:,1)), h);
                %image2 = imrotate(image2, -90);
                clear im1 im2
                
                %inorm and val95 find aspects of the noise which are used to filter the
                %image. The norm_Crop is used to assess difference in pixels frame to
                %frame in an area where no intensity difference should be occurring. A
                %histogram of the differences frame to frame (not absolute) should ideally
                %be centered around 0 with a normal distribution of the pixel difference frame to frame.
                %inorm is the offset of the mean of the distribution from zero. This is
                %subtracted from the difference between the two images. Once the offset is
                %removed, a %95 confidence interval of the noise is determined. This is set
                %to val95. Since any changes of mouth movement that are within the frame to
                %frame subtracted pixel noise can not be distinguished from the noise, all
                %changes in pixel intensity below the 95% confidence value (val95) are
                %removed and the difference matrix is reduced by this amount with any
                %negative values set to zero. 
                %This code finds the difference in pixel activity in the norm_crop area.
                %The difference in the norm_crop area is attributed to temporal visual noise and is
                %subtracted from the absolute difference in the output.
                ioffset = mean(mean(image2(xnleft:xnright, yntop: ynbottom) - image(xnleft:xnright, yntop: ynbottom)));
                %this creates a matrix of the difference frame to frame of the Crop_norm
                %region with the offset subtracted so as to make the mean of the
                %distribution zero. The absolute value is taken so that only one side of
                %the confidence interval needs to be taken (the noise is assumed to be
                %gaussian distributed.
                norm = abs(image2(xnleft:xnright, yntop: ynbottom) - image(xnleft:xnright, yntop: ynbottom) - ioffset);
                
                %this finds the difference value which includes 95% of the noise
                norm_sorted = sort(norm);
                
                %To make sure there are not massive fluctuations in val95, it is averaged
                %over the previous 5 trials
                if k <= 5
                vals(k) = norm_sorted(ceil(0.95*length(norm_sorted)));    
                val95 = mean(vals);
                else
                %circshift permutes the values in vals so that the oldest value moves
                %to the five position and is overwritten
                vals = circshift(vals, [1,-1]);    
                vals(5) = norm_sorted(ceil(0.95*length(norm_sorted)));   
                val95 = mean(vals);
                end
                
                %image3 stores the values of the absolute difference
                image3 = zeros((xright - xleft), (ybottom - ytop));
                
                %the following two loops iterate through each pixel in the image
                %matrix, but only for the cropped area.
                for x = 1:(xright - xleft)
                    
                    %The following if statement allows this for loop to break if a
                    %break is called in the next for loop. Break is the variable
                    %that keeps track of this
                    
                    if Break == 1
                        
                        Break = 0;
                        break
                    end
                    
                    for y = 1:(ybottom - ytop)
                        
                        %The following statements check whether the pixel and
                        %surrounding pixels are white, if so, the frame is dropped
                        %from the analysis, as a paw is probably in the frame.
                        %Pawfilter is the variable that must be set as a threshold
                        %for paw detection
                        
                        if (x > 1 && y > 1) && (((image(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                image(xleft + x - 2,ytop + y - 1,1) >= pawfilter) || ...
                                (image(xleft + x - 1,ytop + y - 1,1) && ...
                                image(xleft + x - 1,ytop + y - 2,1) >= pawfilter)) || ...
                                ((image2(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                image2(xleft + x - 2,ytop + y - 1,1) >= pawfilter) ...
                                || (image2(xleft + x - 1,ytop + y - 1,1) >= pawfilter && ...
                                image2(xleft + x - 1,ytop + y - 2,1) >= pawfilter)))
                            
                            image3 = NaN((xright - xleft), (ybottom - ytop));
                            
                            Break = 1;
                            break
                        end
                        
                        %Each change in pixel intensity is stored in image3, inorm is subtracted
                        %off here for each pixel change as it is the estimate for the temporal noise of the
                        %video.
                        image3(x,y) = abs(image2(xleft + x - 1,ytop + y - 1,1) - image(xleft + x - 1,ytop + y - 1,1) - ioffset);
                        
                        
                        
                    end
                end
                
                %the following two lines removes all values within the noise band by first
                %subtracting off the value marking the 95% CDF of the noise and setting all
                %negative values to zero.
                image3 = image3 - val95;
                
                image3(image3 < 0) = 0;
                
                %The total change in pixel intensity for a single frame is stored
                %in the output matrix trial, which contains trials in rows and
                %times of the frames in columns. This is normalized by the number of 
                %pixels so the total is an average per pixel
                Trial.(Event_Names{i})(m,k) = sum(sum(image3))/pixel_number;
                expe{1}.OroFacial.Trial.(Event_Names{i})(m,k) = sum(sum(image3))/pixel_number;
                
                clear image2 image image3;
            end
            
            %Flag1 needs to be set so that the loop doesn't start all the
            %way back at the beginning of the times variable.
            flag1 = n;
            
            %This ends the n for loop since the event was found
            break
            
        end
    end
    
end

%The following lines remove any values in trials lying three standard
%deviations from the mean by replacing the value with NaN and output this
%matrix as trial_no_out
trial_list = expe{1}.OroFacial.Trial.(Event_Names{i})(:);
meantrial = nanmean(trial_list);
stdtrial = nanstd(trial_list);
thresh_high = meantrial + 1*stdtrial;
[b1, b2] = find(expe{1}.OroFacial.Trial.(Event_Names{i}) >= thresh_high);
expe{1}.OroFacial.Trial_no_out.(Event_Names{i}) = expe{1}.OroFacial.Trial.(Event_Names{i});

for j = 1 : length(b1)
    expe{1}.OroFacial.Trial_no_out.(Event_Names{i})(b1(j), b2(j)) = NaN;
end

clear trial_list meantrial stdtrial thresh_high b1 b2

end

clear I

save('StimInfo','expe');
%Mouth_Move_Plot;
end


