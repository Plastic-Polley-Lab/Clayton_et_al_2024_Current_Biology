%% create representative videos
% load the avi.video
% K = VideoReader('KeC23-Session2-Run1.000DLC_resnet50_PTCHD-010421Jan4shuffle1_1030000_labeled.mp4'); 
clear
% K = VideoReader('KeC22-Session2-Run1.001.avi'); 
K = VideoReader('KeC44-Session6-Run1.001.avi'); 
eventFrame = 372 - 1;
% K = VideoReader('KeC22-Session2-Run1.001DLC_resnet101_PTCHD-010421Jan4shuffle2_1030000_labeled.mp4'); 
% eventFrame = 3836, 7967,  15528 -1; 4rd video end; 14963 starts on 7th video
load('Summary_data.mat');
load('facial_ROI.mat')
%%
% ellipse = pupil_size.pupil_param;
% pupilsize = pupil_size.pupil_area;
ellipse = data(2).pupil.pupil_param;
pupilsize = data(2).pupil.pupil_area;
% K = VideoReader('KeC23-Session2-Run1.000.avi');

%
% create AVI object
vidObj           = VideoWriter('test_forself.mp4', 'MPEG-4')                                                    ;
vidObj.Quality   =100                                                                       ;
vidObj.FrameRate =60                                                                        ;
open(vidObj)                                                                                ;

% set time for orofacial trace
% for i=1:90
%     if i == 1
%         Time(i) = 27;
%     else
%         Time(i) = Time(i-1)+5;
%     end
% end
%load ('TraceLeftRight5s.mat')
% trim the traces 
% % trim the traces 
%  traceLeft = mean(AllLeftST);
%  traceRight= mean(AllRightST);

% time=expeL.OroFacial.edges;
% % invert the sign of the traces
% traceLeft2 = 90-traceLeft; % invert the sign of the y axes for plotting into an image
% max_value = max(traceLeft2); 
% traceLeft2 = 5.2*((traceLeft2-max_value)*4+max_value);clear max_value;
% traceRight2 = 90-traceRight; % invert the sign of the y axes for plotting into an image
% max_value = max(traceRight2); 
% traceRight2 = 5.2*((traceRight2-max_value)*4+max_value);clear max_value;
numFrames = floor(K.FrameRate*K.Duration);
T = 0.030*(1:numFrames);
% create the object for animatedline
%subplot(2,1,2); 
Line = animatedline('Color','r');                                                                         
%time = 1:150;
% f = 164;
%
for f = 1:600
%for f = 1:200
    mov =     read(K,f); 
%     mov2 = rgb2gray(mov);

    imshow(mov); % adjust here for the brightness

%     set(gcf, 'Position',  [0, 0, 480, 640]);    % adjust here for the size of the video
    hold on;
    if isnan(pupilsize(f))
    else
        t = linspace(0,2*pi) ;
        x2 = ellipse{f}(1) + ellipse{f}(3)*cos(t);
        y2 = ellipse{f}(2) + ellipse{f}(4)*sin(t);
        plot(x2,y2, 'r-')
        hold on
        plot(480-pupilsize(1:f)/10, '-r', 'LineWidth', 2)
        drawnow
    end
    text(10,30,[sprintf('%.3f',T(1,f)) ' s'],'FontSize',20,'Color','w')                         ;hold on;

   writeVideo(vidObj, getframe(gca));
   

  end
    
close(vidObj)