%%%%%
% script to try to generate video movie and overlay the trace of orofacial
% movements

%%
% load the avi.video
K = VideoReader('RVKC235_053117_Trim3847-3908_BrightUp.avi'); 


%
% create AVI object
vidObj           = VideoWriter('test')                                                    ;
vidObj.Quality   =100                                                                       ;
vidObj.FrameRate =30                                                                        ;
open(vidObj)                                                                                ;

% set time for orofacial trace
for i=1:90
    if i == 1
        Time(i) = 27;
    else
        Time(i) = Time(i-1)+5;
    end
end
%load ('TraceLeftRight5s.mat')
load('RVKC235_053117.mat');
% trim the traces 
% trim the traces 
 traceLeft = mean(AllLeftST);
 traceRight= mean(AllRightST);

time=expeL.OroFacial.edges;
% invert the sign of the traces
traceLeft2 = 90-traceLeft; % invert the sign of the y axes for plotting into an image
max_value = max(traceLeft2); 
traceLeft2 = 5.2*((traceLeft2-max_value)*4+max_value);clear max_value;
traceRight2 = 90-traceRight; % invert the sign of the y axes for plotting into an image
max_value = max(traceRight2); 
traceRight2 = 5.2*((traceRight2-max_value)*4+max_value);clear max_value;
%
T = 0.030*(1:K.NumFrames);
% create the object for animatedline
%subplot(2,1,2); 
Line = animatedline('Color','r');                                                                         
%time = 1:150;
f = 164;
for f = 1:K.NumFrames
%for f = 1:200
    mov =     (read(K,f)); 
%     mov2 = rgb2gray(mov);
%     h = ones(5,5)/25                                                                        ;
%     image = imfilter(mov2, h);  
%     Image = image(:,:,1);
%     IImage = imadjust(Image);
    % filter it....
    %imshow(IImage(:,:,1),[1 300],'border','tight'); % adjust here for the brightness
    imshow(mov(:,:,1),[10 150],'border','tight'); % adjust here for the brightness

    set(gcf, 'Position',  [50, 100, 500, 500]);    % adjust here for the size of the video
    hold on                                                                                 ;
    rectangle('Position',[320,230,60,90],'EdgeColor','g','LineWidth',2,'LineStyle','--'); hold on; % add the roi
    text(10,30,[sprintf('%.3f',T(1,f)) ' s'],'FontSize',20,'Color','w')                         ;hold on;
    plot([Time(1) Time(15)],[380 380],'w','LineWidth',2); hold on;
    text(35,360,'0.5 s','FontSize',20,'Color','w')                         ;hold on;
   if f<322 
    text(280,30,'Left Trial','FontSize',25,'Color',[0.1 0.4 0.9])                         ;
   else 
    text(280,30,'Right Trial','FontSize',25,'Color',[0.6 0.1 0.1])                         ;
   end
   
   if (f>1) && (f<50) 
    rectangle('Position',[97,51,28,37],'EdgeColor','w','LineWidth',2,'LineStyle','--'); hold on; % add the roi
    text(50,111,'Control ROI','FontSize',20,'Color','w') 
   else 
   end
   
   if (f>88) && (f<106) 
    text(490,30,'Sampling','FontSize',25,'Color','w') 
   elseif (f>106) && (f<162)
    text(490,30,'Delay','FontSize',25,'Color','w') 
    elseif (f>161) && (f<170)
    text(490,30,'Decision','FontSize',25,'Color','w') 
   end
   
   if (f>410) && (f<427) 
    text(490,30,'Sampling','FontSize',25,'Color','w') 
   elseif (f>426) && (f<482)
    text(490,30,'Delay','FontSize',25,'Color','w') 
    elseif (f>481) && (f<490)
    text(490,30,'Decision','FontSize',25,'Color','w') 
   end
   
   % mark on the video when the decision (toungue touching lateral spout)
   % occur
   if (f > 161) && (f<190)
      plot([Time(62) Time(62)],[320 600],'--w','LineWidth',2);hold on;
   else
   end
   
   if (f > 480) && (f<511)
      plot([Time(62) Time(62)],[320 600],'--w','LineWidth',2);hold on;
   else
   end
  % here we need to modify it to have the trace appearing only in the right
  % time
  if (f > 100) && (f<191)
   %addpoints(Line,time(f-109),traceLeft2(1,f-109)) ;
   %plot(Time(1:f-100),5.2*traceLeft2(1:f-100),'LineWidth',5,'Color',[0.1 0.4 0.9]);
   plot(Time(1:f-100),traceLeft2(1:f-100),'LineWidth',5,'Color',[0.1 0.4 0.9]);

   drawnow 
   hold on;
   writeVideo(vidObj, getframe(gca));
   
  elseif (f > 420) && (f<511)
   %addpoints(Line,time(f-430),traceLeft2(1,f-430)) ;
   %plot(Time(1:f-420),5.2*traceRight2(1:f-420),'LineWidth',5,'Color',[0.6 0.1 0.1]);
   plot(Time(1:f-420),traceRight2(1:f-420),'LineWidth',5,'Color',[0.6 0.1 0.1]);
   drawnow 
   hold on
   %plot()
   writeVideo(vidObj, getframe(gca));
   
  else
   writeVideo(vidObj, getframe(gca));
  end
    
end

close(vidObj)



