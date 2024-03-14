function video_analysis_basedROI(path)
cd(path)
files = dir('*.avi');
ff = 1;
load('facial_ROI.mat')
for i = 1:length(files) 
    fprintf('Processing Sub Video # %d of %d\n', i, length(files))
    filename = files(i).name;
    v = VideoReader(filename);
    while hasFrame(v)
        im = readFrame(v);
        im_gray(:,:,ff) = rgb2gray(im);
        roi_posterior(:,:,ff) = imcrop(im_gray(:,:,ff), crop1);
        roi_anterior(:,:,ff) = imcrop(im_gray(:,:,ff), crop2);
        ff = ff+1;
    end
    
end

% roi_posterior = cell(1,length(files));  % high-efficient to speed up the
% processing
% tic
% parfor i = 1:length(files)
%     filename = files(i).name;
%     v = VideoReader(filename);
%     numframes = v.NumberOfFrames;
%     for ii = 1:numframes
%         im = read(v,ii);
%         im_gray = rgb2gray(im);
%         roi_posterior{i}(:,:,ii) = imcrop(im_gray, crop1);
% %         roi_anterior{i}(:,:,ii) = imcrop(im_gray, crop2);
%     end
%     i
% end
% run_time = toc
move_anterior= abs(diff(double(roi_anterior), 1, 3));
Mtrace_anterior = mean(squeeze(mean(move_anterior,1)),1);

move_posterior = abs(diff(double(roi_posterior), 1, 3));
Mtrace_posterior = mean(squeeze(mean(move_posterior,1)),1);
Orofacial.Mtrace_anterior = Mtrace_anterior;
Orofacial.Mtrace_posterior = Mtrace_posterior;
Orofacial.roi.posterior = crop1;
Orofacial.roi.anterior = crop2;
Orofacial.files = files;
animalID = split(filename, '.');
save([char(animalID(1)),'_Orofacial.mat'], 'Orofacial')
figure
imshow(im_gray(:,:,10))
hold on                                                                                 
rectangle('Position',crop1,'EdgeColor','b','LineWidth',2,'LineStyle','--'); hold on; % add the roi
rectangle('Position',crop2,'EdgeColor','r','LineWidth',2,'LineStyle','--'); hold on; % add the roi
export_fig('ROI',  '-png')
close all
clear Mtrace_posterior move_posterior Mtrace_anterior move_anterior roi_anterior roi_posterior im_gray
%%
% % create AVI object
% vidObj           = VideoWriter('Example_videos')                                                    ;
% vidObj.Quality   =100                                                                       ;
% vidObj.FrameRate =30                                                                        ;
% open(vidObj)           
% t = 0.030*(1:length(Mtrace_posterior));
% T = 1:3:200*3;
% ft = 4
% figure
% for i = ft:200
% %     subplot(2,1,1)
%     imshow(im_gray(:,:,i))
%     hold on                                                                                 ;
%     rectangle('Position',crop1,'EdgeColor','b','LineWidth',2,'LineStyle','--'); hold on; % add the roi
%     rectangle('Position',crop2,'EdgeColor','r','LineWidth',2,'LineStyle','--'); hold on; % add the roi
%     text(10,30,[sprintf('%.3f',t(1,i)) ' s'],'FontSize',20,'Color','w')                         ;hold on;
%     title('Raw Videos')
% %      drawnow update;
% %     subplot(2,2,2)
% %     imagesc(move_anterior(:,:,i)); colormap gray
% %     title('Whisker pads')
% %     subplot(2,2,3)
% %     imagesc(move_posterior(:,:,i)); colormap gray
% %     title('Face')
% %     subplot(2,1,2)
%     hold on
%     plot(T(ft:i), (80-Mtrace_anterior(ft:i))*5+60, 'r-', 'LineWidth',2 )
%     hold on
%     plot(T(ft:i), (80-Mtrace_posterior(ft:i))*5+60, 'b-', 'LineWidth',2)
%     ylabel('\Delta Pixel Intensity')
%     xlabel('Time (s)')
%     legend({'Whiskerpad', 'Face'})
%     title('Movement traces')
% %     ylim([0,20])
% %     xlim([0,T(200)])
% %     set(gcf, 'Position', [0 0 800 800]);
%     drawnow update;
%     writeVideo(vidObj, getframe(gca));
% 
% end
% close(vidObj)

% v = VideoReader('xylophone.mp4');
% figure
% imshow(im_gray(:,:,10))
% h = imrect(gca, crop2);
% % subplot(2,2,2)
% imshow(binaryImage)
% subplot(2,2,3)
% blackMaskedImage = allvideo(:,:,10);
% blackMaskedImage(~binaryImage) = 0;
% imshow(blackMaskedImage)
% subplot(2,2,4)
% imshow(imcrop(allvideo(:,:,10), crop1))