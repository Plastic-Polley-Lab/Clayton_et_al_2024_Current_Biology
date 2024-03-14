%% Plot the data for cortical silencing experiments 
clear all 
close all 

%CHANGE THESE PATHS
load('\\apollo\Polley_Lab\Mouse Videography Audiogram\Data\Data Summary\updated_laser_z_unfilt.mat')
addpath(genpath('C:\Users\kx776\Dropbox\codeBase\Videography\'))
pv_inact = squeeze(nanmean(laser_z,4));

cmap=cbrewer('div','RdBu',80);
cmap = flip(cmap)
cmap(cmap<0) = 0; 

%% Plot average RLF 
figure 
subplot(1,2,1)
imagesc(flipud(squeeze(mean(pv_inact(:,1,:,:)))))
xlim([0 120])
colormap(cmap)
caxis([-4 4])
colorbar
yticks([1:2:7]) 
yticklabels(string(95:-20:35))
xticks([1 31 61 91 121]) 
xticklabels(string(-1:1:2))
xlabel('Time re: sound onset')
ylabel('Intensity (dB SPL)') 
box off
 set(gca,'xcolor','k','ycolor','k')
title('Laser ON') 

subplot(1,2,2) 
imagesc(flipud(squeeze(mean(pv_inact(:,2,:,:)))))
xlim([0 120])
colormap(cmap)
caxis([-4 4])
colorbar
yticks([1:2:7]) 
yticklabels(string(95:-20:35))
xticks([1 31 61 91 121]) 
xticklabels(string(-1:1:2))
box off
xlabel('Time re: sound onset')
ylabel('Intensity (dB SPL)') 
title('Laser OFF') 
 set(gca,'xcolor','k','ycolor','k')
%print(gcf,[savedir 'Heat_map_laser_on_vs_off_vs_spl.pdf'],'-dpdf','-fillpage')

%% Get RLFs 
time_win = 30:60;
%For face, get Laser on RLF
for c = 1:10
    for int = 1:size(pv_inact,3)
        clear temp_ind
        [dum temp_ind] = max(pv_inact(c,1,int,time_win));
        temp_ind = temp_ind +30; % Correct for ROI 
        laser_on_rlf(c,int) = squeeze(nanmean(pv_inact(c,1,int,temp_ind-1:temp_ind+1)))% -... 
      %squeeze(nanmean(pv_inact(c,1,int,1:30))) ; 
    end
end


time_win = 30:60;
%For face, get Laser on RLF
for c = 1:10
    for int = 1:size(pv_inact,3)
        clear temp_ind
        [dum temp_ind] = max(pv_inact(c,2,int,time_win));
        temp_ind = temp_ind +30; % Correct for ROI 
        laser_off_rlf(c,int) = squeeze(nanmean(pv_inact(c,2,int,temp_ind-1:temp_ind+1))) %-... 
      %squeeze(nanmean(pv_inact(c,2,int,1:30))) ; 
    end
end


 
 %% Get laser on and off Thresholds 
 for m = 1:size(laser_z,1)
    laser_off_thresh(m,:) =  getThreshold(squeeze(laser_z(m,2,:,:,:)));
    laser_on_thresh(m,:) =  getThreshold(squeeze(laser_z(m,1,:,:,:)));
     
 end
 
 %Combine laser on and off thresholds
  thresh = laser_on_thresh + laser_off_thresh;
 %Set any responsive level to 1
 thresh (thresh>0) = 1; 
 
 %Find the first intensity with 2 consecutive response 
 for m = 1:10
     x = find(thresh(m,:)>0,2);
     if x(2)-x(1) == 1
         threshold(m) = x(1)
     else
         y = find(thresh(1,:)>0,3);
         if y(3)-y(2) == 1
             threshold(m) = y(2);
         end
     end 
 end
 
 %Now NaN out scatters 
 laser_off_rlf_thresh = laser_off_rlf; 
 laser_on_rlf_thresh = laser_on_rlf; 
 for m = 1:10
     laser_off_rlf_thresh(m,1:threshold(m)-1) = NaN; 
     laser_on_rlf_thresh(m,1:threshold(m)-1) = NaN; 
     
 end
 

%% Now make the scatters re: threshold
laser_on_re_thresh_mat = nan(length(threshold),10)
laser_off_re_thresh_mat = nan(length(threshold),10)

for i = 1:10
    laser_on_re_thresh_mat(i,1:(8-threshold(i))) = laser_on_rlf_thresh(i,threshold(i):end); 
    laser_off_re_thresh_mat(i,1:(8-threshold(i))) = laser_off_rlf_thresh(i,threshold(i):end); 
    
end

 figure; 
clear h 
cmap = colormap('hot')
cmap = cmap(20:35:180,:)

for int = 1:5
   h(int) = plot(laser_off_re_thresh_mat(:,int),laser_on_re_thresh_mat(:,int),'o','Color',cmap(int,:))
   hold on
    
end
pbaspect([1 1 1])
plot([-.5 3],[-.5 3])
legend(h,{'0','10','20','30','40'},'Location','Southeast')
xlabel('Laser off')
ylabel('Laser on')
box off 
pbaspect([1 1 1])
set(gca,'xcolor','k','ycolor','k')
%print(gcf,[savedir 'Scatter_laser_on_v_off.pdf'],'-dpdf')

 
%% Plot change as connected line plots 
for m = 1:10
    laser_off_lp(m,:) = mean(squeeze((pv_inact(m,2,threshold(m):end,1:120))))
    laser_on_lp(m,:) = mean(squeeze((pv_inact(m,1,threshold(m):end,1:120))))
    
    
    loff_at_thresh(m) =laser_off_rlf_thresh(m,threshold(m))
    loff_at_thresh_p_10(m) =laser_off_rlf_thresh(m,threshold(m)+1)
    loff_at_thresh_p_0_10(m) =mean(laser_off_rlf_thresh(m,threshold(m):threshold(m)+1))
    
    lon_at_thresh(m) =laser_on_rlf_thresh(m,threshold(m))
    lon_at_thresh_p_10(m) =laser_on_rlf_thresh(m,threshold(m)+1)
     lon_at_thresh_p_0_10(m) =mean(laser_on_rlf_thresh(m,threshold(m):threshold(m)+1))
     
end


figure; 
plot([nanmean(laser_off_rlf_thresh,2)';nanmean(laser_on_rlf_thresh,2)'],'-','Color',[0.5 0.5 0.5])
xlim([0.6 2.4]) 
xticks([1 2])
xticklabels({'Laser OFF','Laser ON'}) 
box off 
hold on
plot([0.85 1],ones(1,2)*nanmean(nanmean(laser_off_rlf_thresh,2)'),'k','LineWidth',2)
plot([2 2.15],ones(1,2)*nanmean(nanmean(laser_on_rlf_thresh,2)'),'k','LineWidth',2)
title('Per mouse (all intensities above threshold') 
ylabel('Movement amplitude (z)') 
[h p ci t]  = ttest(nanmean(laser_off_rlf_thresh,2),nanmean(laser_on_thresh,2))
%set(gcf,'Position',[2208          32         326         297])
set(gca,'xcolor','k','ycolor','k')
title('All responses above threshold') 

%% Plot average line plots 

for m = 1:10
    laser_off_lp(m,:) = mean(squeeze((pv_inact(m,2,threshold(m):end,1:120))))
    laser_on_lp(m,:) = mean(squeeze((pv_inact(m,1,threshold(m):end,1:120))))
end



laser_off_lp = laser_off_lp./max(laser_off_lp,[],2) 
laser_on_lp = laser_on_lp./max(laser_off_lp,[],2) 

figure; 
h1 = shadedErrorBar(1:120,mean(laser_off_lp),std(laser_off_lp)./sqrt(8))
hold on
h2 = shadedErrorBar(1:120,mean(laser_on_lp),std(laser_on_lp)./sqrt(8),'lineProp','b')
xticks([1 31 61 91 121]) 
xticklabels(string(-1:1:2))
xlim([0 91])
ylabel('Facial motion energy (normalized re: laser off)')
box off 
xlabel('Time (s)') 
% box off 
legend([h1.mainLine h2.mainLine],{'Laser off','Laser on'})
set(gca,'xcolor','k','ycolor','k')
%print(gcf,[savedir 'line_plot_normalized_above_threshold_resp.pdf'],'-dpdf')



%% Helper functions 
function thresh = getThreshold(dat)

for i = 1:size(dat,1)
    clear temp_ind temp_mat
    temp_mat = squeeze(dat(i,:,:));
    [x temp_ind] = max(temp_mat(:,30:60),[],2);
    temp_ind = temp_ind +30; % Correct for ROI
    [thresh(i) p] = ttest(mean(temp_mat(:,1:30),2),temp_mat(:,temp_ind:temp_ind));
end

end

