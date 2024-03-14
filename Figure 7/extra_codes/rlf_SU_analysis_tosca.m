function [spls,rlf, spont, psth_summary] = rlf_SU_analysis_tosca(neuron_num, spikedata, ana_window)
% neuron_num = 1; % 20, 25, 33, 35
psth = spikedata(neuron_num).psth;
% rep  = spikedata(1).stimData.numReps;
time = psth.timepoint;
%% get the spontaneous firing rate
% from the impale event to the stimulus onset(delay)
base_time = find(time<0);
spont = mean(sum(psth.scmatrix(:, base_time),2)/length(base_time) * 1000); % firingr rate Hz
%% by default there are 60 trials
% figure
% psth_plot(psth,1, 1:750)
%%
spls = spikedata(neuron_num).stimData.StimulusData.Level; % inner seq
n_spls = length(unique(spls));
spls_id = unique(spls);

% innerIndexs = spikedata(neuron_num).clusterData.stimData.innerIndexes;
n_rep = length(spls)/ n_spls  ;
if mod(length(spls), n_spls) == 0
else
    warning('There is an error, trial numbers not matched')
end
%%
n_rep = floor(n_rep);
rlf = zeros(n_spls, n_rep);
% figure
k = 1;
for j = 1:n_spls
    idx = find(spls == spls_id(j));
    idx = idx(1:n_rep);
    scmatrix = psth.scmatrix(idx,:);
    %         subplot(n_freqs,n_spls,k)
    %         fr_avg = mean(scmatrix);
    %         psth_avg_smooth = smoothts(fr_avg,'b',10);
    %         plot(psth_avg_smooth/0.001)
    %         ylim([0,500])
    %         xlim([0, 100])
    %         k = k+1
    %         title([num2str(freqs(i)), '-', num2str(spls(j))])
%     scmatrix_sum = sum(scmatrix);
%     delay = psth.stimulus.delay;
%     duration = psth.stimulus.width;
    analysis_window = find(time > (ana_window(1)-1)/1000 & time < (ana_window(end)/1000));
    rlf(j,:) = sum(scmatrix(:,analysis_window),2)./(ana_window(end)-ana_window(1) + 1) * 1000; % convert to FR (Hz)
end


%% plot the rate level function
fig =0;
if fig ==1
    subplot(1,2,1)
    data = rlf;
    err = std(data,0,2)/sqrt(size(data,2)); % standard error mean
    CT=cbrewer('div', 'RdYlBu', 6); % for nice color
    e = errorbar(mean(data,2),err, '-s','MarkerSize',6,'MarkerFaceColor', CT(6,:), 'Color',CT(6,:), 'LineWidth',1);
    hold on
    xticks([1:1:length(spls)])
    for i = 1:length(spls)
        labels_spls{i} = num2str(spls(i));
    end
    xticklabels(labels_spls(1:1:end))
    ylabel('Spike Rate (Hz)')
    xlabel('Level (dB SPL)')
    box off
    set(gca,'TickDir','out')
    % set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    set(gcf,'position',[100,200,400,300])
    hold off
end
%% plot example psth
% %%
spike = spikedata(neuron_num).psth;

for j = 1:n_spls % there are 2 different patterns
    idx = find(spls == spls_id(j));
    psth_summary.(['SPL', num2str(spls_id(j))]).raster = spike.spikeraster(idx);
    psth_summary.(['SPL', num2str(spls_id(j))]).scmatrix = spike.scmatrix(idx, :);
    psth_summary.(['SPL', num2str(spls_id(j))]).time = spike.timepoint;
end
if fig ==1
    subplot(1,2,2)
    CT=cbrewer('seq', 'PuRd', n_spls); % for nice color
    %     CT=cbrewer
    hold on
    for i = 1:n_spls
        %         idx_start = 22 + (i-1) * 100; % 100 interval for each pattern
        idx_start = 1 + (i-1) * n_rep; % 100 interval for each pattern
        hr(i) = rectangle('Position',[0,idx_start, 0.05,n_rep],'FaceColor', CT(i,:), 'EdgeColor', CT(i,:));
        hold on
        for j = 1:length(psth_summary.(['SPL', num2str(spls_id(i))]).raster)
            if ~isempty(psth_summary.(['SPL', num2str(spls_id(i))]).raster(j).times)
                %                 scatter(psth_summary.(['SPL', num2str(spls(i))]).raster.ts, j*ones(size(psth_summary.(['SPL', num2str(spls(i))]).raster).ts)), 6, '.','k')
                for k = 1:length(psth_summary.(['SPL', num2str(spls_id(i))]).raster(j).times)
                    plot([psth_summary.(['SPL', num2str(spls_id(i))]).raster(j).times(k), psth_summary.(['SPL', num2str(spls_id(i))]).raster(j).times(k)], [j-1,j] +n_rep *(i-1), 'k', 'LineWidth',1)
                end
            end
        end
    end
    
    % imagesc(spike.scmatrix),colormap(BW)
%     xlim([window(1),window(end)])
    ylim([0,n_spls * n_rep])
    axis on
    ylabel('Trial #')
    xlabel('Time (ms)')
    set(gca,'TickDir','out')
    % set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    set(gcf,'position',[100,100,800,400])
%         legend([hr(1), hr(2)], {'reg', 'rand'})
end
end