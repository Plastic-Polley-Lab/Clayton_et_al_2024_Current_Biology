function [spls,rlf, spont, psth_summary] = rlf_SU_analysis(neuron_num, spikedata, ana_window)
% neuron_num = 1; % 20, 25, 33, 35
psth = spikedata(neuron_num).clusterData.psth;
rep  = spikedata(1).clusterData.stimData.numReps;
plot_window = 1:750;
%% get the spontaneous firing rate
% from the impale event to the stimulus onset(delay)
spont = mean(sum(psth.scmatrix(:, 1:psth.stimulus.delay),2)/psth.stimulus.delay * 1000); % firingr rate Hz
%% by default there are 60 trials
% figure
% psth_plot(psth,1, 1:750)
%%
spls = spikedata(neuron_num).clusterData.stimData.inner_variables; % inner seq
n_spls = length(spls);

innerIndexs = spikedata(neuron_num).clusterData.stimData.innerIndexes;
n_rep = spikedata(1).clusterData.stimData.numTrials/ n_spls  ;
%%
rlf = zeros(n_spls, n_rep);
% figure
k = 1
for j = 1:n_spls
    idx = find(innerIndexs == j);
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
    delay = psth.stimulus.delay;
    duration = psth.stimulus.width;
    analysis_window = ana_window + delay;
    rlf(j,:) = sum(scmatrix(:,analysis_window),2)./(ana_window(end)-ana_window(1) + 1) * 1000; % convert to FR (Hz)
end


%% plot the rate level function
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
%% plot example psth
% %%
subplot(1,2,2)
spike = spikedata(neuron_num).clusterData.psth;

for j = 1:n_spls % there are 2 different patterns
    idx = find(innerIndexs == j);
    psth_summary.(['SPL', num2str(spls(j))]).raster = spike.raster(idx);
    psth_summary.(['SPL', num2str(spls(j))]).scmatrix = spike.scmatrix(idx, :);
    psth_summary.(['SPL', num2str(spls(j))]).stimulus = spike.stimulus;
end
window = plot_window;
% if fig ==1
    
    CT=cbrewer('seq', 'PuRd', n_spls); % for nice color
    %     CT=cbrewer
    hold on
    for i = 1:n_spls
        %         idx_start = 22 + (i-1) * 100; % 100 interval for each pattern
        idx_start = 1 + (i-1) * rep; % 100 interval for each pattern
        hr(i) = rectangle('Position',[spike.stimulus.delay(1),idx_start, spike.stimulus.width(1),rep],'FaceColor', CT(i,:), 'EdgeColor', CT(i,:));
        hold on
        for j = 1:length(psth_summary.(['SPL', num2str(spls(i))]).raster)
            if ~isempty(psth_summary.(['SPL', num2str(spls(i))]).raster(j).ts)
                %                 scatter(psth_summary.(['SPL', num2str(spls(i))]).raster.ts, j*ones(size(psth_summary.(['SPL', num2str(spls(i))]).raster).ts)), 6, '.','k')
                for k = 1:length(psth_summary.(['SPL', num2str(spls(i))]).raster(j).ts)
                    plot([psth_summary.(['SPL', num2str(spls(i))]).raster(j).ts(k), psth_summary.(['SPL', num2str(spls(i))]).raster(j).ts(k)], [j-1,j] +rep *(i-1), 'k', 'LineWidth',1)
                end
            end
        end
    end
    
    % imagesc(spike.scmatrix),colormap(BW)
    xlim([window(1),window(end)])
    ylim([0,n_spls * rep])
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