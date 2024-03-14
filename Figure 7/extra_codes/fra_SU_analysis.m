function [fra, spls, freqs, tuning, rlf, spont] = fra_SU_analysis(neuron_num, spikedata, ana_window, fig)
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
freqs = floor(spikedata(neuron_num).clusterData.stimData.outer_variables/1e2)*1e2; %outer seq %just rounding off to octaves recognizable
n_freqs = length(freqs);
innerIndexs = spikedata(neuron_num).clusterData.stimData.innerIndexes;
outerIndexs = spikedata(neuron_num).clusterData.stimData.outerIndexes;
n_rep = spikedata(1).clusterData.stimData.numTrials/ (n_spls * n_freqs);
%%
fra = zeros(n_freqs, n_spls, n_rep);
% figure
k = 1;
for i = 1: n_freqs
    for j = 1:n_spls
        idx = find(outerIndexs == i & innerIndexs == j);
        scmatrix = psth.scmatrix(idx,:);
        %         subplot(n_freqs,n_spls,k)
        %         fr_avg = mean(scmatrix);
        %         psth_avg_smooth = smoothts(fr_avg,'b',10);
        %         plot(psth_avg_smooth/0.001)
        %         ylim([0,500])
        %         xlim([0, 100])
        %         k = k+1
        %         title([num2str(freqs(i)), '-', num2str(spls(j))])
        scmatrix_sum = sum(scmatrix);
        delay = psth.stimulus.delay;
        duration = psth.stimulus.width;
        analysis_window = ana_window + delay;
        fra(i,j,:) = sum(scmatrix(:,analysis_window),2);
    end
end

if fig == 1
    subplot(2,2,1)
    fra_plotSU = sum(fra,3);
    max_fra = max(max(fra_plotSU))
    imagesc(rot90(fra_plotSU/max_fra))
    colormap(hot)
    x_ticks_indx = 1:1:length(freqs);
    xticks(x_ticks_indx)
    labels = {};
    for i = 1:length(x_ticks_indx)
        labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
    end
    xticklabels(labels)
    yticks([1:1:length(spls)])
    for i = 1:length(spls)
        labels_spls{i} = num2str(spls(i));
    end
    yticklabels(flip(labels_spls(1:1:end)))
    xlabel('Frequency (kHz)')
    ylabel('Level (dB SPL)')
    box off
    set(gca,'TickDir','out')
    % set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    set(gcf,'position',[100,600,400,350])
end
%% plot the frequence tuning
tuning = squeeze(sum(fra, 2));
%
if fig ==1
    subplot(2,2,3)
    data = tuning;
    err = std(data,0,2)/sqrt(size(data,2)); % standard error mean
    CT=cbrewer('div', 'RdYlBu', 6); % for nice color
    e = errorbar(mean(data,2),err, '-s','MarkerSize',6,'MarkerFaceColor', CT(6,:), 'Color',CT(6,:), 'LineWidth',1);
    hold on
    x_ticks_indx = 1:1:length(freqs);
    xticks(x_ticks_indx)
    labels = {};
    for i = 1:length(x_ticks_indx)
        labels{i} = num2str(freqs(x_ticks_indx(i))/1000);
    end
    xticklabels(labels)
    ylabel('Spike Counts')
    xlabel('Frequency (kHz)')
    title(['Frequency', ' Tuning'])
    box off
    set(gca,'TickDir','out')
    % set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    set(gcf,'position',[100,200,400,300])
    hold off
    %
end
%% plot the rate level function
%
[~, freq_indx] = max(mean(tuning,2)); % best frequency tuning
rlf = squeeze(fra(freq_indx,:,:))./(ana_window(end)-ana_window(1) + 1) * 1000;

if fig ==1
    subplot(2,2,4)
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
    title(['Frequency', ' ', num2str(freqs(freq_indx)), ' Hz'])
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
if fig ==1
    subplot(2,2,2)
    spike = spikedata(neuron_num).clusterData.psth;
    
    for j = 1:n_spls % there are 2 different patterns
        idx = find(outerIndexs == freq_indx & innerIndexs == j);
        psth_summary.(['SPL', num2str(spls(j))]).raster = spike.raster(idx);
        psth_summary.(['SPL', num2str(spls(j))]).scmatrix = spike.scmatrix(idx, :);
        psth_summary.(['SPL', num2str(spls(j))]).stimulus = spike.stimulus;
    end
    window = plot_window;
    % if fig ==1
    
    CT=cbrewer('seq', 'PuRd', n_freqs); % for nice color
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
    set(gcf,'position',[100,100,800,800])
    %         legend([hr(1), hr(2)], {'reg', 'rand'})
end