function pattern_result = plot_pattern_raster(neuron_num, spikedata, fig)
% INPUT:
%      fig: 1 is to plot the figure.
spike = spikedata(neuron_num).clusterData.psth;
% for j = 1:10 % there are 10 different patterns
%     idx_start = 22 + (j-1) * 100; % 100 interval for each pattern
%     idx_stop  = 21 + j * 100;
%     psth_summary.(['pattern', num2str(j)]).raster = spike.raster(idx_start:idx_stop);
%     psth_summary.(['pattern', num2str(j)]).scmatrix = spike.scmatrix(idx_start:idx_stop, :);
%     psth_summary.(['pattern', num2str(j)]).stimulus = spike.stimulus;
% end
for j = 1:2 % there are 2 different patterns
    idx_start = 52 + (j-1) * 100; % 100 interval for each pattern
    idx_stop  = 51 + j * 100;
    psth_summary.(['pattern', num2str(j)]).raster = spike.raster(idx_start:idx_stop);
    psth_summary.(['pattern', num2str(j)]).scmatrix = spike.scmatrix(idx_start:idx_stop, :);
    psth_summary.(['pattern', num2str(j)]).stimulus = spike.stimulus;
end
psth_avg = mean(spike.scmatrix,1)/0.001;
psth_avg_smooth = smoothts(psth_avg,'g',5,1);


window = 1:100;
if fig ==1
    figure;
    subaxis(2,2,1 , 'sh', 0.1, 'sv', 0.1)
    CT=cbrewer('seq', 'PuRd', 2); % for nice color
%     CT=cbrewer
    load('E:\Ke_Chen\MATLAB\Ephys\EphysAnalysis\BW.mat')
    hold on
    for i = 1:2
%         idx_start = 22 + (i-1) * 100; % 100 interval for each pattern
        idx_start = 52 + (i-1) * 100; % 100 interval for each pattern

        hr(i) = rectangle('Position',[spike.stimulus.delay(1),idx_start, spike.stimulus.width(1),100],'FaceColor', CT(i,:), 'EdgeColor', CT(i,:))
    end
    % imagesc(spike.scmatrix),colormap(BW)
    hold on
    for i = 1:length(spike.raster)
        if ~isempty(spike.raster(i).ts)
            scatter(spike.raster(i).ts, i*ones(size(spike.raster(i).ts)), 6, '.','k')
            for j = 1:length(spike.raster(i).ts)
                plot([spike.raster(i).ts(j), spike.raster(i).ts(j)], [i-1,i], 'k', 'LineWidth',1)
            end
        end
    end
    hold off
    xlim([window(1),window(end)])
    ylim([0,length(spike.raster)])
    axis on
    ylabel('Trial #')
    xlabel('Time (ms)')
%     legend([hr(1), hr(2)], {'reg', 'rand'})
end

for i = 1: 2
    firingRate(i,:) =smoothts(mean(psth_summary.(['pattern', num2str(i)]).scmatrix,1)/0.001, 'g', 5, 1);
end

if fig ==1
    subaxis(2, 2, 3, 'sh', 0.1, 'sv', 0.1)
    for i = 1: 2
        plot(window -0.5, firingRate(i,:), 'color', CT(i,:))
        hold on
    end
    hold off
    xlim([window(1),window(end)])
    axis on
    ylabel('Firing Rate (Hz)')
    xlabel('Time (ms)')
end

clear normalizedACF
% for i = 1:2
%     [normalizedACF(i,:), lags] = autocorr(firingRate(i,:),99);
% end

for i = 1:2
    [c,lags_temp] = xcorr(firingRate(i,:),'coeff');
    lags(i,:) = lags_temp((length(lags_temp)+1)/2:end);
    normalizedACF(i,:) = c((length(lags_temp)+1)/2:end);
end

if fig ==1
    subaxis(2,2,2,'sh', 0.1, 'sv', 0.1)
    for i = 1:2
        h(i) = plot(lags(i,:), normalizedACF(i,:), 'color', CT(i,:))
        hold on
    end
    plot([0,100], [0,0], '--', 'color', [0.5, 0.5, 0.5])
    xlim([window(1),window(end)])
    axis on
    ylim([-0.2,1])
    ylabel('ACF')
    xlabel('Time (ms)')
    legend([h(1), h(2)], {'reg', 'rand'})
end


% for i = 1:2  
%     idx = find(normalizedACF(i,:)< 0);
%     lags_fit{i} = lags(1:idx(1));
%     normalizedACF_fit{i} = normalizedACF(i,1:idx(1));
%     % title(['Jitter ascending; Unit', num2str(neuron_num)])
%     model(i).f = fit(lags_fit{i}',normalizedACF_fit{i}','exp1');
%     model(i).tau = -1/model(i).f.b;
% end

for i = 1:2
    options = fitoptions('exp1');
    options.StartPoint = [1 -0.5];
    % options.Upper = [Inf 0];
    options.Upper = [1 0];
    options.Lower = [1 -Inf];
    [curve,gof] = fit(lags(i,:)',normalizedACF(i,:)','exp1',options);
    model(i).tau = -1/(curve.b);
    model(i).tau_gof_single = gof.adjrsquare;
    model(i).exp_fit_type = 1;
    if (gof.adjrsquare<0.75)
        model(i).exp_fit_type = 2;
        options = fitoptions('exp2');
        options.StartPoint = [1 -0.5 0.5 -0.3];
        options.Upper = [Inf 0 Inf -0.01];
        %     options.Lower = [0 -Inf 0 -Inf];
        [curve,gof] = fit(lags(i,:)',normalizedACF(i,:)','exp2',options);
        tau1 = -1/(curve.b);
        tau2 = -1/(curve.d);
        model(i).tau_final = (curve.a*tau1+curve.c*tau2)/(curve.a+curve.c);
        if model(i).tau_final>100
            model(i).tau_final = model(i).tau;
            model(i).exp_fit_type = 1;
        end
    else
        model(i).tau_final = model(i).tau;
    end
    model(i).tau_gof_single = gof.adjrsquare; 
    model(i).f = curve;
end




if fig ==1
    subaxis(2,2,4,'sh', 0.1, 'sv', 0.1)
    hold off
    for i = 1:2
        h(i) = scatter(lags(i,:), normalizedACF(i,:), [], CT(i,:), 'filled');
        hold on
        h1= plot(model(i).f);
        set(h1, 'Color', CT(i,:))
    end
%     legend([h(1), h(2), h(3), h(4), h(5), h(6), h(7), h(8), h(9), h(10)], ...
%         {sprintf('tau = %2.2f ms', model(1).tau), sprintf('tau = %2.2f ms', model(2).tau), ...
%         sprintf('tau = %2.2f ms', model(3).tau), sprintf('tau = %2.2f ms', model(4).tau),...
%         sprintf('tau = %2.2f ms', model(5).tau), sprintf('tau = %2.2f ms', model(6).tau), ...
%         sprintf('tau = %2.2f ms', model(7).tau), sprintf('tau = %2.2f ms', model(8).tau), ...
%         sprintf('tau = %2.2f ms', model(9).tau), sprintf('tau = %2.2f ms', model(10).tau)})
    legend([h(1), h(2)], ...
        {sprintf('reg tau = %2.2f ms', model(1).tau_final), sprintf('rand tau = %2.2f ms', model(2).tau_final)})   
    xlabel('Time (s)')
    ylabel('ACF')
    set(gcf,'position',[100,200,800,800])
    suptitle(['Neuron ', num2str(neuron_num)])
    
    if neuron_num < 10
        print(['Jitter-ascending_response_0', num2str(neuron_num)],'-dpdf','-bestfit')
%            set(gcf, 'Color', 'w')
%            export_fig(['Jitter-ascending_response_0', num2str(neuron_num)],  '-png')
    else
        print(['Jitter-ascending_response_', num2str(neuron_num)],'-dpdf','-bestfit')
%            set(gcf, 'Color', 'w')
%            export_fig(['Jitter-ascending_response_', num2str(neuron_num)],  '-png')
    end
    close
end
% save the result
pattern_result.psth_summary = psth_summary;
pattern_result.firingRate   = firingRate;
pattern_result.normalizedACF= normalizedACF;
pattern_result.lags = lags;
% pattern_result.lags_fit = lags_fit;
% pattern_result.normalizedACF_fit = normalizedACF_fit;
pattern_result.model = model;
end