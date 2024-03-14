function psth1 = paired_pulse_preprocess(spikedata, neuron_num, fig)
innerVar = spikedata(1).clusterData.stimData.inner_variables;
innerIndx = spikedata(1).clusterData.stimData.innerIndexes;
num_inner = length(innerVar);

% neuron_num = 1;

for i = 1:max(innerIndx) % each inner index represents a different interval
    indx = find(innerIndx==i);
    psth_idx = indx;
    psth = spikedata(neuron_num).clusterData.psth;
    psth1.(['Interval',num2str(i)]).raster = psth.raster(psth_idx );
    psth1.(['Interval',num2str(i)]).scmatrix = psth.scmatrix(psth_idx ,:);
    psth1.(['Interval',num2str(i)]).stimulus = psth.stimulus;
end

if fig ==1
    iter = 0;
%     figure
    CT=cbrewer('div', 'RdYlBu', 6); % for nice color
    colorIdx = [2,6];
    for i = 1:num_inner
        spike = psth1.(['Interval',num2str(i)]);
        delay = spike.stimulus.delay(1);
        width = spike.stimulus.width(1);
        num_trial = length(spike.raster);
        hold on
        rectangle('Position',[delay,0 + iter, width,num_trial],'FaceColor', [0.5,0.5,0.5], 'EdgeColor', [0.5,0.5,0.5])
        rectangle('Position',[delay + innerVar(i)+ width, 0 + iter, width, num_trial],'FaceColor', [0.5,0.5,0.5], 'EdgeColor', [0.5,0.5,0.5])
        for k = 1:length(spike.raster)
            if ~isempty(spike.raster(k).ts)
                %                 scatter(spike.raster(k).ts, k*ones(size(spike.raster(k).ts)), 6, '.','k')
                for j = 1:length(spike.raster(k).ts)
                    if mod(i,2) ==1
                        plot([spike.raster(k).ts(j), spike.raster(k).ts(j)], [k-1,k] + iter, 'Color', [0.7, 0.1, 0.16], 'LineWidth',1.5)
                        %                     plot([spike.raster(k).ts(j), spike.raster(k).ts(j)], [k-1,k] + iter, 'Color', 'r', 'LineWidth',1)
                        
                    else
                        plot([spike.raster(k).ts(j), spike.raster(k).ts(j)], [k-1,k] + iter, 'Color', [0.13, 0.4, 0.67], 'LineWidth',1.5)
                        %                     plot([spike.raster(k).ts(j), spike.raster(k).ts(j)], [k-1,k] + iter, 'Color', 'b', 'LineWidth',1)
                        
                    end
                end
            end
        end
        iter = iter + num_trial;
    end
    xlim([1, 1000])
    yticks(12.5 + [25:50:475])
    labels = cellstr(num2str(innerVar'));
    yticklabels(labels(2:2:end))
    ylabel('Intervals')
    xlabel('Time (ms)')
    set(gcf, 'position', [100, 100, 600, 600])
end