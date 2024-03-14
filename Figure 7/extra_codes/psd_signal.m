function data = psd_signal(data)
    dFibPhot = data.dFibPhot;
    type = {'REG','RAND'};
    fiberTab={'fiber1','fiber2'};
    avg_signal.cyc  = sort(unique([dFibPhot.trial.cycTab]));
    avg_signal.sets = sort(unique([dFibPhot.trial.sets]));
    
    for j = 1:length(fiberTab)
        for k = 1:length(type)
            data.(type{k}).psd_signal.(fiberTab{j}) = cell(length(avg_signal.cyc), length(avg_signal.sets));
            data.(type{k}).psd_signal.f          = cell(length(avg_signal.cyc), length(avg_signal.sets));
            %     figure;
            for i = 1:length(dFibPhot.trial)
                cyc = dFibPhot.trial(i).cycTab;
                cyc_idx = find(avg_signal.cyc == cyc);
                sets = dFibPhot.trial(i).sets;
                x = dFibPhot.trial(i).(type{k}).(fiberTab{j});
                fs = dFibPhot.fs;
                [pxx, f] = periodogram(x,[], length(x),fs);
                data.(type{k}).psd_signal.(fiberTab{j}){cyc_idx, sets} = pxx;
                data.(type{k}).psd_signal.f{cyc_idx, sets} = f;
%         hold on
                %         scatter(cyc, dFibPhot.trial(i).REG.(['avg_',fiberTab{j}]), 'o', 'filled', 'MarkerFaceColor',color(sets,:))
                %         xlim([2, 14])
            end
            data.(type{k}).psd_signal.cyc = avg_signal.cyc;
            data.(type{k}).psd_signal.sets = avg_signal.sets;
        end
        %     title(fiberTab{j})
    end
    
    % plot the psd
    for k = 1:2
        figure;
        color =cbrewer('div', 'RdYlBu',4);
        for i = 1:5
            subplot(5, 1, i)
            for j = 1:4
                f = data.REG.psd_signal.f{i, j};
                pxx = data.REG.psd_signal.(fiberTab{k}){i, j};
                if j ~= 4
                    plot(f, pxx, 'Color', color(j,:))
                else
                    plot(f, pxx, 'Color', 'k')
                end
                xlim([0.4,6])
                hold on
            end
            legend({'Set1','Set2', 'Set3','RAND-Control'})
            title([data.animal, '', fiberTab{k}, ' Cycle Size', ' ', num2str(avg_signal.cyc(i))])
        end
        set(gcf, 'Position', [200, 200, 600, 800] )    
    end

