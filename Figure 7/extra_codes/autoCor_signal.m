function data = autoCor_signal(data)
    dFibPhot = data.dFibPhot;
    type = {'REG','RAND'};
    fiberTab={'fiber1','fiber2'};
    location = {'rostral', 'caudal'};
    avg_signal.cyc  = sort(unique([dFibPhot.trial.cycTab]));
    avg_signal.sets = sort(unique([dFibPhot.trial.sets]));
    
    for j = 1:length(fiberTab)
        for k = 1:length(type)
            data.(type{k}).autoCor_signal.(fiberTab{j}) = cell(length(avg_signal.cyc), length(avg_signal.sets));
            data.(type{k}).autoCor_signal.lags          = cell(length(avg_signal.cyc), length(avg_signal.sets));
            %     figure;
            for i = 1:length(dFibPhot.trial)
                cyc = dFibPhot.trial(i).cycTab;
                cyc_idx = find(avg_signal.cyc == cyc);
                sets = dFibPhot.trial(i).sets;
                dF = dFibPhot.trial(i).(type{k}).(fiberTab{j});
                fs = dFibPhot.fs;
                [normalizedACF, lags] = autocorr(dF,5000);
%                 unnormalizedACF = normalizedACF*var(dF,1);
                data.(type{k}).autoCor_signal.(fiberTab{j}){cyc_idx, sets} = normalizedACF;
                data.(type{k}).autoCor_signal.lags{cyc_idx, sets} = lags/fs;
%         hold on
                %         scatter(cyc, dFibPhot.trial(i).REG.(['avg_',fiberTab{j}]), 'o', 'filled', 'MarkerFaceColor',color(sets,:))
                %         xlim([2, 14])
            end
            data.(type{k}).autoCor_signal.cyc = avg_signal.cyc;
            data.(type{k}).autoCor_signal.sets = avg_signal.sets;
        end
        %     title(fiberTab{j})
    end
    
    % plot the autocorrelation
    for k = 1:2
        figure;
        color =cbrewer('div', 'RdYlBu',4);
        for i = 1:5
            subplot(5, 1, i)
            for j = 1:4
                lags = data.REG.autoCor_signal.lags{i, j};
                normalizedACF = data.REG.autoCor_signal.(fiberTab{k}){i, j};
                if j ~= 4
                    plot(lags, normalizedACF, 'Color', color(j,:))
                else
                    plot(lags, normalizedACF, 'Color', 'k')
                end
                xlim([0,5])
                ylim([-0.1,1])
                hold on
            end
            legend({'Set1','Set2', 'Set3','RAND-Control'})
            ylabel([' Cycle Size', ' ', num2str(avg_signal.cyc(i))])
            xlabel('Lag (s)')
        end
        suptitle([data.animal, ' ', location{k}, 'AutoCorrelation'])
        set(gcf, 'Position', [200, 200, 600, 800] )    
    end

