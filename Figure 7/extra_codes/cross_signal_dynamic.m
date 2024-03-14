function data = cross_signal_dynamic(data)
dFibPhot = data.dFibPhot;
type = {'RAND1','REG','RAND'};
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
            stimulus = dFibPhot.trial(i).(type{k}).stimuli;
            % add by ke to make the stimuli have the same number
            %                 stimulus = stimulus(1:101);
            Fs = dFibPhot.fs;
            signal = writeaudio_signal(stimulus, Fs);
            [t_on, ~]=Timing_onset_offset(signal, 1:length(signal), 0,30,0);
            t_on = t_on-1;
            
            segments_on = t_on(1:cyc *5:end); % change the size of frame here
            for zz = 1:length(segments_on)
                if zz< length(segments_on)
                    [acor,lags] = xcorr(dF(segments_on(zz):segments_on(zz+1)-1),signal(segments_on(zz):segments_on(zz+1)-1),'coeff', 800);
                else
                    [acor,lags] = xcorr(dF(segments_on(zz):length(signal)),signal(segments_on(zz):length(signal)),'coeff', 800);
                    
                end
                data.(type{k}).cross_signal.(fiberTab{j}){cyc_idx, sets}(zz,:) = acor;
                data.(type{k}).cross_signal.lags{cyc_idx, sets}(zz,:) = lags/Fs;
            end

            %         hold on
            %         scatter(cyc, dFibPhot.trial(i).REG.(['avg_',fiberTab{j}]), 'o', 'filled', 'MarkerFaceColor',color(sets,:))
            %         xlim([2, 14])
        end
        data.(type{k}).cross_signal.cyc = avg_signal.cyc;
        data.(type{k}).cross_signal.sets = avg_signal.sets;
    end
end
    %     title(fiberTab{j})
    
    % plot the autocorrelation
    %     for k = 1:2
    %         figure;
    %         color =cbrewer('div', 'RdYlBu',4);
    %         for i = 1:5
    %             for z = 1:2
    %                 subplot(5, 2, (i-1)*2 + z)
    %                 for j = 1:4
    %                     lags = data.(type{z}).cross_signal.lags{i, j};
    %                     acor = data.(type{z}).cross_signal.(fiberTab{k}){i, j};
    %                     if j ~= 4
    %                         plot(lags, acor, 'Color', color(j,:))
    %                     else
    %                         plot(lags, acor, 'Color', 'k')
    %                     end
    %                     xlim([-0.8,0.8])
    %                     ylim([-0.3,0.3])
    %                     hold on
    %                 end
    %                 legend({'Set1','Set2', 'Set3','RAND-Control'})
    %                 ylabel([' Cycle Size', ' ', num2str(avg_signal.cyc(i))])
    %                 xlabel('Lag (s)')
    %             end
    %         end
    %         suptitle(['REG      ', data.animal, ' ', location{k}, 'Cross-Correlation', '     RAND'])
    %         set(gcf, 'Position', [200, 200, 600, 800] )
    %     end
    
