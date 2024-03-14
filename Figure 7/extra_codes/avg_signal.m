function data = avg_signal(data)
for z = 1:length(data)
    dFibPhot = data(z).dFibPhot;
    type = {'REG','RAND'};
    fiberTab={'fiber1','fiber2'};
    for i = 1:length(dFibPhot.trial)
        for j = 1:length(type)
            dFibPhot.trial(i).(type{j}).avg_fiber1 = mean(dFibPhot.trial(i).(type{j}).fiber1);
            dFibPhot.trial(i).(type{j}).avg_fiber2 = mean(dFibPhot.trial(i).(type{j}).fiber2);
        end
    end
    
    avg_signal.cyc = sort(unique([dFibPhot.trial.cycTab]));
    avg_signal.sets = sort(unique([dFibPhot.trial.sets]));
    
    for j = 1:length(fiberTab)
        for k = 1:length(type)
            data(z).(type{k}).avg_signal.(fiberTab{j}) = zeros(length(avg_signal.cyc), length(avg_signal.sets));
            %     figure;
            for i = 1:length(dFibPhot.trial)
                cyc = dFibPhot.trial(i).cycTab;
                cyc_idx = find(avg_signal.cyc == cyc);
                sets = dFibPhot.trial(i).sets;
                data(z).(type{k}).avg_signal.(fiberTab{j})(cyc_idx, sets) = dFibPhot.trial(i).(type{k}).(['avg_',fiberTab{j}]);
                %         hold on
                %         scatter(cyc, dFibPhot.trial(i).REG.(['avg_',fiberTab{j}]), 'o', 'filled', 'MarkerFaceColor',color(sets,:))
                %         xlim([2, 14])
            end
            data(z).(type{k}).avg_signal.cyc = avg_signal.cyc;
            data(z).(type{k}).avg_signal.sets = avg_signal.sets;
        end
        %     title(fiberTab{j})
    end
end
