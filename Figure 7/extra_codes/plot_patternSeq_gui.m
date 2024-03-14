function plot_patternSeq_gui(t_jitter, CycleLength, jitters, jitter_type)
figure

color =cbrewer('seq', 'YlGnBu', 9); % for nice color;
% color = lines(length(C));
% gscatter(1:length(t),t,group',color)
if strcmp(jitter_type, 'jitter-swap')
    for k = 1:CycleLength
        subplot(size(t_jitter,1),1,k)
        for i = 1:CycleLength
            %             scatter(i:CycleLength:size(t_jitter,2), t_jitter(k,i:CycleLength:end),'o','filled','MarkerEdgeColor', color(i,:), 'MarkerFaceColor', color(i,:))
            %             ylim([150,550])
            %             ylim([100,400])
            %             ylabel('Interval (ms)')
            scatter(i:CycleLength:size(t_jitter,2),i* ones(size(t_jitter(k,i:CycleLength:end))),40, t_jitter(k,i:CycleLength:end),'o','filled','MarkerEdgeColor', [0,0,0])
            colormap(color)
            cb = colorbar;
            caxis([150, 350])
            cb.Ticks = [150, 250, 350] ; %Create 8 ticks from zero to 1
            cb.TickLabels = {'150', '250', '350'} ;
            cb.Label.String = 'Interval (ms)';
            ylim([0, CycleLength+1])
            hold on
        end
        if k ==1
            title('REG')
        elseif k < CycleLength
            title(['FIX', ' ', num2str(CycleLength-k), ' Interval'])
        else
            title('RAND')
        end
    end
    xlabel('Sound sequence #')
    set(gcf,'position',[100,50,800,600])
else
    for k = 1:size(t_jitter,1)
        %         subplot(size(t_jitter,1),1,k)
        %         subaxis(size(t_jitter,1),1,k, 'sh', 0.01, 'sv', 0.01, 'padding', 0, 'margin', 0.1);
        subaxis(size(t_jitter,1),1,k, 'sv', 0.01);
        
        for i = 1:CycleLength
            scatter(i:CycleLength:size(t_jitter,2), t_jitter(k,i:CycleLength:end),'o','filled','MarkerEdgeColor', color(i,:), 'MarkerFaceColor', color(i,:))
            hold on
        end
        ylim([150,550])
        if k == ceil(size(t_jitter,1)/2)
            ylabel('Interval (ms)', 'FontWeight','bold','Color','k')
            set(gca,'xtick',[])
            yyaxis right
            set(gca,'ytick',[])
            ylabel([jitter_type, ' ', num2str(jitters(k)*100), '% '],'FontWeight','bold','Color','k')
        elseif k == size(t_jitter,1)
            yyaxis right
            set(gca,'ytick',[])
            ylabel([num2str(jitters(k)*100), '% '], 'FontWeight','bold','Color','k')
            
        else
            set(gca,'xtick',[])
            yyaxis right
            set(gca,'ytick',[])
            ylabel([num2str(jitters(k)*100), '% '], 'FontWeight','bold','Color','k')
        end
        %         title(['Jitter ',num2str(jitters(k)*100), '% ', jitter_type])
    end
    xlabel('Sound sequence #', 'FontWeight','bold','Color','k')
    set(gcf,'position',[100,50,600,1000])
end