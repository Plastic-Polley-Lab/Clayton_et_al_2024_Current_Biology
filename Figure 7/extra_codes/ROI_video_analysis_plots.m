function ROI_video_analysis_plots(t, data, level, CT, region)
genotype = {'wt', 'ko'};
for j = 1:length(genotype)
    figure
    for i = 1:length(level)
        h(i) = lineplot_error(t,data.(genotype{j}).(['level', num2str(level(i))]), 'k')
        set(h(i),'Color', CT(i,:), 'LineWidth', 1)
        hold on
    end
    xlim([-0.2, 1])
    xlabel('Time(s)')
    ylabel(['\Delta Normalized ', region ' Movement'])
    legend(h, {'35 dB', '45 dB', '55 dB', '65 dB','75 dB', '85 dB','95 dB'})
    set(gca,'TickDir','out')
    set(gca,'fontsize',12)
    set(gca,'TickLengt', [0.015 0.015]);
    set(gca, 'LineWidth',1)
    set(gcf,'position',[100,200,400,400])
    title(['PTCHD1 ', upper(genotype{j})])
    set(gcf, 'Color', 'w')
    ylim([-0.5, 2.5])
%     export_fig([upper(genotype{j}),'_',region, 'Change_Tosca'],  '-png', '-pdf')
end