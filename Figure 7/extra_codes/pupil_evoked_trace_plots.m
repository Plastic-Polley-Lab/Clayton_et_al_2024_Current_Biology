function pupil_evoked_trace_plots(summaryData, fieldName)

% data = summary
% fieldName = 'evoke_avg';
plots.n_level = size(summaryData.wt(1).(fieldName), 1);
plots.n_frame = size(summaryData.wt(1).(fieldName), 2);
plots.wt_evoked = reshape([summaryData.wt.(fieldName)], plots.n_level, plots.n_frame, []) ;
plots.ko_evoked = reshape([summaryData.ko.(fieldName)], plots.n_level, plots.n_frame, []) ;

figure;
plots.t = summaryData.wt(1).t + 0.033; % fix the issue that event occurs between 60 and 61th frame
CT=cbrewer('seq', 'PuRd', 7);
for i = 1:plots.n_level
    %     h(i) = lineplot_error(t,squeeze(wt_evoked(i,:,:))', 'k');
    h(i) = line_sem_plot(plots.t, squeeze(plots.wt_evoked(i,:,:))', 'k')
    set(h(i),'Color', CT(i,:), 'LineWidth', 1)
    hold on
end
xlim([-1, 3])
xlabel('Time(s)')
ylabel('\Delta Normalize Pupil Size')
legend(h, {'35 dB', '45 dB', '55 dB', '65 dB','75 dB', '85 dB','95 dB'})
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
title('PTCHD1 WT')
set(gcf, 'Color', 'w')

figure;
CT=cbrewer('seq', 'PuRd', 7);
for i = 1:plots.n_level
    %     h(i) = lineplot_error(t,squeeze(wt_evoked(i,:,:))', 'k');
    h(i) = line_sem_plot(plots.t, squeeze(plots.ko_evoked(i,:,:))', 'k')
    set(h(i),'Color', CT(i,:), 'LineWidth', 1)
    hold on
end
xlim([-1, 3])
xlabel('Time(s)')
ylabel('\Delta Normalize Pupil Size')
legend(h, {'35 dB', '45 dB', '55 dB', '65 dB','75 dB', '85 dB','95 dB'})
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
title('PTCHD1 KO')
set(gcf, 'Color', 'w')


%% Let's also get the Maximal pupil change
t = plots.t;
indx = find(t>0 & t<2); % use the interval between [0, 2] seconds
for i = 1:plots.n_level
    temp = squeeze(plots.wt_evoked(i,:,:))';
    plots.wt_Max_pupil(:,i) = max(temp(:,indx),[],2);
    %     ko.pupi_area(:,i) = sum(temp(:,indx),2)*step;
end

for i = 1:plots.n_level
    temp = squeeze(plots.ko_evoked(i,:,:))';
    plots.ko_Max_pupil(:,i) = max(temp(:,indx),[],2);
    %     wt.pupi_area(:,i) = sum(temp(:,indx),2)*step;
    
end
figure
line_errorbar_drc(plots.wt_Max_pupil, plots.ko_Max_pupil)
% line_errorbar_drc(wt.pupi_area./max(wt.pupi_area,[],2), ko.pupi_area./max(ko.pupi_area,[],2))
% line_errorbar_drc(wt.pupi_area./wt.pupi_area(:,end), ko.pupi_area./ko.pupi_area(:,end))

xticks([0, 1:7, 8])
xlim([0,8])
xticklabels({'','35', '45', '55', '65','75', '85','95',''})
ylabel(['\Delta Normalized Pupil Size'])
xlabel('Sound Intensity (dB)')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
% title('Noise induced Whisker Movement')
% title(['Noise induced ', label, ' Movement'])
% plot(plots.wt_Max_pupil', 'color', [0.5, 0.5, 0.5])
% plot(plots.ko_Max_pupil', 'color', [0.9766, 0.6211, 0.7070])
set(gcf, 'Color', 'w')