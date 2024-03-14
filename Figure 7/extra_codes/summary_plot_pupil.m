function plots = summary_plot_pupil(summaryData, additional_fig)
%% plot the example normalized trace of pupil fluctuation
example_id_wt = 9;
example_id_ko = 8;
wt_time = summaryData.wt(example_id_wt).time;
wt_pupil_n = summaryData.wt(example_id_wt).pupil_area_n;
wt_event_t = summaryData.wt(example_id_wt).event_time;


ko_time = summaryData.ko(example_id_ko).time;
ko_pupil_n = summaryData.ko(example_id_ko).pupil_area_n;
ko_event_t = summaryData.ko(example_id_ko).event_time;

figure;
h1 = plot(wt_time, wt_pupil_n, '-k');
hold on
for i = 1:length(wt_event_t)
    plot([wt_event_t(i), wt_event_t(i)], [0, 1], '-k')
end


h2 = plot(summaryData.ko(example_id_ko).time, summaryData.ko(example_id_ko).pupil_area_n, '-r');
hold on
for i = 1:length(ko_event_t)
    plot([ko_event_t(i), ko_event_t(i)], [0, 1], '-r')
end
xlim([120,180])
xlabel('Time (s)')
ylabel('Norma. pupil size')
legend([h1, h2], {'WT', 'KO'})
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,200])
set(gcf, 'Color', 'w')

%% plot the spontaneous pupil fluctuation
% 1st just plot one sessio

figure;
title_label = 'Group Data';
plots.wt_spont = [summaryData.wt.sponta_n_avg];
plots.ko_spont = [summaryData.ko.sponta_n_avg];
[fig1, fig2] = ecdf_bar_plot(plots.wt_spont, plots.ko_spont, title_label);
set(get(fig1,'XLabel'), 'String', 'Spontaneous Pupil Size');
set(get(fig2,'YLabel'), 'String', 'Spontaneous Pupil Size');
%export_fig('Spontaneous Pupil Size',  '-png')

%% plot the sound evoke pupil size (Not normalized to the baseline)
genotype = {'wt', 'ko'};
for i = 1 : length(genotype)
    temp = summaryData.(genotype{i});
    for j = 1:length(temp)
        temp_data = temp(j).evoke_baseline;
        
        for k = 1 : length(temp_data)
            summaryData.(genotype{i})(j).evoke_avg_raw(k, :) = mean(temp_data{k});
            
        end
    end
end
pupil_evoked_trace_plots(summaryData, 'evoke_avg_raw')
end
