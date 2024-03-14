function plots = summary_plot_orofacial(summaryData, additional_fig, cutOff, repeats)
% 1st, let's dealing with the repeats;
if ~isempty(repeats)
    summaryData_raw = summaryData;
    genotype = {'wt', 'ko'};
    for i = 1: 2
        label = genotype{i};
        repeats_indx = repeats.(label);
        remove_row = sort(cell2mat(repeats_indx));
        
        data_label = summaryData.(label);
        data_label(remove_row) =[];
        
        for j = 1:length(repeats_indx);
            temp_data = summaryData.(label)(repeats_indx{j});
            avg_data(j).filename = [temp_data.filename];
            avg_data(j).genotype = temp_data(1).genotype;
            temp_baseline = [];
            for k = 1:length(temp_data)
                avg_data(j).Oro_post{k} = temp_data(k).Oro_post;
%                 avg_data(j).Oro_ant{k} = temp_data(k).Oro_ant;

                avg_data(j).video_summary(k) = temp_data(k).video_summary;
                avg_data(j).sponta{k}         = temp_data(k).sponta;
                avg_data(j).sponta_n{k}       = temp_data(k).sponta_n;
                avg_data(j).evoke_baseline(k,:) = temp_data(k).evoke_baseline;
                avg_data(j).evoke_n(k,:) = temp_data(k).evoke_n;
                avg_data(j).evoke_avg(:,:,k) = temp_data(k).evoke_avg;
                avg_data(j).t = temp_data(k).t;
                avg_data(j).orofacial_n{k} = temp_data(k).orofacial_n;
                avg_data(j).time{k}        = temp_data(k).time;
                avg_data(j).event_time{k}  = temp_data(k).event_time;
                avg_data(j).level          = temp_data(k).level;
                avg_data(j).spont_move{k}  = temp_data(k).spont_move;
                
                
            end
            avg_data(j).sponta_n_avg = mean([temp_data.sponta_n_avg]);
            avg_data(j).spont_move_avg = mean(reshape([temp_data.spont_move_avg],2, []),2);
            avg_data(j).evoke_baseline = cellfun(@(a, b) cat(1, a, b),...
                avg_data(j).evoke_baseline(1,:), avg_data(j).evoke_baseline(2,:), 'UniformOutput',false);
            
            avg_data(j).evoke_n = cellfun(@(a, b) cat(1, a, b),...
                avg_data(j).evoke_n(1,:), avg_data(j).evoke_n(2,:), 'UniformOutput',false);
            
            avg_data(j).evoke_avg = squeeze(mean(avg_data(j).evoke_avg, 3));
        end
        
        data_label = [data_label, avg_data];
        
        summaryData.(label) = data_label;
    end
end
%%
%% plot the spontaneous orofacial movements
example_id_wt = 10;
example_id_ko = 6;
wt_time = summaryData.wt(example_id_wt).time;
wt_pupil_n = summaryData.wt(example_id_wt).orofacial_n;
wt_event_t = summaryData.wt(example_id_wt).event_time;
wt_deviation = std(wt_pupil_n);


ko_time = summaryData.ko(example_id_ko).time;
ko_pupil_n = summaryData.ko(example_id_ko).orofacial_n;
ko_event_t = summaryData.ko(example_id_ko).event_time;
ko_deviation = std(ko_pupil_n);

figure;
h1 = plot(wt_time, wt_pupil_n, '-k');
hold on
for i = 1:length(wt_event_t)
    plot([wt_event_t(i), wt_event_t(i)], [0, 10], '-k')
end
plot(wt_time, cutOff* wt_deviation * ones(size(wt_time)), '--k')

h2 = plot(ko_time, ko_pupil_n, '-r');
hold on
for i = 1:length(ko_event_t)
    plot([ko_event_t(i), ko_event_t(i)], [0, 10], '-r')
end
plot(ko_time, cutOff * ko_deviation * ones(size(ko_time)), '--r')

xlim([120,180])
xlabel('Time (s)')
ylabel('Norm. Orofacial')
legend([h1, h2], {'WT', 'KO'})
box off
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,200])
set(gcf, 'Color', 'w')


%% plot the sound evoked orofacial movements

plots.n_level = size(summaryData.wt(1).evoke_avg, 1);
plots.n_frame = size(summaryData.wt(1).evoke_avg, 2);
plots.wt_evoked = reshape([summaryData.wt.evoke_avg], plots.n_level, plots.n_frame, []) ;
plots.ko_evoked = reshape([summaryData.ko.evoke_avg], plots.n_level, plots.n_frame, []) ;

figure;
plots.t = summaryData.wt(1).t + 0.033; % 0.033 was added because the event occurs between 60 and 61 frames()
CT=cbrewer('seq', 'PuRd', 7);
for i = 1:plots.n_level
    %     h(i) = lineplot_error(t,squeeze(wt_evoked(i,:,:))', 'k');
    h(i) = line_sem_plot(plots.t, squeeze(plots.wt_evoked(i,:,:))', 'k');
    set(h(i),'Color', CT(i,:), 'LineWidth', 1)
    hold on
end
xlim([-0.5, 1.5])
xlabel('Time(s)')
ylabel('\Delta Orofacial Movements')
legend(h, {'35 dB', '45 dB', '55 dB', '65 dB','75 dB', '85 dB','95 dB'})
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
title('PTCHD1 WT')
set(gcf, 'Color', 'w')
% export_fig('WT_Orofacial_Movements',  '-png', '-pdf')


figure;
CT=cbrewer('seq', 'PuRd', 7);
for i = 1:plots.n_level
    %     h(i) = lineplot_error(t,squeeze(wt_evoked(i,:,:))', 'k');
    h(i) = line_sem_plot(plots.t, squeeze(plots.ko_evoked(i,:,:))', 'k');
    set(h(i),'Color', CT(i,:), 'LineWidth', 1)
    hold on
end
xlim([-0.5, 1.5])
xlabel('Time(s)')
ylabel('\Delta Orofacial Movements')
legend(h, {'35 dB', '45 dB', '55 dB', '65 dB','75 dB', '85 dB','95 dB'})
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
title('PTCHD1 KO')
set(gcf, 'Color', 'w')
% export_fig('KO_Orofacial_Movements',  '-png', '-pdf')

%% Let's also get the Maximal orofacial Movement
t = plots.t;
indx = find(t>0 & t< 0.2); % use the interval between [0, 0.2] seconds
for i = 1:plots.n_level
    temp = squeeze(plots.wt_evoked(i,:,:))';
    plots.wt_Max_Oro(:,i)  = max(temp(:,indx),[],2);
    plots.wt_Oro_area(:,i) = sum(temp(:,indx),2);
end

for i = 1:plots.n_level
    temp = squeeze(plots.ko_evoked(i,:,:))';
    plots.ko_Max_Oro(:,i) = max(temp(:,indx),[],2);
    plots.ko_Oro_area(:,i) = sum(temp(:,indx),2);
    
end
figure
line_errorbar_drc(plots.wt_Max_Oro, plots.ko_Max_Oro)


xticks([0, 1:7, 8])
xlim([0,8])
xticklabels({'','35', '45', '55', '65','75', '85','95',''})
ylabel(['\Delta Orofacial Movements'])
xlabel('Sound Intensity (dB)')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,400,400])
set(gcf, 'Color', 'w')
end
