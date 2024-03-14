clear
load('video_roi_burst_sumamry.mat')
group_num = 10;
genotype = {'wt', 'ko'};
for i = 1:2 % genotype
    for j = 1:length(video_sum.(genotype{i})) % animal
        oro_movement_posterior = video_sum.(genotype{i})(j).video_summary.Mtrace_Posterior_norm;
        oro_movement_anterior = video_sum.(genotype{i})(j).video_summary.Mtrace_anterior_norm;
        
        
        iter = floor(size(oro_movement_posterior, 1)/group_num);
        for kk = 1:iter
            indx = (1:group_num)+ group_num*(kk-1);
            video_sum.(genotype{i})(j).group_movment_p_m(kk,:) = mean(oro_movement_posterior(indx, :),1);
            video_sum.(genotype{i})(j).group_movment_p_sem(kk,:) = std(oro_movement_posterior(indx, :),0,1)/sqrt(group_num);
            video_sum.(genotype{i})(j).group_movment_a_m(kk,:) = mean(oro_movement_anterior(indx, :),1);
            video_sum.(genotype{i})(j).group_movment_a_sem(kk,:) = std(oro_movement_anterior(indx, :),0,1)/sqrt(group_num);
            
            
        end
    end
end
%%
genotype = {'ko', 'wt'};
select_group = 10; % Only analyze the first 10 groups of data (e.g. from trial 1 to trial 5 * group_num);
for kk = 1:length(genotype)
    for i = 1:length(video_sum.(genotype{kk}))
        summary.(genotype{kk}).Mtrace_anterior_g(:,:,i) = video_sum.(genotype{kk})(i).group_movment_a_m(1:select_group, :);
        summary.(genotype{kk}).Mtrace_posterior_g(:,:,i) = video_sum.(genotype{kk})(i).group_movment_p_m(1:select_group, :);
    end

end
%%
figure;
step = 0.033;
t= -0.033 * 9: step: -0.033; % the event is between 10 and 11 th frame
t_post = 0:step:17*step;
t = [t, t_post];
group_indx = 1

h(1) = lineplot_error(t,(squeeze(summary.wt.Mtrace_posterior_g(group_indx,:,:)))', 'k')
set(h(1), 'LineWidth', 1)

hold on
h(2) = lineplot_error(t,(squeeze(summary.ko.Mtrace_posterior_g(group_indx,:,:)))', 'r')

set(h(2), 'LineWidth', 1)
xlim([-0.1, 0.5])
xlabel('Time(s)')
ylabel('\Delta Normalize Movement')
%% 
figure;
step = 0.033;
t= -0.033 * 9: step: -0.033; % the event is between 10 and 11 th frame
t_post = 0:step:17*step;
t = [t, t_post];

indx_bin = find(t > 0.1 & t<0.3);
% max(video_sum.wt(1).video_summary.Mtrace_Posterior_norm(:, indx_bin), [],2)
% level = 35:10:95;
genotype = {'ko', 'wt'}
for kk = 1:length(genotype)
        summary_Mtrace_anterior.(genotype{kk}).peak = squeeze(max(summary.(genotype{kk}).Mtrace_anterior_g(:, indx_bin,:),[],2));
        summary_Mtrace_posterior.(genotype{kk}).peak = squeeze(max(summary.(genotype{kk}).Mtrace_posterior_g(:, indx_bin,:),[],2));
end

figure;
line_errorbar_drc((summary_Mtrace_posterior.wt.peak./summary_Mtrace_posterior.wt.peak(1,:))', (summary_Mtrace_posterior.ko.peak./summary_Mtrace_posterior.ko.peak(1,:))')
% line_errorbar_drc(summary_Mtrace_posterior.wt.peak', summary_Mtrace_posterior.ko.peak')
xticks([0, 1:10, 11])
xlim([0,11])
xticklabels({'','1', '2', '3', '4','5', '6','7','8', '9', '10', ''})
ylabel('\Delta Normalize Face Movement')
% ylabel('\Delta Normalize Whiskerpad Movement')
xlabel('Every 10 trials')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,600])



figure;
line_errorbar_drc((summary_Mtrace_anterior.wt.peak./summary_Mtrace_anterior.wt.peak(1,:))', (summary_Mtrace_anterior.ko.peak./summary_Mtrace_anterior.ko.peak(1,:))')
% line_errorbar_drc(summary_Mtrace_anterior.wt.peak', summary_Mtrace_anterior.ko.peak')
xticks([0, 1:10, 11])
xlim([0,11])
xticklabels({'','1', '2', '3', '4','5', '6','7','8', '9', '10', ''})
% ylabel('\Delta Normalize Face Movement')
ylabel('\Delta Normalize Whiskerpad Movement')
xlabel('Every 10 trials')
set(gca,'TickDir','out')
set(gca,'fontsize',12)
set(gca,'TickLengt', [0.015 0.015]);
set(gca, 'LineWidth',1)
set(gcf,'position',[100,200,600,600])